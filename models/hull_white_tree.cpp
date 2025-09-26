#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp>

// implement hull-white 1 factor model
// dr = (theta(t) - r) dt + sigma * dW

enum class EuropeanCalcType
{
    EXPIRY_TREE,
};

struct ZCBOptionPrice
{
    double call;
    double put;
};

class HullWhiteTree
{
public:
    // helper functions
    //  cumulative density function of a standard normal distribution
    static double NormalCdf(double x)
    {
        boost::math::normal dist(0.0, 1.0);
        return boost::math::cdf(dist, x);
    }
    // probability density function of a standard normal distribution
    static double NormalPdf(double x)
    {
        boost::math::normal dist(0.0, 1.0);
        return boost::math::pdf(dist, x);
    }

    // linear interpolation
    static double interpolate(double t, const std::vector<double> &times, const std::vector<double> &values)
    {
        if (t <= times.front())
            return values.front();
        if (t >= times.back())
            return values.back();

        for (size_t i = 0; i < times.size() - 1; ++i)
        {
            if (t >= times[i] && t <= times[i + 1])
            {
                double w = (t - times[i]) / (times[i + 1] - times[i]);
                return values[i] * (1 - w) + values[i + 1] * w;
            }
        }
        throw std::runtime_error("Interpolation failed.");
    }

    HullWhiteTree(double sigma, double a, int numTimeSteps, const std::vector<double> &treeTimes) : sigma_{sigma}, a_{a}, numTimesSteps_{numTimeSteps}, treeTimes_{treeTimes}, treeBuilt_{false} {}

    // analytical solution of option on a zero coupon bond
    ZCBOptionPrice optionOnZCB(double tExp, double tMat, double strike, double faceAmount, const std::vector<double> &dfTimes, const std::vector<double> &dfValues) const
    {
        if (tExp > tMat)
        {
            throw std::invalid_argument("Option expiry is after bond maturity.");
        }
        if (tExp < 0.0)
        {
            throw std::invalid_argument("Option expiry time cannot be negative");
        }

        double ptExp{interpolate(tExp, dfTimes, dfValues)};
        double ptMat{interpolate(tMat, dfTimes, dfValues)};

        double sigma{sigma_};
        double a{std::max(a_, 1e-6)};

        double sigmap{(sigma / a) * (1.0 - std::exp(-a * (tMat - tExp)))};
        sigmap *= std::sqrt((1.0 - std::exp(-2.0 * a * tExp)) / (2.0 * a));
        sigmap = std::max(sigmap, 1e-6);

        double h{std::log((faceAmount * ptMat) / (strike * ptExp)) / sigmap + 0.5 * sigmap};

        double call{faceAmount * ptMat * NormalCdf(h) - strike * ptExp * NormalCdf(h - sigmap)};
        double put{strike * ptExp * NormalCdf(-h + sigmap) - faceAmount * ptMat * NormalCdf(-h)};
        return ZCBOptionPrice{call, put};
    }

    void buildTrinomialTree(const std::vector<double> &dfTimes, const std::vector<double> &dfValues)
    {
        // save market input curves
        dfTimes_ = dfTimes;
        dfs_ = dfValues;

        // extend tree maturity slightly beyond given tree maturity
        double treeMaturity{treeTimes_.back() * (numTimesSteps_ + 1.0) / numTimesSteps_};
        treeTimes_.resize(numTimesSteps_ + 2);

        for (int i = 0; i <= numTimesSteps_ + 1; ++i)
        {
            treeTimes_[i] = i * treeMaturity / (numTimesSteps_ + 1.0);
        }
        // dfTree is the interpolated discount factors match the treeTimes
        std::vector<double> dfTree(numTimesSteps_ + 2);
        dfTree[0] = 1.0;
        for (int i = 1; i <= numTimesSteps_ + 1; ++i)
        {
            dfTree[i] = interpolate(treeTimes_[i], dfTimes, dfValues);
        }

        // start building tree
        const double dt{treeMaturity / (numTimesSteps_ + 1)};
        dt_ = dt;
        const double dR{sigma_ * std::sqrt(3.0 * dt)};
        const int jMax{static_cast<int>(std::ceil(0.1835 / (a_ * dt)))};
        const int N{jMax};

        Q_.assign(numTimesSteps_ + 2, std::vector<double>(2 * N + 1, 0.0));
        pu_.assign(2 * N + 1, 0.0);
        pm_.assign(2 * N + 1, 0.0);
        pd_.assign(2 * N + 1, 0.0);
        r_.assign(numTimesSteps_ + 2, std::vector<double>(2 * N + 1, 0.0));
        r_t_.assign(numTimesSteps_ + 1, 0.0);

        Q_[0][N] = 1.0;

        // fill in the truncated probabilities
        for (int j = -jMax; j <= jMax; ++j)
        {
            const int jN{j + N};
            double ajdt{a_ * j * dt};

            if (j == jMax)
            {
                pu_[jN] = 7.0 / 6.0 + 0.5 * (ajdt * ajdt - 3.0 * ajdt);
                pm_[jN] = -1.0 / 3.0 - ajdt * ajdt + 2.0 * ajdt;
                pd_[jN] = 1.0 / 6.0 + 0.5 * (ajdt * ajdt - ajdt);
            }
            else if (j == -jMax)
            {
                pu_[jN] = 1.0 / 6.0 + 0.5 * (ajdt * ajdt + ajdt);
                pm_[jN] = -1.0 / 3.0 - ajdt * ajdt - 2.0 * ajdt;
                pd_[jN] = 7.0 / 6.0 + 0.5 * (ajdt * ajdt + 3.0 * ajdt);
            }
            else
            {
                pu_[jN] = 1.0 / 6.0 + 0.5 * (ajdt * ajdt - ajdt);
                pm_[jN] = 2.0 / 3.0 - ajdt * ajdt;
                pd_[jN] = 1.0 / 6.0 + 0.5 * (ajdt * ajdt + ajdt);
            }
        }

        // build tree
        for (int m = 0; m <= numTimesSteps_; ++m)
        {
            int nm = std::min(m, jMax);
            double sum_qz = 0.0;

            for (int j = -nm; j <= nm; ++j)
            {
                int jN = j + N;
                double rdt = j * dR * dt;
                sum_qz += Q_[m][jN] * std::exp(-rdt);
            }

            r_t_[m] = std::log(sum_qz / dfTree[m + 1] / dt);

            for (int j = -nm; j <= nm; ++j)
            {
                int jN = j + N;
                r_[m][jN] = r_t_[m] + j * dR;
                double rdt = r_[m][jN] * dt;
                double z = std::exp(-rdt);

                if (j == jMax)
                {
                    Q_[m + 1][jN] += Q_[m][jN] * pu_[jN] * z;
                    Q_[m + 1][jN - 1] += Q_[m][jN] * pm_[jN] * z;
                    Q_[m + 1][jN - 2] += Q_[m][jN] * pd_[jN] * z;
                }
                else if (j == -jMax)
                {
                    Q_[m + 1][jN] += Q_[m][jN] * pd_[jN] * z;
                    Q_[m + 1][jN + 1] += Q_[m][jN] * pm_[jN] * z;
                    Q_[m + 1][jN + 2] += Q_[m][jN] * pu_[jN] * z;
                }
                else
                {
                    Q_[m + 1][jN + 1] += Q_[m][jN] * pu_[jN] * z;
                    Q_[m + 1][jN] += Q_[m][jN] * pm_[jN] * z;
                    Q_[m + 1][jN - 1] += Q_[m][jN] * pd_[jN] * z;
                }
            }
        }

        treeBuilt_ = true;
    }

private:
    double sigma_;                      // volatility
    double a_;                          // speed of mean reversion
    int numTimesSteps_;                 // number of time steps in the tree, not include t = 0 and at maturity
    EuropeanCalcType EuropeanCalcType_; // calculation type for european option

    std::vector<std::vector<double>> Q_; // value matrix including option or bond values
    std::vector<std::vector<double>> r_; // short rate tree, r_[i][j] is the short rate at time i, state j
    std::vector<double> treeTimes_;      // time steps

    // trinomial trees
    std::vector<double> pu_; // upward probability
    std::vector<double> pm_; // middle probability
    std::vector<double> pd_; // downward probability

    // discount curve
    std::vector<double> dfTimes_; // time steps for discount factors inputs
    std::vector<double> dfs_;     // corresponding discount factors

    std::vector<double> r_t_; // theta(t) or draft adjustments

    bool treeBuilt_; // flag if tree is built
    double dt_;      // time step size = totalTime / numTimeSteps
};