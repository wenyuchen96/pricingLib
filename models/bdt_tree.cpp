#include <iostream>
#include <string>
#include <vector>
#include <cmath>

struct BondOptionResult
{
    double call;
    double put;
};

enum class OptionType
{
    EUROPEAN = 1,
    BERMUDAN = 2,
    AMERICAN = 3
};

BondOptionResult americanBondOption(
    double t_exp, double t_mat, double strike,
    double faceAmount,
    const std::vector<double> &cpnTimes,
    const std::vector<double> &cpnFlows,
    OptionType exerciseType,
    const std::vector<double> &dfTimes,
    const std::vector<double> &dfValues,
    const std::vector<double> &treeTimes,
    const std::vector<std::vector<double>> &Q,
    const std::vector<std::vector<double>> &rt,
    double dt)
{
    int nLevels = treeTimes.size();
    int expIndex = -1, matIndex = -1;

    for (int i = 0; i < nLevels; ++i)
    {
        if (std::abs(treeTimes[i] - t_exp) < 1e-6)
            expIndex = i;
        if (std::abs(treeTimes[i] - t_mat) < 1e-6)
            matIndex = i;
    }
    if (expIndex == -1 || matIndex == -1)
        throw std::runtime_error("Could not align tree time steps with expiry or maturity.");

    // Lambda to run backward induction with early exercise logic
    auto buildValueTree = [&](bool isCallable) -> double
    {
        std::vector<std::vector<double>> valueTree(nLevels);
        for (int i = 0; i < nLevels; ++i)
            valueTree[i].resize(i + 1, 0.0);

        // Add coupons
        for (size_t k = 0; k < cpnTimes.size(); ++k)
        {
            double t = cpnTimes[k];
            double cpn = cpnFlows[k];
            for (int i = 0; i < nLevels; ++i)
            {
                if (std::abs(treeTimes[i] - t) < 1e-6)
                {
                    for (int j = 0; j <= i; ++j)
                        valueTree[i][j] += cpn;
                }
            }
        }

        // Add face value
        for (int j = 0; j <= matIndex; ++j)
            valueTree[matIndex][j] += faceAmount;

        // Backward induction with optional early exercise
        for (int i = matIndex - 1; i >= 0; --i)
        {
            for (int j = 0; j <= i; ++j)
            {
                double contVal = 0.5 * (valueTree[i + 1][j] + valueTree[i + 1][j + 1]);
                contVal *= std::exp(-rt[i][j] * dt);

                double nodeVal = contVal;

                if (treeTimes[i] <= t_exp)
                {
                    if (exerciseType == OptionType::AMERICAN ||
                        (exerciseType == OptionType::EUROPEAN && std::abs(treeTimes[i] - t_exp) < 1e-6))
                    {
                        if (isCallable)
                        {
                            nodeVal = std::min(strike, contVal); // issuer redeems at strike
                        }
                        else
                        {
                            nodeVal = std::max(strike, contVal); // holder sells at strike
                        }
                    }
                }

                valueTree[i][j] = nodeVal;
            }
        }

        return valueTree[0][0];
    };

    double callValue = buildValueTree(true); // callable bond
    double putValue = buildValueTree(false); // putable bond

    return BondOptionResult{callValue, putValue};
}

// Short-rate process under the Black-Derman-Toy (BDT) model:
// The logarithm of the short rate follows a stochastic process:
// d(log r) = theta(t) * dt + sigma * dW_t
class BDTTree
{
public:
    BDTTree(double sigma, int numTimeSteps, const std::vector<double> &marketTimes, const std::vector<double> &marketDfs) : sigma_{sigma},
                                                                                                                            numTimeSteps_{numTimeSteps},
                                                                                                                            marketTimes_{marketTimes},
                                                                                                                            marketDfs_{marketDfs}
    {
        if (sigma_ < 0.0)
        {
            throw std::invalid_argument("Sigma must be non-negative");
        }
        if (numTimeSteps_ < 3)
        {
            throw std::invalid_argument("Must have at least 3 time steps");
        }
    }
    // helper function to calculate ZCB prices for calibration
    double computeZCBPrice(int m, double alpha, double dt, double sigma, const std::vector<std::vector<double>> &Q)
    {
        double sum = 0.0;
        for (int j = 0; j <= m; ++j)
        {
            double r = alpha * std::exp(2.0 * j * sigma * std::sqrt(dt));
            sum += Q[m][j] * std::exp(-r * dt);
        }
        return sum;
    }
    // calibrate the alpha parameter in the short rate lattice
    double calibrateAlpha(int m, double marketDf, double dt, double sigma, const std::vector<std::vector<double>> &Q)
    {
        double low{1e-6}, high{1.0}, tol{1e-8};
        double mid;

        for (int iter = 0; iter < 100; ++iter)
        {
            mid = 0.5 * (low + high);
            double price = computeZCBPrice(m, mid, dt, sigma, Q);

            if (std::abs(price - marketDf) < tol)
                return mid;
            if (price > marketDf)
                low = mid;
            else
                high = mid;
        }
        return mid;
    }
    // build the calibrated short rate lattice
    void buildTree()
    {
        double treeMaturity = marketTimes_.back();   // last time point from market
        dt_ = treeMaturity / (numTimeSteps_ + 1);    // +1 for step at expiry
        treeTimes_.resize(numTimeSteps_ + 2);        // +1 for t = 0 and +1 for maturity
        for (int i = 0; i <= numTimeSteps_ + 1; ++i) // initialize time values for tree
        {
            treeTimes_[i] = i * dt_;
        }

        marketDfTree_.resize(treeTimes_.size()); // resize time values for market tree
        // interpolate the observed discount factors to match the resized market tree
        for (size_t i = 0; i < treeTimes_.size(); ++i)
        {
            double t = treeTimes_[i];
            // linear interpolation
            if (t <= marketTimes_.front())
            {
                marketDfTree_[i] = marketDfs_.front();
            }
            else if (t >= marketTimes_.back())
            {
                marketDfTree_[i] = marketDfs_.back();
            }
            else
            {
                for (size_t j = 0; j < marketTimes_.size() - 1; ++j)
                {
                    if (marketTimes_[j] <= t && t < marketTimes_[j + 1])
                    {
                        double t0 = marketTimes_[j], t1 = marketTimes_[j + 1];
                        double df0 = marketDfs_[j], df1 = marketDfs_[j + 1];
                        marketDfTree_[i] = df0 + (df1 - df0) * (t - t0) / (t1 - t0);
                        break;
                    }
                }
            }
        }
        // resize Q_ and rt_
        int levels = numTimeSteps_ + 2; // t=0 to t=T
        Q_.resize(levels);
        rt_.resize(levels);
        for (int i = 0; i < levels; ++i)
        {
            Q_[i].resize(i + 1, 0.0);
            rt_[i].resize(i + 1, 0.0);
        }
        // set initial short rate rt_ and elementary prices Q_
        double r0 = -std::log(marketDfTree_[1]) / dt_;
        rt_[0][0] = r0;
        Q_[0][0] = 1.0;
        Q_[1][0] = 0.5 * std::exp(-rt_[0][0] * dt_);
        Q_[1][1] = 0.5 * std::exp(-rt_[0][0] * dt_);

        // start the calibration
        for (int m = 1; m <= numTimeSteps_; ++m)
        {
            double alpha = calibrateAlpha(m, marketDfTree_[m + 1], dt_, sigma_, Q_);
            for (int j = 0; j <= m; ++j)
            {
                rt_[m][j] = alpha * std::exp(2.0 * j * sigma_ * std::sqrt(dt_));
            }
            // update Q_ lattice based on the calibrated rt_
            Q_[m + 1][0] = 0.5 * Q_[m][0] * std::exp(-rt_[m][0] * dt_);
            for (int j = 1; j <= m; ++j)
            {
                Q_[m + 1][j] = 0.5 * Q_[m][j - 1] * std::exp(-rt_[m][j - 1] * dt_) + 0.5 * Q_[m][j] * std::exp(-rt_[m][j] * dt_);
            }
            Q_[m + 1][m + 1] = 0.5 * Q_[m][m] * std::exp(-rt_[m][m] * dt_);
        }
    }
    // value american callable and putable bonds
    BondOptionResult bondOption(
        double t_exp, double strikePrice, double faceAmount,
        const std::vector<double> &cpnTimes,
        const std::vector<double> &cpnFlows,
        OptionType exerciseType)
    {
        double bondMat = cpnTimes.back(); // t_exp is option expiry and bondMat is bond maturity

        if (t_exp > bondMat)
        {
            throw std::invalid_argument("Option expiry after bond matures.");
        }
        if (t_exp < 0.0)
        {
            throw std::invalid_argument("Option expiry time negative.");
        }

        return americanBondOption(
            t_exp, bondMat, strikePrice, faceAmount,
            cpnTimes, cpnFlows, exerciseType,
            marketTimes_, marketDfs_, treeTimes_, Q_, rt_, dt_);
    }

private:
    double sigma_;                        // Constant volatility of the short rate process
    int numTimeSteps_;                    // Number of time steps in the binomial tree (exlcude t = 0 and t = maturity)
    std::vector<std::vector<double>> Q_;  // Risk-neutral probability tree (discounted probabilities)
    std::vector<std::vector<double>> rt_; // Short rates at each node of the tree
    std::vector<double> treeTimes_;       // Time values corresponding to each tree level (0, dt, 2dt, ..., T)
    std::vector<double> marketTimes_;     // Input times for discount factor curve
    std::vector<double> marketDfs_;       // Input discount factor values (e.g., from market zero curve)
    std::vector<double> marketDfTree_;    // Interpolated discount factors at tree times
    double uProbs{0.5};                   // Probability of upward movement in the binomial tree (default 0.5)
    double dProbs{0.5};                   // Probability of downward movement in the binomial tree (default 0.5)
    double dt_;                           // Time step size = tree maturity / (numTimeSteps + 1)
};