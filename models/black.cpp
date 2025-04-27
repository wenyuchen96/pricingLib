#include <iostream>
#include <string>
#include <cmath>
#include <boost/math/distributions/normal.hpp>

// Black’s model prices European options on forwards or futures by discounting the expected payoff under the forward measure using a known discount factor.
class Black
{
public:
    Black(double f, double k, double t, double df, double sigma, std::string type) : f_{f}, k_{k}, t_{t}, df_{df}, sigma_{sigma}, type_{type} {};

    // cumulative density function of a standard normal distribution
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

    double calculate()
    {
        // pricing a european derivative using Black model
        auto [d1, d2] = d1_d2();
        if (type_ == "call")
        {
            return df_ * (f_ * NormalCdf(d1) - k_ * NormalCdf(d2));
        }
        else if (type_ == "put")
        {
            return df_ * (k_ * NormalCdf(-d2) - f_ * NormalCdf(-d1));
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

    std::pair<double, double> d1_d2() const
    {
        double d1 = (log(f_ / k_) + sigma_ * sigma_ * t_ / 2.0) / (sigma_ * sqrt(t_));
        double d2 = d1 - sigma_ * sqrt(t_);
        return std::make_pair(d1, d2);
    }

private:
    double f_;         // forward rate f
    double k_;         // strike price k
    double t_;         // time to expiry (years)
    double df_;        // df rfr to expiry date
    double sigma_;     // volatility
    std::string type_; // call or put
};