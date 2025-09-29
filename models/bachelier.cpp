#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/math/distributions/normal.hpp>

// implement Bachelier's Model to value a european call and a put option

/*
Bachelier model inputs
f: forward price (F)
k: strike price (K)
t: time to expiration (in years)
df: discount factor (e.g., exp(-r * t))
sigma: absolute volatility (in price units)
type: option type "call" or "put"
*/
class Bachelier
{
public:
    Bachelier(double f, double k, double t, double df, double sigma, std::string type) : f_{f}, k_{k}, t_{t}, df_{df}, sigma_{sigma}, type_{type} {}

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
        double root_t = std::sqrt(t_);
        double d = (f_ - k_) / (sigma_ * root_t);

        if (type_ == "call")
        {
            return df_ * ((f_ - k_) * NormalCdf(d) + sigma_ * root_t * NormalPdf(d));
        }
        else if (type_ == "put")
        {
            return df_ * ((k_ - f_) * NormalCdf(-d) + sigma_ * root_t * NormalPdf(d));
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

private:
    double f_;
    double k_;
    double t_;
    double df_;
    double sigma_;
    std::string type_;
};
