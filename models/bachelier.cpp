#include <iostream>
#include <string>
#include <cmath>
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

int main()
{
    // BSM inputs
    double s = 10.5;
    double k = 5.5;
    double t = 2.0;
    double r = 0.045;
    double q = 0.0;
    double sigma = 0.6;
    std::string type = "call";

    // Convert to Bachelier-compatible inputs
    double f = s * std::exp((r - q) * t); // Forward price
    double df = std::exp(-r * t);         // Discount factor
    double sigma_bachelier = sigma * s;   // Approx absolute vol
    // Create Bachelier instance
    Bachelier bachelier(f, k, t, df, sigma_bachelier, type);

    std::cout << "Bachelier " << type << " option price is $" << bachelier.calculate() << "\n";
}