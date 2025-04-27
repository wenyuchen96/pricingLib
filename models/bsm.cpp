#include <iostream>
#include <string>
#include <cmath>
#include <boost/math/distributions/normal.hpp>

// implement Black - Scholes - Merton ("BSM") Option Pricing Model to value a european call and a put option

/*
BSM model inputs
s: stock price
k: strike price
t: time to expiration
r: risk-free rate
q: dividend yield
sigma: volatility
type: option type "call" or "put"
*/

class BlackScholes
{
public:
    BlackScholes(double s, double k, double t, double r, double q, double sigma, std::string type) : s_(s), k_(k), t_(t), r_(r), q_(q), sigma_(sigma), type_(type) {};

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

    std::pair<double, double> d1_d2() const
    {
        double d1 = (log(s_ / k_) + (r_ - q_ + pow(sigma_, 2) / 2.0) * t_) / (sigma_ * sqrt(t_));
        double d2 = d1 - sigma_ * sqrt(t_);
        return std::make_pair(d1, d2);
    }

    // black scholes formula to calculate european call or put option value
    double calculate()
    {
        // define the variables
        auto [d1, d2] = d1_d2();
        double c;
        double p;

        if (type_ == "call")
        {
            c = s_ * exp(-q_ * t_) * NormalCdf(d1) - k_ * exp(-r_ * t_) * NormalCdf(d2);
            return c;
        }
        else if (type_ == "put")
        {
            p = k_ * exp(-r_ * t_) * NormalCdf(-d2) - s_ * exp(-q_ * t_) * NormalCdf(-d1);
            return p;
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

    // greek letters
    double delta()
    {
        auto [d1, d2] = d1_d2();
        if (type_ == "call")
        {
            return exp(-q_ * t_) * NormalCdf(d1);
        }
        else if (type_ == "put")
        {
            return -exp(-q_ * t_) * NormalCdf(-d1);
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

    double gamma()
    {
        auto [d1, d2] = d1_d2();
        return (exp(-q_ * t_) * NormalPdf(d1)) / (s_ * sigma_ * sqrt(t_));
    }

    double vega()
    {
        auto [d1, d2] = d1_d2();
        return s_ * exp(-q_ * t_) * NormalPdf(d1) * sqrt(t_);
    }

    double theta()
    {
        auto [d1, d2] = d1_d2();
        double term1 = -(s_ * exp(-q_ * t_) * NormalPdf(d1) * sigma_) / (2 * sqrt(t_));
        if (type_ == "call")
        {
            return term1 - r_ * k_ * exp(-r_ * t_) * NormalCdf(d2) + q_ * s_ * exp(-q_ * t_) * NormalCdf(d1);
        }
        else if (type_ == "put")
        {
            return term1 + r_ * k_ * exp(-r_ * t_) * NormalCdf(-d2) - q_ * s_ * exp(-q_ * t_) * NormalCdf(-d1);
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

    double rho()
    {
        auto [d1, d2] = d1_d2();
        if (type_ == "call")
            return k_ * t_ * exp(-r_ * t_) * NormalCdf(d2);
        else if (type_ == "put")
            return -k_ * t_ * exp(-r_ * t_) * NormalCdf(-d2);
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

private:
    double s_;
    double k_;
    double t_;
    double r_;
    double q_;
    double sigma_;
    std::string type_;
};

int main()
{

    // BSM model inputs
    double s = 10.5;
    double k = 5.5;
    double t = 2.0;
    double r = 0.045;
    double q = 0.0;
    double sigma = 0.6;
    std::string type = "call";

    // Create an instance of the BlackScholes class
    BlackScholes bs(s, k, t, r, q, sigma, type);

    std::cout << type << " option price is $" << bs.calculate() << "\n";
    std::cout << "delta: " << bs.delta() << "\n";
    std::cout << "gamma: " << bs.gamma() << "\n";
    std::cout << "vega: " << bs.vega() << "\n";
    std::cout << "theta: " << bs.theta() << "\n";
    std::cout << "rho: " << bs.rho() << "\n";
}