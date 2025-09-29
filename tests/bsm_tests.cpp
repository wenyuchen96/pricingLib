#include <catch2/catch_all.hpp>

#include "../models/bsm.cpp"

namespace
{
    struct Inputs
    {
        double s;
        double k;
        double t;
        double r;
        double q;
        double sigma;
    };

    std::pair<double, double> analytic_d1d2(const Inputs &in)
    {
        const double denom = in.sigma * std::sqrt(in.t);
        const double d1 = (std::log(in.s / in.k) + (in.r - in.q + 0.5 * in.sigma * in.sigma) * in.t) / denom;
        const double d2 = d1 - denom;
        return {d1, d2};
    }

    double analytic_price(const Inputs &in, const std::string &type)
    {
        const auto [d1, d2] = analytic_d1d2(in);
        if (type == "call")
        {
            return in.s * std::exp(-in.q * in.t) * BlackScholes::NormalCdf(d1) -
                   in.k * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(d2);
        }
        else
        {
            return in.k * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(-d2) -
                   in.s * std::exp(-in.q * in.t) * BlackScholes::NormalCdf(-d1);
        }
    }

    double analytic_delta(const Inputs &in, const std::string &type)
    {
        const auto [d1, _] = analytic_d1d2(in);
        if (type == "call")
            return std::exp(-in.q * in.t) * BlackScholes::NormalCdf(d1);
        return -std::exp(-in.q * in.t) * BlackScholes::NormalCdf(-d1);
    }

    double analytic_gamma(const Inputs &in)
    {
        const auto [d1, _] = analytic_d1d2(in);
        return std::exp(-in.q * in.t) * BlackScholes::NormalPdf(d1) / (in.s * in.sigma * std::sqrt(in.t));
    }

    double analytic_vega(const Inputs &in)
    {
        const auto [d1, _] = analytic_d1d2(in);
        return in.s * std::exp(-in.q * in.t) * BlackScholes::NormalPdf(d1) * std::sqrt(in.t);
    }

    double analytic_theta(const Inputs &in, const std::string &type)
    {
        const auto [d1, d2] = analytic_d1d2(in);
        const double term1 = -(in.s * std::exp(-in.q * in.t) * BlackScholes::NormalPdf(d1) * in.sigma) /
                              (2.0 * std::sqrt(in.t));
        if (type == "call")
        {
            return term1 - in.r * in.k * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(d2) +
                   in.q * in.s * std::exp(-in.q * in.t) * BlackScholes::NormalCdf(d1);
        }
        else
        {
            return term1 + in.r * in.k * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(-d2) -
                   in.q * in.s * std::exp(-in.q * in.t) * BlackScholes::NormalCdf(-d1);
        }
    }

    double analytic_rho(const Inputs &in, const std::string &type)
    {
        const auto [_, d2] = analytic_d1d2(in);
        if (type == "call")
            return in.k * in.t * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(d2);
        return -in.k * in.t * std::exp(-in.r * in.t) * BlackScholes::NormalCdf(-d2);
    }
}

TEST_CASE("Black-Scholes price and Greeks match analytic formulas")
{
    const Inputs params{100.0, 95.0, 1.25, 0.03, 0.01, 0.24};

    BlackScholes callPricer(params.s, params.k, params.t, params.r, params.q, params.sigma, "call");
    BlackScholes putPricer(params.s, params.k, params.t, params.r, params.q, params.sigma, "put");

    const double expected_call = analytic_price(params, "call");
    const double expected_put = analytic_price(params, "put");

    CHECK(callPricer.calculate() == Catch::Approx(expected_call).margin(1e-12));
    CHECK(putPricer.calculate() == Catch::Approx(expected_put).margin(1e-12));

    CHECK(callPricer.delta() == Catch::Approx(analytic_delta(params, "call")).margin(1e-12));
    CHECK(putPricer.delta() == Catch::Approx(analytic_delta(params, "put")).margin(1e-12));

    CHECK(callPricer.gamma() == Catch::Approx(analytic_gamma(params)).margin(1e-12));
    CHECK(callPricer.vega() == Catch::Approx(analytic_vega(params)).margin(1e-12));

    CHECK(callPricer.theta() == Catch::Approx(analytic_theta(params, "call")).margin(1e-12));
    CHECK(putPricer.theta() == Catch::Approx(analytic_theta(params, "put")).margin(1e-12));

    CHECK(callPricer.rho() == Catch::Approx(analytic_rho(params, "call")).margin(1e-12));
    CHECK(putPricer.rho() == Catch::Approx(analytic_rho(params, "put")).margin(1e-12));
}

TEST_CASE("Black-Scholes throws on invalid type")
{
    BlackScholes pricer(100.0, 100.0, 1.0, 0.05, 0.0, 0.2, "strangle");
    CHECK_THROWS_AS(pricer.calculate(), std::invalid_argument);
    CHECK_THROWS_AS(pricer.delta(), std::invalid_argument);
    CHECK_THROWS_AS(pricer.theta(), std::invalid_argument);
    CHECK_THROWS_AS(pricer.rho(), std::invalid_argument);
}
