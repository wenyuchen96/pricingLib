#include <catch2/catch_all.hpp>

#include "../models/black.cpp"

namespace
{
    double call_price_closed_form(double f, double k, double t, double df, double sigma)
    {
        const double d1 = (std::log(f / k) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t));
        const double d2 = d1 - sigma * std::sqrt(t);
        return df * (f * Black::NormalCdf(d1) - k * Black::NormalCdf(d2));
    }

    double put_price_closed_form(double f, double k, double t, double df, double sigma)
    {
        const double d1 = (std::log(f / k) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t));
        const double d2 = d1 - sigma * std::sqrt(t);
        return df * (k * Black::NormalCdf(-d2) - f * Black::NormalCdf(-d1));
    }
}

TEST_CASE("Black call and put recover analytic price")
{
    const double f = 105.0;
    const double k = 100.0;
    const double t = 1.25;
    const double df = std::exp(-0.03 * t);
    const double sigma = 0.25;

    const double expected_call = call_price_closed_form(f, k, t, df, sigma);
    const double expected_put = put_price_closed_form(f, k, t, df, sigma);

    REQUIRE(expected_call > 0.0);
    REQUIRE(expected_put > 0.0);

    Black callPricer(f, k, t, df, sigma, "call");
    Black putPricer(f, k, t, df, sigma, "put");

    CHECK(callPricer.calculate() == Catch::Approx(expected_call).margin(1e-12));
    CHECK(putPricer.calculate() == Catch::Approx(expected_put).margin(1e-12));
}

TEST_CASE("Black d1 and d2 maintain no-arbitrage relationship")
{
    const double f = 98.0;
    const double k = 101.0;
    const double t = 0.75;
    const double df = std::exp(-0.01 * t);
    const double sigma = 0.2;

    Black pricer(f, k, t, df, sigma, "call");
    const auto [d1, d2] = pricer.d1_d2();

    CHECK(d2 == Catch::Approx(d1 - sigma * std::sqrt(t)).margin(1e-12));
}

TEST_CASE("Black throws on unsupported option type")
{
    Black pricer(100.0, 100.0, 1.0, 0.95, 0.3, "straddle");
    CHECK_THROWS_AS(pricer.calculate(), std::invalid_argument);
}
