#include <catch2/catch_all.hpp>

#include "../models/black_shifted.cpp"

namespace
{
    double shifted_call(double f, double shift, double k, double t, double df, double sigma)
    {
        const double d1 = (std::log(f / k) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t));
        const double d2 = d1 - sigma * std::sqrt(t);
        return df * ((f + shift) * BlackShifted::NormalCdf(d1) - (k + shift) * BlackShifted::NormalCdf(d2));
    }

    double shifted_put(double f, double shift, double k, double t, double df, double sigma)
    {
        const double d1 = (std::log(f / k) + 0.5 * sigma * sigma * t) / (sigma * std::sqrt(t));
        const double d2 = d1 - sigma * std::sqrt(t);
        return df * ((k + shift) * BlackShifted::NormalCdf(-d2) - (f + shift) * BlackShifted::NormalCdf(-d1));
    }
}

TEST_CASE("Shifted Black reproduces unshifted price when shift is zero")
{
    const double f = 50.0;
    const double k = 45.0;
    const double t = 0.5;
    const double df = std::exp(-0.02 * t);
    const double sigma = 0.35;

    BlackShifted pricer(f, 0.0, k, t, df, sigma, "call");
    const double shiftedPrice = pricer.calculate();

    const double baseline = shifted_call(f, 0.0, k, t, df, sigma);
    CHECK(shiftedPrice == Catch::Approx(baseline).margin(1e-12));
}

TEST_CASE("Shifted Black handles non-zero shift consistently")
{
    const double f = 0.015;
    const double shift = -0.01;
    const double k = 0.01;
    const double t = 1.5;
    const double df = std::exp(-0.015 * t);
    const double sigma = 0.01;

    BlackShifted callPricer(f, shift, k, t, df, sigma, "call");
    BlackShifted putPricer(f, shift, k, t, df, sigma, "put");

    const double call = callPricer.calculate();
    const double put = putPricer.calculate();

    const double expectedCall = shifted_call(f, shift, k, t, df, sigma);
    const double expectedPut = shifted_put(f, shift, k, t, df, sigma);

    CHECK(call == Catch::Approx(expectedCall).margin(1e-12));
    CHECK(put == Catch::Approx(expectedPut).margin(1e-12));
}

TEST_CASE("Shifted Black rejects invalid type")
{
    BlackShifted pricer(0.02, 0.03, -0.01, 1.0, 0.98, 0.01, "digital");
    CHECK_THROWS_AS(pricer.calculate(), std::invalid_argument);
}
