#include <catch2/catch_all.hpp>

#include "../models/bachelier.cpp"

TEST_CASE("Bachelier pricing matches benchmark values for spot inputs")
{
    const double s = 10.5;
    const double k = 5.5;
    const double t = 2.0;
    const double r = 0.045;
    const double q = 0.0;
    const double sigma = 0.6;

    const double f = s * std::exp((r - q) * t);
    const double df = std::exp(-r * t);
    const double sigma_bachelier = sigma * s;

    auto price = [&](const std::string &type)
    {
        Bachelier pricer(f, k, t, df, sigma_bachelier, type);
        return pricer.calculate();
    };

    SECTION("Call price matches closed-form benchmark")
    {
        CHECK(price("call") == Catch::Approx(6.6926036883).margin(1e-9));
    }

    SECTION("Put-call parity holds under Bachelier assumptions")
    {
        const double call = price("call");
        const double put = price("put");
        const double parity = df * (f - k);
        CHECK(call - put == Catch::Approx(parity).margin(1e-12));
    }

    SECTION("Unknown option type triggers an exception")
    {
        CHECK_THROWS_AS(price("straddle"), std::invalid_argument);
    }
}
