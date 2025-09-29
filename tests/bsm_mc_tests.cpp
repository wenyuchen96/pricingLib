#include <catch2/catch_all.hpp>

#include "../models/bsm.cpp"
#include "../models/bsm_mc.cpp"

namespace
{
    struct MCInputs
    {
        double s;
        double k;
        double t;
        double r;
        double q;
        double sigma;
        std::string type;
        int paths;
        int seed;
    };
}

TEST_CASE("Monte Carlo engines agree for deterministic seed")
{
    const MCInputs inputs{100.0, 95.0, 1.0, 0.02, 0.0, 0.3, "call", 200'000, 99};

    const double loopPrice = euro_mc_for_loop(inputs.s, inputs.k, inputs.t, inputs.r, inputs.q,
                                              inputs.sigma, inputs.type, inputs.paths, inputs.seed);
    const double valarrayPrice = euro_mc_valarray(inputs.s, inputs.k, inputs.t, inputs.r, inputs.q,
                                                  inputs.sigma, inputs.type, inputs.paths, inputs.seed);
    const double simdPrice = euro_mc_xsimd(inputs.s, inputs.k, inputs.t, inputs.r, inputs.q,
                                           inputs.sigma, inputs.type, inputs.paths, inputs.seed);

    CHECK(valarrayPrice == Catch::Approx(loopPrice).margin(1e-12));
    CHECK(simdPrice == Catch::Approx(loopPrice).margin(1e-9));
}

TEST_CASE("Monte Carlo results align with analytic Black-Scholes price")
{
    const MCInputs callInputs{110.0, 100.0, 2.0, 0.025, 0.01, 0.2, "call", 250'000, 123};
    const MCInputs putInputs{110.0, 100.0, 2.0, 0.025, 0.01, 0.2, "put", 250'000, 123};

    BlackScholes callAnalytic(callInputs.s, callInputs.k, callInputs.t, callInputs.r,
                              callInputs.q, callInputs.sigma, "call");
    BlackScholes putAnalytic(putInputs.s, putInputs.k, putInputs.t, putInputs.r,
                             putInputs.q, putInputs.sigma, "put");

    const double referenceCall = callAnalytic.calculate();
    const double referencePut = putAnalytic.calculate();

    const double loopCall = euro_mc_for_loop(callInputs.s, callInputs.k, callInputs.t, callInputs.r,
                                             callInputs.q, callInputs.sigma, callInputs.type,
                                             callInputs.paths, callInputs.seed);
    const double multithreadCall = euro_mc_multithread(callInputs.s, callInputs.k, callInputs.t,
                                                       callInputs.r, callInputs.q, callInputs.sigma,
                                                       callInputs.type, callInputs.paths, callInputs.seed, 4);
    const double simdMTCall = euro_mc_simd_multithread(callInputs.s, callInputs.k, callInputs.t,
                                                       callInputs.r, callInputs.q, callInputs.sigma,
                                                       callInputs.type, callInputs.paths, callInputs.seed, 4);

    const double loopPut = euro_mc_for_loop(putInputs.s, putInputs.k, putInputs.t, putInputs.r,
                                            putInputs.q, putInputs.sigma, putInputs.type,
                                            putInputs.paths, putInputs.seed);

    CHECK(loopCall == Catch::Approx(referenceCall).margin(0.25));
    CHECK(multithreadCall == Catch::Approx(referenceCall).margin(0.35));
    CHECK(simdMTCall == Catch::Approx(referenceCall).margin(0.35));
    CHECK(loopPut == Catch::Approx(referencePut).margin(0.25));
}

TEST_CASE("Monte Carlo pricers reject invalid option type")
{
    CHECK_THROWS_AS(euro_mc_for_loop(100.0, 100.0, 1.0, 0.02, 0.0, 0.2, "straddle", 10, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(euro_mc_valarray(100.0, 100.0, 1.0, 0.02, 0.0, 0.2, "straddle", 10, 1),
                    std::invalid_argument);
}
