#include <catch2/catch_all.hpp>

#include "../models/equity_crr_tree.cpp"

TEST_CASE("American call without dividends matches European call")
{
    const double s = 100.0;
    const double r = 0.02;
    const double q = 0.0;
    const double sigma = 0.25;
    const int stepsPerYear = 2;
    const double t = 1.0;
    const double k = 100.0;

    CRRTree euroCall(s, r, q, sigma, stepsPerYear, t, "EUROPEAN_CALL", k, true);
    CRRTree amerCall(s, r, q, sigma, stepsPerYear, t, "AMERICAN_CALL", k, true);

    const double euro = euroCall.calculate();
    const double american = amerCall.calculate();

    CHECK(american == Catch::Approx(euro).margin(1e-8));
}

TEST_CASE("American put is at least as valuable as European put")
{
    const double s = 90.0;
    const double r = 0.01;
    const double q = 0.0;
    const double sigma = 0.3;
    const int stepsPerYear = 2;
    const double t = 2.0;
    const double k = 100.0;

    CRRTree euroPut(s, r, q, sigma, stepsPerYear, t, "EUROPEAN_PUT", k, false);
    CRRTree amerPut(s, r, q, sigma, stepsPerYear, t, "AMERICAN_PUT", k, false);

    const double euro = euroPut.calculate();
    const double american = amerPut.calculate();

    CHECK(american >= euro);
}

TEST_CASE("Average even/odd tree matches explicit average")
{
    const double s = 120.0;
    const double r = 0.015;
    const double q = 0.005;
    const double sigma = 0.2;
    const int stepsPerYear = 3;
    const double t = 1.5;
    const double k = 110.0;

    CRRTree evenTree(s, r, q, sigma, stepsPerYear, t, "EUROPEAN_CALL", k, true);
    CRRTree oddTree(s, r, q, sigma, stepsPerYear, t, "EUROPEAN_CALL", k, false);
    CRRTree averaged(s, r, q, sigma, stepsPerYear, t, "EUROPEAN_CALL", k, true);

    const double evenVal = evenTree.calculate();
    const double oddVal = oddTree.calculate();
    const double expectedAverage = 0.5 * (evenVal + oddVal);

    CHECK(averaged.calculateAvgEvenOddTree() == Catch::Approx(expectedAverage).margin(1e-8));
}
