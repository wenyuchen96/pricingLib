#include <catch2/catch_all.hpp>

#include "../models/hull_white_tree.cpp"

namespace
{
    struct Inputs
    {
        double sigma;
        double a;
        double tExp;
        double tMat;
        double strike;
        double face;
        std::vector<double> dfTimes;
        std::vector<double> dfValues;
    };

    ZCBOptionPrice analytic_price(const Inputs &in)
    {
        const double ptExp = HullWhiteTree::interpolate(in.tExp, in.dfTimes, in.dfValues);
        const double ptMat = HullWhiteTree::interpolate(in.tMat, in.dfTimes, in.dfValues);
        const double a = std::max(in.a, 1e-6);
        double sigmap = (in.sigma / a) * (1.0 - std::exp(-a * (in.tMat - in.tExp)));
        sigmap *= std::sqrt((1.0 - std::exp(-2.0 * a * in.tExp)) / (2.0 * a));
        sigmap = std::max(sigmap, 1e-6);
        const double h = std::log((in.face * ptMat) / (in.strike * ptExp)) / sigmap + 0.5 * sigmap;

        const double call = in.face * ptMat * HullWhiteTree::NormalCdf(h) -
                            in.strike * ptExp * HullWhiteTree::NormalCdf(h - sigmap);
        const double put = in.strike * ptExp * HullWhiteTree::NormalCdf(-h + sigmap) -
                           in.face * ptMat * HullWhiteTree::NormalCdf(-h);
        return ZCBOptionPrice{call, put};
    }
}

TEST_CASE("Hull-White analytic ZCB option matches manual computation")
{
    const Inputs params{0.01, 0.05, 1.5, 4.0, 0.98, 1.0,
                        {0.0, 1.0, 2.0, 3.0, 4.0},
                        {1.0, 0.99, 0.975, 0.955, 0.94}};

    HullWhiteTree tree(params.sigma, params.a, 20, params.dfTimes);
    const ZCBOptionPrice priced = tree.optionOnZCB(params.tExp, params.tMat, params.strike,
                                                   params.face, params.dfTimes, params.dfValues);
    const ZCBOptionPrice expected = analytic_price(params);

    CHECK(priced.call == Catch::Approx(expected.call).margin(1e-12));
    CHECK(priced.put == Catch::Approx(expected.put).margin(1e-12));
}

TEST_CASE("Hull-White option satisfies put-call parity for ZCBs")
{
    const Inputs params{0.015, 0.03, 0.75, 3.0, 0.97, 1.0,
                        {0.0, 0.5, 1.5, 3.0},
                        {1.0, 0.995, 0.985, 0.96}};

    HullWhiteTree tree(params.sigma, params.a, 10, params.dfTimes);
    const ZCBOptionPrice priced = tree.optionOnZCB(params.tExp, params.tMat, params.strike,
                                                   params.face, params.dfTimes, params.dfValues);

    const double ptExp = HullWhiteTree::interpolate(params.tExp, params.dfTimes, params.dfValues);
    const double ptMat = HullWhiteTree::interpolate(params.tMat, params.dfTimes, params.dfValues);

    const double parity = priced.call - priced.put - (params.face * ptMat - params.strike * ptExp);
    CHECK(parity == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("Hull-White input validation")
{
    HullWhiteTree tree(0.01, 0.05, 5, {0.0, 1.0});
    const std::vector<double> dfValues{1.0, 0.98};

    CHECK_THROWS_AS(tree.optionOnZCB(2.0, 1.0, 1.0, 1.0, {0.0, 1.0}, dfValues), std::invalid_argument);
    CHECK_THROWS_AS(tree.optionOnZCB(-0.1, 1.0, 1.0, 1.0, {0.0, 1.0}, dfValues), std::invalid_argument);
}
