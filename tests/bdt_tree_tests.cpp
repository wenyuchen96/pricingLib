#include <catch2/catch_all.hpp>

#include "../models/bdt_tree.cpp"

TEST_CASE("americanBondOption handles immediate exercise payoffs")
{
    const std::vector<double> cpnTimes{1.0};
    const std::vector<double> cpnFlows{0.0};
    const std::vector<double> dfTimes{0.0, 1.0};
    const std::vector<double> dfValues{1.0, 1.0};
    const std::vector<double> treeTimes{0.0, 1.0};
    const std::vector<std::vector<double>> Q{{1.0}, {0.5, 0.5}};
    const std::vector<std::vector<double>> rt{{0.0}, {0.0, 0.0}};
    const double dt = 1.0;

    const auto result = americanBondOption(1.0, 1.0, 95.0, 100.0,
                                           cpnTimes, cpnFlows, OptionType::AMERICAN,
                                           dfTimes, dfValues, treeTimes, Q, rt, dt);

    CHECK(result.call == Catch::Approx(95.0));
    CHECK(result.put == Catch::Approx(100.0));
}

TEST_CASE("BDTTree bond option prices are stable for calibrated curve")
{
    const double sigma = 0.18;
    const int steps = 4;
    const std::vector<double> marketTimes{0.5, 1.0, 1.5, 2.0, 2.5};
    const std::vector<double> marketDfs{0.99, 0.975, 0.955, 0.94, 0.92};

    BDTTree tree(sigma, steps, marketTimes, marketDfs);
    tree.buildTree();

    const std::vector<double> cpnTimes{0.5, 1.0, 1.5, 2.0, 2.5};
    const std::vector<double> cpnFlows{2.0, 2.0, 2.0, 2.0, 2.0};

    const BondOptionResult result = tree.bondOption(1.5, 100.0, 100.0,
                                                    cpnTimes, cpnFlows, OptionType::AMERICAN);

    CHECK(result.call == Catch::Approx(93.8399992904).margin(1e-6));
    CHECK(result.put == Catch::Approx(100.0).margin(1e-9));
}

TEST_CASE("BDTTree bondOption validates inputs")
{
    const double sigma = 0.15;
    const int steps = 3;
    const std::vector<double> marketTimes{0.5, 1.0, 1.5};
    const std::vector<double> marketDfs{0.99, 0.97, 0.94};

    BDTTree tree(sigma, steps, marketTimes, marketDfs);
    tree.buildTree();

    const std::vector<double> cpnTimes{0.5, 1.0, 1.5};
    const std::vector<double> cpnFlows{2.0, 2.0, 102.0};

    CHECK_THROWS_AS(tree.bondOption(2.0, 100.0, 100.0, cpnTimes, cpnFlows, OptionType::AMERICAN),
                    std::invalid_argument);
}
