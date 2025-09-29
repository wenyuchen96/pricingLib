#include <catch2/catch_all.hpp>

#include <numeric>

#include "../models/blm.cpp"

namespace
{
    long long binomial_coeff(int n, int k)
    {
        if (k < 0 || k > n)
            return 0;
        if (k == 0 || k == n)
            return 1;
        long long res = 1;
        for (int i = 1; i <= k; ++i)
        {
            res = res * (n - (k - i));
            res /= i;
        }
        return res;
    }

    double closed_form_binomial_call(double s, double k, double t, double r, double sigma, int n)
    {
        const double dt = t / static_cast<double>(n);
        const double u = std::exp(sigma * std::sqrt(dt));
        const double d = 1.0 / u;
        const double p = (std::exp(r * dt) - d) / (u - d);
        double price = 0.0;
        for (int j = 0; j <= n; ++j)
        {
            const double sT = s * std::pow(u, n - j) * std::pow(d, j);
            const double payoff = std::max(0.0, sT - k);
            const double prob = binomial_coeff(n, j) * std::pow(p, n - j) * std::pow(1 - p, j);
            price += prob * payoff;
        }
        return std::exp(-r * t) * price;
    }
}

TEST_CASE("createZeroMatrix returns correct shape and zeros")
{
    const auto mat = createZeroMatrix(3);
    REQUIRE(mat.size() == 3);
    for (const auto &row : mat)
    {
        REQUIRE(row.size() == 3);
        for (double value : row)
        {
            CHECK(value == 0.0);
        }
    }
}

TEST_CASE("Binomial lattice matches closed-form binomial pricing")
{
    const double s = 100.0;
    const double k = 95.0;
    const double t = 0.75;
    const double r = 0.03;
    const double sigma = 0.2;
    const int n = 5;

    const double expected = closed_form_binomial_call(s, k, t, r, sigma, n);
    const double lattice = blm_euro(s, k, t, r, sigma, n, "call");

    CHECK(lattice == Catch::Approx(expected).margin(1e-9));
}

TEST_CASE("Binomial lattice handles put-call parity for puts")
{
    const double s = 100.0;
    const double k = 105.0;
    const double t = 1.0;
    const double r = 0.04;
    const double sigma = 0.25;
    const int n = 6;

    const double call = blm_euro(s, k, t, r, sigma, n, "call");
    const double put = blm_euro(s, k, t, r, sigma, n, "put");
    const double parity = call - put - s + k * std::exp(-r * t);

    CHECK(parity == Catch::Approx(0.0).margin(1e-9));
}
