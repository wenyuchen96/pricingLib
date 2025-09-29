#include <iostream>
#include <vector>
#include <cmath>

// implement CRR BLM in a 2D vector
/*
binomial lattice model inputs
s: stock price
k: strike price
t: time to expiration
r: risk-free rate
sigma: volatility
n: number of steps
*/

// create a n by n matrix of 0s
std::vector<std::vector<double>> createZeroMatrix(double n)
{
    return std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
}

// binomial lattice model to value a european call option
double blm_euro(double s, double k, double t, double r, double sigma, int n, std::string type)
{
    // set up parameters
    double dt = t / static_cast<double>(n);
    double u = exp(sigma * sqrt(dt));       // size of upward movement
    double d = 1 / u;                       // size of downward movement
    double p = (exp(r * dt) - d) / (u - d); // risk-neutral probability of upward movement

    // create a binomial lattice
    auto stockPrice = createZeroMatrix(n + 1);
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            stockPrice[j][i] = s * pow(u, (i - j)) * pow(d, j);
        }
    }
    // calculation option values at maturity
    auto optionPrice = createZeroMatrix(n + 1);
    if (type == "call")
    {
        for (int j = 0; j <= n; j++)
        {
            optionPrice[j][n] = std::max(0.0, stockPrice[j][n] - k);
        }
    }
    else if (type == "put")
    {
        for (int j = 0; j <= n; j++)
        {
            optionPrice[j][n] = std::max(0.0, k - stockPrice[j][n]);
        }
    }
    // backward induction to calculate option value
    for (int i = (n - 1); i >= 0; i--)
    {
        for (int j = 0; j <= i; j++)
        {
            optionPrice[j][i] = exp(-r * dt) * (p * optionPrice[j][i + 1] + (1 - p) * optionPrice[j + 1][i + 1]);
        }
    }
    return optionPrice[0][0];
}
