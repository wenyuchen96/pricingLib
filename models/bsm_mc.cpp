#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <valarray>
#include <thread>
#include <chrono>
#include <algorithm>
#include <xsimd/xsimd.hpp>
#include <numeric>

/*
model inputs
s: stock price
k: strike price
t: time to expiration
r: risk-free rate
q: dividend yield
sigma: volatility
type: option type "call" or "put"
numOfPaths: number of paths
seed: seed number for ranodm number generation
*/

double euro_mc_for_loop(double s, double k, double t, double r, double q, double sigma, std::string type, int numOfPaths, int seed)
{
    std::mt19937 gen(seed);
    std::normal_distribution<> norm(0.0, 1.0);
    double payoffSum{};

    for (int i = 0; i < numOfPaths; i++)
    {
        double z{norm(gen)};
        double st = s * std::exp((r - q - 0.5 * sigma * sigma) * t + sigma * std::sqrt(t) * z);

        if (type == "call")
        {
            payoffSum += std::max(st - k, 0.0);
        }
        else if (type == "put")
        {
            payoffSum += std::max(k - st, 0.0);
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }
    }

    return (payoffSum / numOfPaths) * std::exp(-r * t);
}

double euro_mc_valarray(double s, double k, double t, double r, double q, double sigma,
                        const std::string &type, int numPaths, int seed)
{
    std::valarray<double> z(numPaths);
    std::mt19937 gen(seed);
    std::normal_distribution<> norm(0.0, 1.0);
    for (int i = 0; i < numPaths; ++i)
        z[i] = norm(gen);

    double drift = (r - q - 0.5 * sigma * sigma) * t;
    double vol_sqrt_t = sigma * std::sqrt(t);
    std::valarray<double> st = s * std::exp(drift + vol_sqrt_t * z);

    std::valarray<double> payoff(numPaths);
    if (type == "call")
    {
        std::transform(std::begin(st), std::end(st), std::begin(payoff), [=](double si)
                       { return std::max(si - k, 0.0); });
    }
    else if (type == "put")
    {
        std::transform(std::begin(st), std::end(st), std::begin(payoff), [=](double si)
                       { return std::max(k - si, 0.0); });
    }
    else
    {
        throw std::invalid_argument("Option type must be 'call' or 'put'");
    }

    double discount = std::exp(-r * t);
    return discount * (payoff.sum() / numPaths);
}

// Monte Carlo using xsimd
double euro_mc_xsimd(double s, double k, double t, double r, double q, double sigma,
                     const std::string &type, int numPaths, int seed)
{
    using b_type = xsimd::batch<double>;
    constexpr std::size_t simd_size = b_type::size;

    std::vector<double> z(numPaths);
    std::mt19937 gen(seed);
    std::normal_distribution<> norm(0.0, 1.0);
    for (auto &zi : z)
        zi = norm(gen);

    double drift = (r - q - 0.5 * sigma * sigma) * t;
    double vol_sqrt_t = sigma * std::sqrt(t);
    std::vector<double> payoff(numPaths);

    b_type zero(0.0);
    for (std::size_t i = 0; i + simd_size <= numPaths; i += simd_size)
    {
        b_type zb = xsimd::load_unaligned(&z[i]);
        b_type st = s * xsimd::exp(b_type(drift) + vol_sqrt_t * zb);
        b_type result;

        if (type == "call")
        {
            result = xsimd::max(st - b_type(k), zero);
        }
        else if (type == "put")
        {
            result = xsimd::max(b_type(k) - st, zero);
        }
        else
        {
            throw std::invalid_argument("Option type must be 'call' or 'put'");
        }

        result.store_unaligned(&payoff[i]);
    }

    // Handle any remaining paths not divisible by simd_size
    for (std::size_t i = (numPaths / simd_size) * simd_size; i < numPaths; ++i)
    {
        double st = s * std::exp(drift + vol_sqrt_t * z[i]);
        if (type == "call")
            payoff[i] = std::max(st - k, 0.0);
        else
            payoff[i] = std::max(k - st, 0.0);
    }

    double sum = std::accumulate(payoff.begin(), payoff.end(), 0.0);
    return std::exp(-r * t) * (sum / numPaths);
}

double worker_mc(double s, double k, double t, double r, double q, double sigma,
                 const std::string &type, int numPaths, int seed_offset)
{
    std::mt19937 gen(seed_offset);
    std::normal_distribution<> norm(0.0, 1.0);
    double drift = (r - q - 0.5 * sigma * sigma) * t;
    double vol_sqrt_t = sigma * std::sqrt(t);
    double local_sum = 0.0;

    for (int i = 0; i < numPaths; ++i)
    {
        double z = norm(gen);
        double st = s * std::exp(drift + vol_sqrt_t * z);
        if (type == "call")
            local_sum += std::max(st - k, 0.0);
        else if (type == "put")
            local_sum += std::max(k - st, 0.0);
    }

    return local_sum;
}

double euro_mc_multithread(double s, double k, double t, double r, double q, double sigma,
                           const std::string &type, int numPaths, int seed, int numThreads)
{
    std::vector<std::thread> threads;
    std::vector<double> partial_sums(numThreads, 0.0);
    int pathsPerThread = numPaths / numThreads;

    for (int i = 0; i < numThreads; ++i)
    {
        threads.emplace_back([&, i]()
                             { partial_sums[i] = worker_mc(s, k, t, r, q, sigma, type, pathsPerThread, seed + i); });
    }

    for (auto &th : threads)
        th.join();

    double total = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
    double discount = std::exp(-r * t);
    return discount * (total / numPaths);
}

// xsimd + multi-thread

double simd_worker(double s, double k, double t, double r, double q, double sigma,
                   const std::string &type, int numPaths, int seed)
{
    using b_type = xsimd::batch<double>;
    constexpr std::size_t simd_size = b_type::size;

    std::vector<double> z(numPaths);
    std::mt19937 gen(seed);
    std::normal_distribution<> norm(0.0, 1.0);
    for (auto &zi : z)
        zi = norm(gen);

    double drift = (r - q - 0.5 * sigma * sigma) * t;
    double vol_sqrt_t = sigma * std::sqrt(t);
    std::vector<double> payoff(numPaths);

    b_type zero(0.0);
    for (std::size_t i = 0; i + simd_size <= numPaths; i += simd_size)
    {
        b_type zb = xsimd::load_unaligned(&z[i]);
        b_type st = s * xsimd::exp(b_type(drift) + vol_sqrt_t * zb);
        b_type result;

        if (type == "call")
            result = xsimd::max(st - b_type(k), zero);
        else
            result = xsimd::max(b_type(k) - st, zero);

        result.store_unaligned(&payoff[i]);
    }

    // Handle leftovers
    for (std::size_t i = (numPaths / simd_size) * simd_size; i < numPaths; ++i)
    {
        double st = s * std::exp(drift + vol_sqrt_t * z[i]);
        if (type == "call")
            payoff[i] = std::max(st - k, 0.0);
        else
            payoff[i] = std::max(k - st, 0.0);
    }

    return std::accumulate(payoff.begin(), payoff.end(), 0.0);
}

double euro_mc_simd_multithread(double s, double k, double t, double r, double q, double sigma,
                                const std::string &type, int numPaths, int seed, int numThreads)
{
    std::vector<std::thread> threads;
    std::vector<double> partial_sums(numThreads, 0.0);
    int pathsPerThread = numPaths / numThreads;

    for (int i = 0; i < numThreads; ++i)
    {
        threads.emplace_back([&, i]()
                             { partial_sums[i] = simd_worker(s, k, t, r, q, sigma, type, pathsPerThread, seed + i); });
    }

    for (auto &th : threads)
        th.join();

    double total = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
    return std::exp(-r * t) * (total / numPaths);
}

int main()
{
    double s = 10.5;
    double k = 5.5;
    double t = 2.0;
    double r = 0.045;
    double q = 0.0;
    double sigma = 0.6;
    std::string type = "call";
    int numOfPaths = 100'000;
    int seed = 50;

    // For loops
    //  Start timing
    auto start1 = std::chrono::high_resolution_clock::now();
    double price1 = euro_mc_for_loop(s, k, t, r, q, sigma, type, numOfPaths, seed);
    // End timing
    auto end1 = std::chrono::high_resolution_clock::now();
    // Compute elapsed time in milliseconds
    std::chrono::duration<double, std::milli> duration1 = end1 - start1;
    // Output
    std::cout << "For Loop Monte Carlo Price: " << price1 << '\n';
    std::cout << "Elapsed time: " << duration1.count() << " ms\n";

    // Valarray
    //  Start timing
    auto start2 = std::chrono::high_resolution_clock::now();
    double price2 = euro_mc_valarray(s, k, t, r, q, sigma, type, numOfPaths, seed);
    // End timing
    auto end2 = std::chrono::high_resolution_clock::now();
    // Compute elapsed time in milliseconds
    std::chrono::duration<double, std::milli> duration2 = end2 - start2;
    // Output
    std::cout << "Vectorized Valarray Monte Carlo Call Price: " << price2 << '\n';
    std::cout << "Elapsed time: " << duration2.count() << " ms\n";

    // SIMD
    //  Start timing
    auto start3 = std::chrono::high_resolution_clock::now();
    double price3 = euro_mc_xsimd(s, k, t, r, q, sigma, type, numOfPaths, seed);
    // End timing
    auto end3 = std::chrono::high_resolution_clock::now();
    // Compute elapsed time in milliseconds
    std::chrono::duration<double, std::milli> duration3 = end3 - start3;
    // Output
    std::cout << "SIMD Monte Carlo Call Price: " << price3 << '\n';
    std::cout << "Elapsed time: " << duration3.count() << " ms\n";

    // multi-thread mc
    int numThreads = std::thread::hardware_concurrency(); // auto-detect cores
    auto start4 = std::chrono::high_resolution_clock::now();
    double price4 = euro_mc_multithread(s, k, t, r, q, sigma, type, numOfPaths, seed, numThreads);
    auto end4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration4 = end4 - start4;
    std::cout << "Multithreaded MC Price: " << price4 << "\n";
    std::cout << "Elapsed time: " << duration4.count() << " ms\n";

    // xsimd + multi-thread
    // int numThreads = std::thread::hardware_concurrency();
    auto start5 = std::chrono::high_resolution_clock::now();
    double price5 = euro_mc_simd_multithread(s, k, t, r, q, sigma, type, numOfPaths, seed, numThreads);
    auto end5 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration5 = end5 - start5;
    std::cout << "SIMD + Multithreaded MC Price: " << price5 << '\n';
    std::cout << "Elapsed time: " << duration5.count() << " ms\n";
}