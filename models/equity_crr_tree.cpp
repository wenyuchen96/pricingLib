#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// value european and american option based on the Cox-Ross-Rubinstein ("CRR") binomial model

class CRRTree
{
public:
    CRRTree(double stockPrice, double interestRate, double dividendRate, double volatility, int numStepsPerYear, double timeToExpiry, std::string optionType, double strikePrice, bool isEven) : stockPrice_{stockPrice}, interestRate_{interestRate}, dividendRate_{dividendRate}, volatility_{volatility}, numStepsPerYear_{numStepsPerYear}, timeToExpiry_{timeToExpiry}, optionType_{optionType}, strikePrice_{strikePrice}, isEven_{isEven} {};

    double calculate()
    {
        // number of steps over the term of the underlying instrument
        int numSteps{static_cast<int>(numStepsPerYear_ * timeToExpiry_)};
        numSteps = std::max(numSteps, 30);

        if (numSteps % 2 == 0 && isEven_ == false)
        {
            numSteps += 1; // was even, but we want odd so add 1
        }
        else if (numSteps % 2 == 1 && isEven_ == true)
        {
            numSteps += 1; // was odd, but we want even so add 1
        }

        // index formula to locate nodes in a 1D vector
        //  index (i,j) = i(i+1)/2 + j
        double dt = timeToExpiry_ / numSteps;              // dt is the size of a step
        int numNodes{(numSteps + 1) * (numSteps + 2) / 2}; // total number of nodes
        std::vector<double> stockValues(numNodes, 0.0);    // stock value matches the size of nodes
        stockValues[0] = stockPrice_;                      // initalize the stock price

        std::vector<double> optionValues(numNodes, 0.0);
        double u = exp(volatility_ * sqrt(dt));
        double d = 1.0 / u;
        double sLow{stockPrice_}; // lowest stock price, or bottom node at each timestep
        std::vector<double> uProbs(numSteps, 0.0);
        std::vector<double> periodDiscountRate(numSteps, 0.0);
        // populate probability tree and discount rate tree
        for (int iTime = 0; iTime < numSteps; iTime++)
        {
            double a = exp((interestRate_ - dividendRate_) * dt);
            uProbs[iTime] = (a - d) / (u - d);
            periodDiscountRate[iTime] = exp(-interestRate_ * dt);
        }
        // populate stock price tree
        for (int iTime = 1; iTime < numSteps + 1; iTime++)
        {
            sLow *= d;
            double stockPriceHolder{sLow};
            for (int iNode = 0; iNode < iTime + 1; iNode++)
            {
                int index{iTime * (iTime + 1) / 2};
                stockValues[static_cast<int>(index + iNode)] = stockPriceHolder;
                stockPriceHolder = stockPriceHolder * (u * u);
            }
        }

        // starting backward induction process
        // at the expiration date
        int index{numSteps * (numSteps + 1) / 2};
        for (int iNode = 0; iNode < numSteps + 1; iNode++)
        {
            double stockPriceHolder{stockValues[index + iNode]};

            if (optionType_ == "EUROPEAN_CALL")
            {
                optionValues[index + iNode] = std::max(stockPriceHolder - strikePrice_, 0.0);
            }
            else if (optionType_ == "EUROPEAN_PUT")
            {
                optionValues[index + iNode] = std::max(stockPrice_ - stockPriceHolder, 0.0);
            }
            else if (optionType_ == "AMERICAN_CALL")
            {
                optionValues[index + iNode] = std::max(stockPriceHolder - strikePrice_, 0.0);
            }
            else if (optionType_ == "AMERICAN_PUT")
            {
                optionValues[index + iNode] = std::max(stockPrice_ - stockPriceHolder, 0.0);
            }
        }
        // backward induction before the expiration date
        for (int iTime = numSteps - 1; iTime >= 0; --iTime)
        {
            int bIndex{iTime * (iTime + 1) / 2};

            for (int iNode = 0; iNode < iTime + 1; iNode++)
            {
                double stockPriceHolder{stockValues[bIndex + iNode]};
                double exerciseValue{0.0};

                // value of pre-exercise the option at timestep bIndex
                if (optionType_ == "EUROPEAN_CALL")
                {
                    exerciseValue = 0.0;
                }
                else if (optionType_ == "EUROPEAN_PUT")
                {
                    exerciseValue = 0.0;
                }
                else if (optionType_ == "AMERICAN_CALL")
                {
                    exerciseValue = std::max(stockPriceHolder - strikePrice_, 0.0);
                }
                else if (optionType_ == "AMERICAN_PUT")
                {
                    exerciseValue = std::max(strikePrice_ - stockPriceHolder, 0.0);
                }

                int bNextIndex{(iTime + 1) * (iTime + 2) / 2};
                double nextDownValue{optionValues[bNextIndex + iNode]};
                double nextUpValue{optionValues[bNextIndex + iNode + 1]};
                double futureEV{uProbs[iTime] * nextUpValue + (1.0 - uProbs[iTime]) * nextDownValue};
                double holdValue{periodDiscountRate[iTime] * futureEV};

                // select max(pre-exercise value, hold value)
                if (optionType_ == "EUROPEAN_CALL")
                {
                    optionValues[bIndex + iNode] = holdValue;
                }
                else if (optionType_ == "EUROPEAN_PUT")
                {
                    optionValues[bIndex + iNode] = holdValue;
                }
                else if (optionType_ == "AMERICAN_CALL")
                {
                    optionValues[bIndex + iNode] = std::max(exerciseValue, holdValue);
                }
                else if (optionType_ == "AMERICAN_PUT")
                {
                    optionValues[bIndex + iNode] = std::max(exerciseValue, holdValue);
                }
            }
        }
        return optionValues[0];
    }

    double calculateAvgEvenOddTree()
    {
        // Run with even number of steps
        CRRTree evenTree(stockPrice_, interestRate_, dividendRate_, volatility_,
                         numStepsPerYear_, timeToExpiry_, optionType_,
                         strikePrice_, true);
        double value1 = evenTree.calculate();

        // Run with odd number of steps
        CRRTree oddTree(stockPrice_, interestRate_, dividendRate_, volatility_,
                        numStepsPerYear_, timeToExpiry_, optionType_,
                        strikePrice_, false);
        double value2 = oddTree.calculate();

        return (value1 + value2) / 2.0;
    }

private:
    double stockPrice_;
    double interestRate_;
    double dividendRate_;
    double volatility_;
    int numStepsPerYear_;
    double timeToExpiry_;
    std::string optionType_;
    double strikePrice_;
    bool isEven_;
};

int main()
{
    // Option and model inputs
    double s = 10.5;                    // Stock price
    double k = 5.5;                     // Strike price
    double t = 2.0;                     // Time to expiry (in years)
    double r = 0.045;                   // Risk-free rate (continuously compounded)
    double q = 0.0;                     // Dividend yield
    double sigma = 0.6;                 // Volatility
    std::string type = "AMERICAN_CALL"; // Option type for CRR model
    int steps_per_year = 1000;          // Binomial steps per year

    // Create a CRRTree instance (even vs. odd controlled inside avg function)
    CRRTree tree(s, r, q, sigma, steps_per_year, t, type, k, true); // isEven doesn't matter here

    // Calculate option value using average of even/odd trees
    double price = tree.calculateAvgEvenOddTree();

    std::cout << type << " option price (CRR tree avg) is $" << price << "\n";

    return 0;
}