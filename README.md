# pricingLib

Collection of standalone C++ pricing model prototypes for equity and rates derivatives. Each program focuses on a specific model.

## Available models

- `models/bachelier.cpp`: Analytic Bachelier pricing for calls and puts, including conversion from spot inputs.
- `models/black.cpp`: Black 1976 forward or futures option pricer with analytic Greeks.
- `models/black_shifted.cpp`: Shifted Black variant for negative rates with a built-in demo.
- `models/bsm.cpp`: Black-Scholes-Merton analytic pricer returning value and Greeks.
- `models/bsm_mc.cpp`: Monte Carlo engine showing scalar, valarray, SIMD (via xsimd), and multithreaded pricing pipelines.
- `models/blm.cpp`: Cox-Ross-Rubinstein binomial lattice for European and American vanilla options.
- `models/equity_crr_tree.cpp`: Flexible CRR tree supporting even or odd step counts and American exercise.
- `models/bdt_tree.cpp`: Black-Derman-Toy short-rate tree calibration plus callable and putable bond valuation utilities.
- `models/hull_white_tree.cpp`: One-factor Hull-White lattice with zero coupon bond analytics and tree-based pricing.
- `models/bk_tree.cpp`: Placeholder for a Black-Karasinski implementation (currently empty).
