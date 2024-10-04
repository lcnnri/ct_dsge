function S = paramstructbase

% -------------------------------------------------------------------------
% Set all the parameters values here
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Annual calibration
% -------------------------------------------------------------------------
S = struct('n_v', 1,                 ... % Number of controls
           'n_g', 2,                 ... % Number of states (both endogenous and exogenous)
           'n_p', 1,                 ... % Number of static equations
           'nEErrors', 1,            ... % Number of expectational errors
           'nShocks',2,              ... % Number of structural shocks
           'gamma', 2.6860,          ... % Weight of labor in utility (consistent with a 1/3 of working in the steady state)
           'alpha', 0.3,             ... % Capital share in output (consistent with the annual labor share of income of 0.7 in steady state)
           'delta', 0.06,            ... % Depreciation rate of capital (consistent with the steady state investment share in output, 0.2, and the annual steady state capital-output ratio, 2.5)
           'rho',   0.03,            ... % Subjective discount rate (consistent with a net return of 4% per annum in steady state)
           'sigz',  0.007*sqrt(4),   ... % Volatility of TFP (consistent with the quarterly volatility in post-war period for the US)
           'rhoz',  -log(0.95)*4,    ... % Persistence of TFP (consistent with the quarterly persistence in post-war period for the US)
           'eta',   0.02,            ... % Labor-augmenting growth rate (consistent with the average annual growth in post-war period for the US)
           'sigk',  0.0052*sqrt(4),  ... % Volatility depreciation shocks (from Ambler and Paquet, 1994)
           'freq',  1/4,             ... % Sampling frequency, e.g., 1/12 is monthly, 1/4 is quarterly and 1 is annualy
           'Tsim',  250000,          ... % Number of years to simulate
           'prec',  1/120            ... % Sampling precision (should be a multiple of the observation frequency
            );
return