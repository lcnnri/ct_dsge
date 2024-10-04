% configurations for estimation ------------------------------------------%
if ~exist('knitro_toolbox')
    knitro_toolbox=false;
end
if ~exist('Tol')
    Tol=1e-9;
end
if ~exist("is_data_stationary")
    is_data_stationary=true;
end
if ~exist("growth_rates")
    growth_rates=true;
end
if ~exist("diffusive_start")
    diffusive_start=false;
end
if ~exist('bootstrap_inference')
    bootstrap_inference=false;
end
if ~exist('compute_hessians')
    compute_hessians=false;
end
if ~exist('which_meas');
    which_meas = {'c' 1; 'n' 1; 'y' 0; 'r' 0};
end
if ~exist('what_freq');
    what_freq = ['q'; 'q'; 'q'; 'q'];
end
if ~exist('meas_error_std')
    meas_error_std = [0 0 0 0]';
end
if ~exist('include_measurement_error','var')
    include_measurement_error = false;
end
if ~exist('est_std_meas_error','var')
    est_std_meas_error = [false false false false]';
end

if knitro_toolbox
    options = knitro_options(   'hessopt',	2, ... % pseudo-newton BFGS
        'algorithm', 1,... % 0: Auto; 1: Interior/Direct; 2: Interior/CG; 3: Active set; 4: SQP; 5: Multiple algos;
        'xtol', Tol, 'ftol', Tol,...
        'linesearch', 2,... % cubic interpolation scheme
        'outlev',0,... % output
        'linsolver', 4); %  HSL MA27 sparse symmetric indefinite solver
    options_no_verbose = options;
    options_no_verbose.outlev = 0;
else % fminunc
    options = optimoptions('fmincon',...
        'Display','iter',  ...
        'OptimalityTolerance',Tol, 'StepTolerance', Tol);
    options_no_verbose = options;
    options_no_verbose.Display = 'off';
end
warning('off')

% Construct InfoTable ====================================================%
%=========================================================================%
S = paramstructbase;

for ny = 1:size(which_meas,1)
    cfg.measurement(ny,1) = logical(which_meas{ny,2});
    ixm(ny,1) = strcmpi(what_freq(ny),'m');
    ixq(ny,1) = strcmpi(what_freq(ny),'q');
end

mixed_frequency=false;
if any(what_freq(cfg.measurement)=='m') % get the highest frequency available
    if any(what_freq(cfg.measurement)=='q')
        mixed_frequency=true;
    end
    S.freq = 1/12;
elseif any(what_freq(cfg.measurement)=='d')
    S.freq = 1/260;
end

if ~cfg.measurement(4) % if r is not used don't run the mixed sampling model
    if exist('models','var')
        idx = strcmpi(models(:),'MX-SSR');
        models = models(~idx);
        clear idx
    end
end
% Measurement Info Table
string_measurement = which_meas(:,1);
details = {'consumption'; 'hours worked'; 'output'; 'interest rate'};
take_differences= [true; false; true; false];
sampling = {'flow'; 'flow'; 'flow';'stock'};
used_in_estimation = cfg.measurement;
with_meas_error    = logical(meas_error_std~=0);
cfg.info = table(string_measurement, ...
    details, ...
    take_differences, ...
    sampling, ...
    used_in_estimation, ...
    with_meas_error);

% info about measurement errors

meas_error_std = include_measurement_error * meas_error_std(cfg.measurement);
est_std_meas_error = include_measurement_error * est_std_meas_error;
theta0 = [theta0; meas_error_std];
trspec(numel(structpars_to_estimate)+1:numel(theta0)) = 2;
structpars_to_estimate = [structpars_to_estimate;
    est_std_meas_error(cfg.measurement).*(cfg.measurement(cfg.measurement))];

% store cfg
cfg.measErr         = true; %always include, for construction of cov matrix
cfg.par_ident       = structpars_to_estimate;
cfg.pars            = theta0;
cfg.h               = S.freq;    % frequency of observations
cfg.freq            = what_freq; % monthly or quarterly
cfg.is_stationary   = is_data_stationary;
cfg.growth_rates    = growth_rates;
cfg.t0              = 0; % starting point for likelihood
cfg.is_diffusive    = diffusive_start;
cfg.trspec          = trspec;
cfg.bootstrap_inference = bootstrap_inference;
cfg.compute_hessians = compute_hessians;
h = cfg.h; %cfg.Delta=S.freq;
pd = makedist('Normal');
cfg.z_ = abs(icdf(pd,alpha/2)); % z score (critical value of standard normal)
cfg.trspec = trspec;
