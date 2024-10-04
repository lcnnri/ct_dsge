switch Table
    case '4'
        if strcmpi(Panel,'b')
            include_measurement_error = true;
        else
            include_measurement_error = false;
        end
        %%
        models = ... models to use in estimation
            {... % EDM for mixed sampling of measurements
            'F-SSR';   % EDM for flow measurements
            'S-SSR';   % EDM for stock measurements
            'EM-SSR'}; % Euler-Maruyama discretization
        % models={'S-SSR'};
        bootstrap_inference = ... construct std errors; if false, se are way too large
            true;

        which_meas = ...  data analogues of which model measurements to use
            {'c' 1 ;  % consumption
            'n' 1 ;  % hours worked
            'y' 0 ;  % output
            'r' 0};  % interest rate

        what_freq = ... what frequency of data for estimation \in {'m','q'}; if both 'm' and 'q' -> mixed frequency
            ['q'; % c
            'q'; % n
            'q'; % y
            'q'];% r



        est_std_meas_error = ... 1 estimate std err. of meas.; 0 calibrates it
            [false;
            false;
            false;
            true];
        %% parameters settings
        structpars_to_estimate = ... structural parameters to estimate
            [  0 % psi (gamma)
            0 % alpha
            0 % delta
            0 % rho
            1 % sigz
            1 % rhoz
            0 % eta
            1 % sigk
            ];

        trspec = ... 0:(-Inf,+Inf); 1:[0,1]; 2:[0,+Inf)
            [2,1,1,1,1,1,0,1];       % bounds of parameters

        meas_error_std = ... starting/calibrated value std of measurement error, scalar or column vector
            [0 0 0 0]';
        meas_error_std= [0.0061;
            0.0073; %73
            0;
            0.0034];%

        minimal_system_repr = ... use minimal system representation for estimation and bootstrap
            true;

        %% estimation settings
        knitro_toolbox=false;

        Tol = ...  optimization tolerance: XTol and fTol
            1e-9;

        is_data_stationary = ... true if data provided is already stationary
            false;

        growth_rates   = ... use differenced data
            false;

        diffusive_start = ... use the diffusive start algorithm (unused)
            false;

        plot_data = ... plot the data for sanity?
            true;

        %% bootstrap settings
        BS = ... number of bootstrap samples (1000 in the paper)
            100;

        num_workers = ... max number of workers for bootstrap
            8;

        compute_hessians = ... % compute se using asymptotic BS se
            false;

        alpha = ... 1-alpha: nominal size of ci
            0.05;

        bs_type = ... bootstrap type: {'wild'; 'nonparametric'}
            'wild';
    case '5'
        if strcmpi(Panel,'b')
            include_measurement_error = true;
        else
            include_measurement_error = false;
        end

        models = ... models to use in estimation
            {... % EDM for mixed sampling of measurements
            'F-SSR';   % EDM for flow measurements
            'S-SSR';   % EDM for stock measurements
            'EM-SSR'}; % Euler-Maruyama discretization
        % models={'S-SSR'};
        bootstrap_inference = ... construct std errors; if false, se are way too large
            true;

        which_meas = ...  data analogues of which model measurements to use
            {'c' 1 ;  % consumption
            'n' 1 ;  % hours worked
            'y' 0 ;  % output
            'r' 0};  % interest rate

        what_freq = ... what frequency of data for estimation \in {'m','q'}; if both 'm' and 'q' -> mixed frequency
            ['m'; % c
            'q'; % n
            'q'; % y
            'q'];% r

        est_std_meas_error = ... 1 estimate std err. of meas.; 0 calibrates it
            [false;
            false;
            false;
            false];
        %% parameters settings
        structpars_to_estimate = ... structural parameters to estimate
            [  0 % psi (gamma)
            0 % alpha
            0 % delta
            0 % rho
            1 % sigz
            1 % rhoz
            0 % eta
            1 % sigk
            ];

        trspec = ... 0:(-Inf,+Inf); 1:[0,1]; 2:[0,+Inf)
            [2,1,1,1,1,1,0,1];       % bounds of parameters

        meas_error_std= [0.0061;
            0.0073; %73
            0;
            0.0034];%

        minimal_system_repr = ... use minimal system representation for estimation and bootstrap
            true;

        %% estimation settings
        knitro_toolbox = ... use knitro or fminunc
            false;

        Tol = ...  optimization tolerance: XTol and fTol
            1e-9;

        is_data_stationary = ... true if data provided is already stationary
            false;

        growth_rates   = ... use differenced data
            false;

        diffusive_start = ... use the diffusive start algorithm (unused)
            false;

        plot_data = ... plot the data for sanity?
            true;

        %% bootstrap settings
        BS = ... number of bootstrap samples (1000 in the paper)
            1000;

        num_workers = ... max number of workers for bootstrap
            8;

        compute_hessians = ... % compute se using asymptotic BS se
            false;

        alpha = ... 1-alpha: nominal size of ci
            0.05;

        bs_type = ... bootstrap type: {'wild'; 'nonparametric'}
            'wild';
    case '7'
        include_measurement_error = true;
        if strcmpi(Panel,'b')
            est_std_meas_error = ... 1 estimate std err. of meas.; 0 calibrates it
                [false;
                false;
                false;
                true];
        else
            est_std_meas_error = ... 1 estimate std err. of meas.; 0 calibrates it
                [true;
                true;
                true;
                true];
        end

        models = ... models to use in estimation
            {'MX-SSR'; % EDM for mixed sampling of measurements
            'F-SSR';   % EDM for flow measurements
            'S-SSR';   % EDM for stock measurements
            'EM-SSR'}; % Euler-Maruyama discretization

        bootstrap_inference = ... construct std errors; if false, se are way too large
            true;

        which_meas = ...  data analogues of which model measurements to use
            {'c' 1 ;  % consumption
            'n' 1 ;  % hours worked
            'y' 0 ;  % output
            'r' 1};  % interest rate

        what_freq = ... what frequency of data for estimation \in {'m','q'}; if both 'm' and 'q' -> mixed frequency
            ['q'; % c
            'q'; % n
            'q'; % y
            'q'];% r

        % include_measurement_error = ...
        %     true;

        % est_std_meas_error = ... 1 estimate std err. of meas.; 0 calibrates it
        %     [true;
        %      true;
        %      false;
        %      true];
        %% parameters settings
        structpars_to_estimate = ... structural parameters to estimate
            [  0 % psi (gamma)
            0 % alpha
            0 % delta
            0 % rho
            1 % sigz
            1 % rhoz
            0 % eta
            1 % sigk
            ];

        trspec = ... 0:(-Inf,+Inf); 1:[0,1]; 2:[0,+Inf)
            [2,1,1,1,1,1,0,1];       % bounds of parameters

        meas_error_std = ... starting/calibrating value std of measurement error, scalar or column vector
            [0 0 0 0]';
        meas_error_std =  [0.0061;
            0.0073; %73
            0;
            0.000];%

        minimal_system_repr = ... use minimal system representation for estimation and bootstrap
            true;

        %% estimation settings
        knitro_toolbox = ... use knitro or fminunc
            true;

        Tol = ...  optimization tolerance: XTol and fTol
            1e-9;

        is_data_stationary = ... true if data provided is already stationary
            false;

        growth_rates   = ... use differenced data
            false;

        diffusive_start = ... use the diffusive start algorithm (unused)
            false;

        plot_data = ... plot the data for sanity?
            true;

        %% bootstrap settings
        BS = ... number of bootstrap samples (1000 in the paper)
            100;

        num_workers = ... max number of workers for bootstrap
            8;

        compute_hessians = ... % compute se using asymptotic BS se
            false;

        alpha = ... 1-alpha: nominal size of ci
            0.05;

        bs_type = ... bootstrap type: {'wild'; 'nonparametric'}
            'wild';
end