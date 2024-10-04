% ESTIMATE_THETA     Iterative procedure for structural estimation
% Inputs:
%   y       (double)    array, cointains data
%   cfg     (struct)    contains configuration settings.
%   theta   (double)    vector of structural parameters.
%
% Outputs:
%   negLL (double)    scalar, joint negative log-likelihood to minimize
%   x_hat (double)    array, filtered states
%   Sigma (double)    3-D array, filtered MSE
%   alpha (double)    array, smoothed states
%   ll    (double)    array, log-likelihod functions
%   ct    (struct)    struct of continuous-time infos
%% =======================================================================%
function [negLL, x, P,ll, struct_shocks,v,K,xs,invF,minSYS] = estimate_theta2(y,cfg,theta, ...
    use_minimal_system_repr, ...
    computation_of_se, ...
    is_bootstrapped_data)

if nargin<4
    use_minimal_system_repr=true;
end
if nargin<5
    computation_of_se=false;
end
if nargin<6
    is_bootstrapped_data=false;
end
if numel(theta)>numel(cfg.par_ident)
    error(['number of parameters to estimate is larger ...' ...
        'than the number of parameters specified in cfg.par_ident.'])
end
[cfg]   = check_sampling(cfg); % check info about data sampling
[theta] = ... transform parameters back into the parameters' space
    transform_parameters(theta,cfg,computation_of_se); % 

h = cfg.h;
Sy = cfg.Sy;
scale = cfg.scale;
eta   = theta(7);
etah  = eta*h;

err = 0;
try
    [A, b, C, d, B, Sig_eps, eigviol] = ... get the system matrices
        prepareFilter_loglin(theta,[],cfg); % continuous-time system matrices
    err = err + eigviol;
catch % for some parameter values there could be no solution
    err =1; % return
    negLL = 999999999999999;
    x=0; P=0; ll=0;
    struct_shocks=0;v=0;K=0;xs=0;invF=0;minSYS=0;
end

if ~err
    which_y_in_transition = cfg.info.take_differences;      % Choice of variables in log-differences
    idx = find(Sy.*which_y_in_transition*cfg.growth_rates); % which measurement gets included in the transition if growth_rates is true
    n_add_states = numel(idx);                              % of how much we should expand the transition equation
    cfg.n_add_states=n_add_states;
    d(idx) =  etah;                                 % add deterministic trend component
    d(4) = log(1+exp(d(4)));                                % we use gross interest rate


    %%%
    
    % cfg.minimal_representation=false;
% A = U*A/U; B = U*B; C=C/U;
    [Ah, Sig_eta, C,~,~,~,~,~,y,cfg,sys] = ... discretize the model
        DiscretizeSystem_wp4(y,cfg,A,B,C,d,Sig_eps,h); %
                         % rescale and resize C
    Sig_eps = cfg.Sig_eps;  
    sys.Bct = B; sys.Cct = cfg.Ct; sys.Act = A;
    minSYS=sys;
%     [idFlag] = isIdentified(theta,A,B,cfg.Ct,minSYS.A,minSYS.B,minSYS.C,cfg,y)

    %store
    if use_minimal_system_repr
        [CheckCO,minns,minSYS] = get_minimal_state_representation(minSYS, 0,false);
        Ah      = minSYS.A;
        C       = minSYS.C;
        minSYS.SIGMA = sys.SIGMA; 
        Sig_eps = minSYS.D*minSYS.SIGMA*minSYS.D';
        Sig_eta = minSYS.B*minSYS.SIGMA(1:size(minSYS.B,2),1:size(minSYS.B,2))*minSYS.B';
    end
    x_0 = zeros(size(Ah,1),1);
    nx  = size(Ah,1);
    P_0 = reshape((eye(nx^2) - kron(Ah,Ah))\Sig_eta(:),nx,nx);

    %% Kalman filter iteration
    if use_minimal_system_repr % delayed filter
        [x, P, loglik,ll, v, K, invF] = ...
            KalmanFilter_mixed_freq_ABCD(y, Ah, C, Sig_eta, Sig_eps, x_0, P_0,cfg);
    else % current filter
        [x, P, loglik,ll, v, K, invF] = ...
            KalmanFilter_mixed_freq(y, Ah, C, Sig_eta, Sig_eps, x_0, P_0,cfg);
        minSYS.A = Ah;
        minSYS.C = C;
        minSYS.Sig_eps = Sig_eps;
        minSYS.Sig_eta = Sig_eta;

    end
    if nargout>=8
%         [~, ~, ~,~, struct_shocks, ~, ~, ~, xs] =...
%             KalmanFilter_Smoother(y, Ah, C, Sig_eta, Sig_eps, x_0, P_0,cfg);
        [~, ~, ~,~, struct_shocks, ~, ~, ~, xs] = ...
            KalmanFilter_Smoother_MixedFreq_ABCD(y, Ah, C, Sig_eta, Sig_eps, x_0, P_0,cfg, use_minimal_system_repr);
    else
        xs=0;
        struct_shocks =0;
    end
    negLL = -loglik;

else
    negLL = 9999999999; x=0; P=0;ll=0; struct_shocks=0;v=0;K=0;xs=0;invF=0;minSYS=0;
end

end

function CheckCO = check_minimality(a,b,c)
%CHECKMINIMALITY Check minimality of the system
% This function computes the controllability and the observability matrices
%  of the ABCD system and checks if the system is minimal
%
% Inputs:
%   (a, b, c) Solution matrices A,B,C of ABCD representation of a DSGE model
% Outputs:
%   CheckCO     (bool)  equals 1, if both rank conditions for observability
%                       and controllability are fullfilled
%=========================================================================%
n = size(a,1);
CC = b; % Initialize controllability matrix
OO = c; % Initialize observability matrix
if n >= 2
    for indn = 1:1:n-1
        CC = [CC, (a^indn)*b]; % Set up controllability matrix
        OO = [OO; c*(a^indn)]; % Set up observability matrix
    end
end
CheckC = (rank(full(CC))==n);   % Check rank of controllability matrix
CheckO = (rank(full(OO))==n);   % Check rank of observability matrix
CheckCO = CheckC&CheckO;        % equals 1 if minimal state space
end % eof

function [theta] = transform_parameters(theta,cfg,is_for_ses)
tmp = zeros(size(cfg.pars));
tmp(logical(cfg.par_ident)) = theta;
tmp(~(cfg.par_ident)) = cfg.pars(~cfg.par_ident);
theta = tmp;
if ~is_for_ses % transform the parameters
    % transform parameters
    for i = 1:numel(theta)
        theta(i) = invtrans(theta(i),cfg.trspec(i));
    end
else % parameters are already transformed
    for i = 1:numel(theta)
        if cfg.par_ident(i)~=1
            theta(i) = invtrans(theta(i),cfg.trspec(i));
        else
            theta(i) = theta(i);
        end
    end
end
end

function [cfg] = check_sampling(cfg)
if ~isfield(cfg,'mixed_sampling')
    cfg.mixed_sampling=false;
end
if ~cfg.mixed_sampling % homogenous sampling of measurements
    if cfg.flows
        sampling = 'flow';
    else
        sampling = 'stock';
    end
    for i = 1:size(cfg.info,1)
        cfg.info.sampling{i}= sampling;
    end
end
end