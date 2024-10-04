% BOOTSTRAP_Y     Iterative procedure for structural estimation
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
function [YY,e,XX] = bootstrap_Y2(Yobs,cfg,theta,use_minimal_system,bs_type)

[negLL, x, P,ll, struct_shocks,v,K,xs,invF,minSYS]  = ...
    estimate_theta2(Yobs,cfg,theta,use_minimal_system,false,false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [A, b, C, D, B, Sig_eps] = ... get the system matrices
%         prepareFilter_loglin(theta,y,cfg); % continuous-time system matrices
%     cfg.growth_rates = false; % get repr for stationary data
%     [Ah, Sig_eta, C,D, x_0, T,nx,ny,y,cfg,sys] = ... discretize the model
%         DiscretizeSystem_wp4(y,cfg,A,B,C,D,Sig_eps,h); %
%     K = K(1:length(x_0),:,:);
%     if use_minimal_system
%         % rescale and store the ABCD system
%         sys.A = sys.A; sys.B=sys.B;
%         sys.C=scale(Sy,:)*sys.C; sys.D = scale(Sy,:)*sys.D;
%         sys.Sig_eps = scale(Sy,Sy).^2 .* sys.Sig_eps(Sy,Sy);
%         sys.Sig_eta = sys.Sig_eta;
% 
%         % find the minimal representation of the ABCD system
%         if ~check_minimality(sys.A,sys.B,sys.C) % check if system is already minimal
%             [CheckCO,minns,minSYS] = get_minimal_state_representation(sys, 0,false); % get the minimal representation of the system, get the prediction errors
%             minSYS.Sig_eta = minSYS.B  * minSYS.B';
%             minSYS.Sig_eps = minSYS.D  * minSYS.D';
%         else
%             minSYS = sys;
%         end
%     else
%                 C = scale(Sy,:)*C; Sig_eps = scale(Sy,Sy).^2 .* Sig_eps(Sy,Sy);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t0 = 1; %size(x,1)^2; 
    v = v(:,t0:end);
    [n,T]=size(v);
    missing = (v==0);  missing = logical([missing(:,1:end-1)]); v(missing)=nan;
    vbar = nanmean(v,2); v(isnan(v))=0;
    switch bs_type
        case 'nonparametric'
            CF = zeros(n,n,T); CiF=CF;
            for t=1:T-1
                missing = (v(:,t))==0;
                CF(~missing,~missing,t)  = chol(F(~missing,~missing,t),'lower');        % sqrt(F)
                CiF(~missing,~missing,t) = eye(sum(~missing))/CF(~missing,~missing,t);              % inv(sqrt(F))
                vz(:,t)    = CiF(:,:,t) * (v(:,t) - vbar);  %z score
            end
            vc = vz(:,randi(T-1,T-1,1));                    % sample errors
            for t=1:T-1
                e(:,t) = CF(:,:,t) * vc(:,t);               % bootstraped errors
            end
            if use_minimal_system
                e(:,T)=zeros(n,1);
            else
                e = [zeros(n,1), e];
            end
        case 'wild'
            vc = v - vbar; 
            vc(missing) = 0; % demean pred error - unnecessary
            w = randn(T,1)'; w = repmat(w,size(vc,1),1);    % random variable: Ew=0 Eww'=I
            e = vc .* w;                                    % Wild bootstrap
        otherwise
            error('Bootstrap type not recognised.')
    end
%     T = T-1; e=e(:,2:end); %K = K(:,:,2:end);
    e=e(:,1+ ~use_minimal_system:end-use_minimal_system); T=T-1;
    XX = zeros(length(x(:,1)),T);
    YY = zeros(size(e,1),T);
    XX(:,1) = x(:,t0); YY(:,1) = Yobs(:,t0); % starting values
    
    if use_minimal_system
        for t=1:T-1 % delayed: y[t|t-1] = C * x[t-1|t-1]; XX(:,t) is actually XX(:,t-1)
            XX(:,t+1) =  minSYS.A * XX(:,t) + K(:,:,t)*e(:,t);
            YY(:,t+1) =  minSYS.C * XX(:,t+1) + e(:,t+1);
        end
%         YY(:,end) = y(:,end);
%         YY = YY(:,1:end-1);
    else
        for t=1:T-1 % current: y[t|t] = C * x[t|t]
            XX(:,t+1) = minSYS.A * XX(:,t) + K(:,:,t)*e(:,t); 
            YY(:,t+1) = minSYS.C  * XX(:,t+1) + e(:,t+1);
        end

    end
    YY(missing) = nan;
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

