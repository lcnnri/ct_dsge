function [x, P, loglik,ll, v, K, invF,FF] = KalmanFilter_mixed_freq_ABCD(y, A, C, Sig_eta, Sig_eps, x_init, P_init,cfg)
% KALMANFILTER_MIXED_FREQ Runs the Kalman filter for mixed frequency data
% given the initialization. 
%       Inputs:
%               y           (array)
%               A           (array)
%               C           (array)
%               Sig_eta     (array)
%               Sig_eps     (array)
%               x_init      (array)
%               P_init      (array) 
%               cfg         (struct)
%       
%       Outputs:
%               x           (array)
%               P           (array)
%               loglik      (array)
%               ll          (scalar)
%               v           (array)
%               K           (array)
%               invF        (array)
%=========================================================================%

if ~isfield(cfg, 't0')
    cfg.t0 = 0;
end

[ny,nx] = size(C);
T       = numel(y(1,:));
%% KALMAN FILTER
%init
x    = zeros(nx,T);        P    = zeros(nx,nx,T); 
v    = zeros(ny,T);        invF = zeros(ny,ny,T);     
K    = zeros(nx,ny,T);     ll   = zeros(T,1);

% if size(P_init,1) == 5
%     wait_a_minute=0;
% end

x(:,1)      = x_init;
P(:,:,1)    = P_init;
rescale_prediction_error_covariance0=false; badly_conditioned_F=false;
rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
kalman_tol = 1e-10;
if size(y,1)> sum(cfg.measurement)
    y = y(cfg.measurement,:);
end
%% main recursion
for t =1:T-1
%     
%     % forecast   
    W = eye(ny);
    missing = isnan(y(:,t+1)); y(isnan(y(:,t+1)),t+1) = 0;
    ix = find(~missing);
    W = W(~missing,:);
    c      = W * C;
    v(~missing,t) = y(~missing,t+1) - c*x(:,t);
    F      = c*P(:,:,t)*c' + W * Sig_eps * W';
    FF(~missing,~missing,t) = F;

    if rcond(F)<kalman_tol
        sig=sqrt(diag(F));
        if any(diag(F)<kalman_tol) || rcond(F./(sig*sig'))<kalman_tol
            badly_conditioned_F = true;
        else
            rescale_prediction_error_covariance=true;
        end
        %                 badly_conditioned_F = true;
    end
    
    if badly_conditioned_F
        x = Inf;
        P = Inf;
        loglik = -1e9;
        ll = Inf;
        v = Inf;
        K = Inf;
        invF = Inf;
        if ~all(abs(F(:))<kalman_tol)
            % Use univariate filter.
            return
        else
            % Pathological case, discard draw
            return
        end
    else
        F_singular = false;
        if rescale_prediction_error_covariance
            logdetF = log(det(F./(sig*sig')))+2*sum(log(sig));
            invF = inv(F./(sig*sig'))./(sig*sig');
            rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
        else
            logdetF = log(det(F));
            if isnan(rcond(F)) || rcond(F)<eps
                invF = sqrt(-1); % if rcond is very low penalize ll
            else
                invF = inv(F);
            end
%             invF = inv(F);
%             if any(~isfinite(F))
%                 wth =0;
%             end
        end
    end
    K(:,~missing,t) = A*P(:,:,t)*c'*invF;
    
    %loglikelihood - prediction error decomposition
    if t> cfg.t0 %m^3 % it takes some observations for correctly estimate the state
        ll(t)       =  - (1/2)*sum(~missing)*log(2*pi) - (1/2)*logdetF - (1/2)*(v(~missing,t)'*invF*v(~missing,t));
    end
    if t<T %update
        x(:,t+1)    = A*x(:,t) + K(:,~missing,t)*v(~missing,t);
        P(:,:,t+1)  = A*P(:,:,t)*(A - K(:,~missing,t)*c)' + Sig_eta;
    end
end
% APA'- APc'K + Sig_eta

% penalize if imaginary
if sum(abs(imag(ll))) > 0
    ll = real(ll) - 1e66;
end

loglik =  sum(ll) / (T-cfg.t0);
