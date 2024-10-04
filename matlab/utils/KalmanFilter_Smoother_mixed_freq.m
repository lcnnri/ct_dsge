function [x, P, loglik,ll, struct_shock, v, K, invF, alpha] = KalmanFilter_Smoother(y, A, C, Sig_eta, Sig_eps, x_init, P_init,cfg)
% KALMANFILTER_SMOOTHER runs the Kalman smoother and recovers structural
% shocks from a countinuous-time structural model
% Inputs
%   y           :       (ny-by-T)  measurements
%   A           :       (nx-by-nx) autoregressive matrix of states
%   C           :       (ny-by-nx) loadings
%   Sig_eta     :       (nx-by-nx) covariance matrix reduced-form innovations
%   Sig_eps     :       (ny-by-ny) covariance matrix measurement errors
%   x_init      :       (nx-by-1)  vector of initial values of states
%   P_init      :       (nx-by-nx) prediction covariance matrix init values
%   cfg         :       (struct)   configurations
%
% Outputs
%   x           :       (nx-by-T) filtered states
%   P           :       (nx-by-nx) MSE of states
%   loglik      :       (scalar) log-likelihood
%   struct_shock:       (struct) containing matrices of smoothed structural shocks 
%   v           :       (nx-by-T) smoothed errors
%   K           :       (ny-by-nx-by-T) Kalman gain
%   invF        :       (ny-by-ny-by-T) inverse of covariance matrix
%   alpha       :       (nx-by-T) smoothed states
%
% Author: ln 2021, Aarhus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
L    = zeros(nx,nx,T);
% if size(P_init,1) == 5
%     wait_a_minute=0;
% end

x(:,1)      = x_init;
P(:,:,1)    = P_init;
rescale_prediction_error_covariance0=false; badly_conditioned_F=false;
rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
kalman_tol = 1e-10;
%% main recursion
for t =1:T
%     
%     % forecast   
    W = eye(ny);
    missing = isnan(y(:,t)); y(isnan(y(:,t)),t) = 0;
    ix = find(~missing);
    W = W(~missing,:);
    c      = W * C;
    v(~missing,t) = W * y(:,t) - c*x(:,t);
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
            invF(~missing,~missing,t) = inv(F./(sig*sig'))./(sig*sig');
            rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
        else
            logdetF = log(det(F));
            if isnan(rcond(F)) || rcond(F)<eps
                invF(~missing,~missing,t) = sqrt(-1); % if rcond is very low penalize ll
            else
                invF(~missing,~missing,t) = inv(F);
            end
%             invF = inv(F);
%             if any(~isfinite(F))
%                 wth =0;
%             end
        end
    end
    K(:,~missing,t) = A*P(:,:,t)*c'*invF(~missing,~missing,t);
    L(:,:,t)        = A - K(:,~missing,t)*c;

    %loglikelihood - prediction error decomposition
    if t> cfg.t0 %m^3 % it takes some observations for correctly estimate the state
        ll(t)       =  - (1/2)*sum(~missing)*log(2*pi) - (1/2)*logdetF - (1/2)*(v(~missing,t)'*invF(~missing,~missing,t)*v(~missing,t));
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISTURBANCE SMOOTHER ==================================================%
% Durbin & Koopman 2012
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init
r(:,T)  = zeros(nx,1);        N(:,:,T)  = zeros(nx,nx);
e       = zeros(ny,T);               n  = zeros(nx,T);
D       = zeros(ny,ny,T); 
Vare    = zeros(ny,ny,T);          Varn = zeros(nx,nx,T);
alpha   = zeros(nx,T); 
alpha(:,T) = x(:,T); Ps = P;

for t=T:-1:2
    r(:,t-1)   = C'*invF(:,:,t)*v(:,t) + L(:,:,t)' * r(:,t);
    N(:,:,t-1) = C'*invF(:,:,t)*C + L(:,:,t)'*N(:,:,t)*L(:,:,t);
    alpha(:,t) = x(:,t) + P(:,:,t)*r(:,t-1);
    V(:,:,t)   = P(:,:,t) - P(:,:,t)*N(:,:,t-1)*P(:,:,t);
end

V(abs(V)<1e-15) = 0;

%
% for t=T:-1:2  
%     u(:,t)   = invF(:,:,t)*v(:,t) - K(:,:,t)'*r(:,t);
%     
%      %disturbances
%     e(:,t)   = Sig_eps*u(:,t);
%     n(:,t)   = Sig_eta*r(:,t);
%     
%     D(:,:,t) = invF(:,:,t) + K(:,:,t)'*N(:,:,t)*K(:,:,t);
%     
%     %variance of disturbances given data
%     Vare(:,:,t) = Sig_eps - Sig_eps*D(:,:,t)*Sig_eps; Vare(:,:,t) = 0.5*(Vare(:,:,t)+Vare(:,:,t)');
%     Varn(:,:,t) = Sig_eta - Sig_eta*N(:,:,t)*Sig_eta; Varn(:,:,t) = 0.5*(Varn(:,:,t)+Varn(:,:,t)');
%     
%     N(:,:,t-1)  = C'*D(:,:,t)*C ...
%                     + A'*N(:,:,t)*A ...
%                     - C'*K(:,:,t)'*N(:,:,t)*A ...
%                     - A'*N(:,:,t)*K(:,:,t)*C;
%         %auxiliary disturbance
%     r(:,t-1)    = C'*u(:,t) + A'*r(:,t);
%         %state smoothed
%     alpha(:,t)  = x(:,t) + P(:,:,t)*r(:,t-1);
% 
%     %MSE smoothed
% %     [L,U] = OuterLU(P(:,:,t));
% %     L=sparse(L); U=sparse(U);
% %     invPt = U\(L\eye(nx));
%     invPt = inv(P(:,:,t));
%     Gt          = P(:,:,t-1)*A'*invPt;
%     Ps(:,:,t)   = P(:,:,t-1) + Gt*(Ps(:,:,t) - P(:,:,t))*Gt'; Ps(:,:,t) = 0.5*(Ps(:,:,t)+Ps(:,:,t)');
% end
% r0    = C'*u(:,1) + A'*r(:,1);
% alpha(:,1)  = x(:,1) + P(:,:,1)*r0;
% n(:,1)   = Sig_eta*r(:,1);

eta(:,1) = zeros(nx,1);
for t=2:T
    eta(:,t) = alpha(:,t) - A*alpha(:,t-1);
end
if isfield(cfg, 'h')
    h = cfg.h;
else
    h = 0.25; % quarterly data;
end
if isfield(cfg, 'A2')
%     Ah2= expm(cfg.A*h/2);
    Ah2 = cfg.A2;
    Act = cfg.Act;
else
    Ah2 = A;
end
if isfield(cfg, 'B')
    B = cfg.B;
else
    error('<KalmanFilter_Smoother>: Matrix B not provided for recovery of structural shocks')
end
if cfg.flows 
    eta2 = eta*(h/2); 
else
    eta2 = eta;
end
% nct = size(Act,1);
nrow= size(A,1); 
ncol= size(A,2);
if cfg.flows ;%|| (~cfg.flows && cfg.growth_rates)
   ncol = size(B,2);
    
    Ah2 = [expm(Act*(h/2)), zeros(2,2);
          (cfg.C_original)*(Act\(expm(Act*(h/2))-eye(size(Act,2))))/(h/2),zeros(2,2)]; %B = B/sqrt(h/2);
    
    nrow = size(Ah2,1);
%     nrow=2;
    if cfg.growth_rates
        Ah2temp = A;
        Ah2temp(1:nrow,1:nrow) = Ah2;
        Ah2 = Ah2temp;
        clear Ah2temp;
    end
%     Ah2(3:end,3:end) = 0; B(3:end,:) = 0; eta(3:end,:) = 0; %H%
    Ah2 = Ah2(1:2,1:2); B = B(1:2,:); eta2 = eta2(1:2,:); %Hs
%     Ah2 = Ah2(3:4,1:2); B = B(1:2,:); eta = eta(3:4,:);%Hf
end
    
[struct_shock.epsil_edr, struct_shock.H]= StructuralShocks_EDR(Ah2, B, eta2,h,cfg.flows,[]);
epsil0 = struct_shock.epsil_edr(:,1);
nct = size(Act,1);
[struct_shock.epsil_em, struct_shock.H_em]= StructuralShocks_EM(alpha(1:nct,:),B(1:nct,:),Act,h,epsil0);
if cfg.flows 
    struct_shock.epsil_edr = struct_shock.epsil_edr/(h/2);
    struct_shock.epsil_em = struct_shock.epsil_em;
end
struct_shock.eta = eta;
end

function [epsil, H] = StructuralShocks_EDR(A, B, eta,h,is_flows,varargin)
% STRUCTURALSHOCKS back out structural shocks
    H = A*B;
    [nx,q] = size(H);
    if q<nx
        invHH = pinv(H'*H);
        sol = invHH*H';        
    elseif q==nx
        sol = eye(nx)/H;
    end
    T = size(eta,2);
    epsil = zeros(q,T);
    for t=1:T % recover structural shocks
          epsil(:,t)=  sol * eta(:,t) ;  
    end
    epsil = epsil/sqrt(h);
    H = H*sqrt(h);

%     if  ~isempty(varargin) && q<nx % FGLS
% %          H = H/sqrt(h);
% %          epsil = epsil*sqrt(h);
%          Sigeta = varargin{1};
%          H = H/sqrt(h);
%          uhat = eta - H*epsil;
%          Sigu =  uhat*uhat'; %H*H' - Sigeta;
%          invSigu = pinv(Sigu);
%          invHH = inv(H'*invSigu*H);
%          sol = invHH*H'*invSigu;
% 
%          epsil = zeros(q,T);
%          for t=1:T % recover structural shocks
%              epsil(:,t)=  sol * eta(:,t) ;
%          end
%         epsil = epsil/sqrt(h);
%         H = H*sqrt(h);
%     end
    
end

function [epsil,H] = StructuralShocks_EM(alpha,B,A,h,epsil0)
% STRUCTURALSHOCKS_EM_APPROX back out structural shocks using the
% Euler-Maruyama approximation
    T = size(alpha,2);
    [nx,q] = size(B);
    if nx>q
        invBB = inv(B'*B);
        sol = invBB *B';
    elseif nx==q
        sol = eye(nx)/(B);
    end
    w = zeros(q,T);
    for t=T-1:-1:1 
        w(:,t) = w(:,t+1) ... 
                 - sol * (alpha(:,t+1)-alpha(:,t)) ...
                 + sol * A*alpha(:,t)*h;
    end
    for t=T:-1:1
        w(:,t) = w(:,t) - w(:,1);
    end
    epsil = zeros(q,T);
    epsil(:,1) = epsil0;
    for t=1:T-1
       epsil(:,t+1) = w(:,t+1) - w(:,t);
    end
    epsil = epsil/sqrt(h);
    H = B*sqrt(h);
%     H = H*sqrt(h);
end