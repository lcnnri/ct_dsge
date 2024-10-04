%function [Ah, Sig_eta, C, D, x_init, T,nx,ny,y] = DiscretizeSystem(y,cfg,A,B,C,D,Sig_eps,h,x_init)
function [Ah, Sig_eta, C, D, x_init, T,nx,ny,y,cfg,ABCD] = DiscretizeSystem_wp4(y,cfg,A,B,C,D,Sig_eps,h)
% DISCRETIZESYSTEM discretize the linear continuous-time time-invariant
% system 
% dx(t) = A * x(t) dt + B * dW(t);  W(t) ~ N(0,t*I)
%  y(t) = C * x(t) + meas_err.      meas_err ~ iid(0, dt*Sig_eps)
% The program allows:
%   o- mixed frequency of the measurements (monthly, quarterly)
%   o- mixed sampling of the measurements  (stocks, flows)
%   o- exact discretization (EDM)
%   o- Euler-Maruyama 'naive' discretization (EM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg, 'em')
    cfg.em = false;
end
if cfg.em
    cfg.flows=false;
end
if ~isfield(cfg, 'growth_rates') % use log growth rates?
    cfg.growth_rates = false;
end
if ~isfield(cfg, 'mixed_sampling');
    cfg.mixed_sampling=false;
end
if ~isfield(cfg, 'minimal_representation')
    cfg.minimal_representation=false;
end
cfg.Act = A;
cfg.Ct = C;
% dimensions of the system
T       = numel(y(1,:));
[ny,nx] = size(C);
n_latent_states =size(A,1);
n_obs_states = sum(cfg.info.used_in_estimation);
%% Construct Discrete Model
if cfg.em % Euler-Maruyama approximation
    [Ah, Sig_eta,cfg] = Sys_em(A,B,h,cfg);
else % Exact Discrete Representation
    [Ah, Sig_eta,BB,cfg] =  Sys_edm(A,B,C,h,cfg); % stocks
    if cfg.flows  % if MEASUREMENTS ARE FLOWS
        [Ah,  C,Sig_eta,cfg] = Sys_edm_flows(A,Ah,BB,C,Sig_eta,h,cfg);
    end
end


if cfg.growth_rates
    cfg.n_latent_states = n_latent_states;
    [Ah, C, Sig_eta] = add_differences_in_transition_eq(Ah,C,Sig_eta,cfg);

end

mixed_frequency= any(cfg.freq(logical(cfg.measurement))=='m') && any(cfg.freq(logical(cfg.measurement))=='q');
hf_vars=(cfg.freq=='m').*cfg.info.used_in_estimation; lf_vars=(cfg.freq=='q').*cfg.info.used_in_estimation;

if cfg.flows && mixed_frequency 
    nx = size(Ah,1);
    flow_lf_vars = lf_vars.*strcmpi(cfg.info.sampling,'flow');
    Ah = blkdiag(Ah, zeros(2*sum(flow_lf_vars)));
    ix_lf_vars = find(flow_lf_vars);
    Ah(nx+1,n_latent_states+ix_lf_vars) = 1;
    Ah(nx+2,n_latent_states+ix_lf_vars+1) = 1;
    
    ny = size(C,1);
    C = [C zeros(ny,2*sum(flow_lf_vars))];
    C(logical(flow_lf_vars),end-1:end) = 1;

    Sig_eta = blkdiag(Sig_eta, zeros(2*sum(flow_lf_vars)));

end    

x_init = zeros(size(Ah,1),1);
scale = cfg.scale; Sy=cfg.Sy;
if mixed_frequency
    slow_freq_flow =(cfg.freq(cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'flow')) == 'q');
    scaletmp = scale(Sy,Sy);
    scaletmp(slow_freq_flow,slow_freq_flow) =  scaletmp(slow_freq_flow,slow_freq_flow)/3;
    scale(Sy,Sy) = scaletmp;
    clear scaletmp
%     C(slow_freq_flow,:) = C(slow_freq_flow,:)/3;
%     Sig_eps(slow_freq_flow, slow_freq_flow)/(3^2)
end
Sig_eps = scale(Sy,Sy).^2 * Sig_eps(Sy, Sy);  % resize measurement error covariance matrix
C       = scale(Sy,:) * C; 

cfg.Sig_eps=Sig_eps;
cfg.scale  = scale;


% C(C~=0)= 1;


%
[ABCD] = getABCD(Ah,B,C,Sig_eta,Sig_eps,cfg);

% [CheckCO,minns,minSYS] = get_minimal_state_representation(ABCD, 0,false);
%     
% ct = ss(A,eye(2),4*cfg.C_original,zeros(2),0);
% dt = ss(ABCD.A,ABCD.B,ABCD.C, ABCD.D,0.25);
% t = 0:0.25:240;
% [Yct,Tct]  = step(ct,t);
% [Ydt, Tdt] = step(dt,t);
% close all
% plot(t,Yct(:,1)); hold on; plot(t,Ydt(:,1)); hold off;
% plot(t,Yct(:,2)); hold on; plot(t,Ydt(:,2)); hold off;
%

if cfg.minimal_representation && ~mixed_frequency
    slow_freq_flow =(cfg.freq(cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'flow')) == 'q');
    cfg.slow_freq_flow=slow_freq_flow;
   [minABCD] = get_structural_minimal_state_representation(A,cfg.Ct,Sig_eta,Sig_eps,cfg,mixed_frequency);
   ABCD=minABCD;
end
end

%%

function [Ah, Sig_eta,cfg] = Sys_em(A,B,h,cfg)
% EULER-MARUYAMA discretization
    nx = size(A,1);
    Ah = (eye(nx) + A*h);  % I + A*h
    Sig_eta = B*B'*h;      % h * BB
    cfg.A2 = A; %
    cfg.B = B*sqrt(h); %
    cfg.H = cfg.B; % impact matrix of structural shocks
end

function [Ah, Sig_eta,BB,cfg] = Sys_edm(A,B,C,h,cfg)
% EXACT DISCRETE MODEL for stock data
    Ah = expm(A*h); % exp(A*h)
    BB = B*B';
    [Sig_eta, ~] = EXPMcomputeIntegral(A, BB, C, h, cfg.flows);
    cfg.A2 = expm(A*h/2); % exp(A*h)
    cfg.B = full(B);
%     cfg.H = (1/sqrt(h))*HComputeIntegral(A,B,h);
    cfg.H = sqrt(h)*cfg.A2*cfg.B;
end

function [Ah, C,Sig_eta,cfg] = Sys_edm_flows(A,Ah,BB,C,Sig_eta, h,cfg)
% EXACT DISCRETE MODEL for flow data
    [ny,nx] = size(C);
    ix_meas = cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'flow');
    n_meas  = sum(ix_meas); %cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'flow') );
    Wh = A\(Ah-eye(nx))/h;
    
    Ah = [Ah, zeros(nx,n_meas);
            C(ix_meas,:)*Wh, zeros(n_meas)]; 
    Sig_eta = [Sig_eta(1:nx,1:nx), ...
               Sig_eta(1:nx,nx+1:end)*C(ix_meas,:)';
               C(ix_meas,:)*Sig_eta(nx+1:end,1:nx), ... 
               C(ix_meas,:)*Sig_eta(nx+1:end,nx+1:end)*C(ix_meas,:)']; % h
    cfg.C_original = C(cfg.info.used_in_estimation,:);

    cfg.A2 = [expm(A*h/2), zeros(nx,n_meas);
                C(ix_meas,:)*(A\(expm(A*(h/2))-eye(nx))), zeros(n_meas)];
%     cfg.H = cfg.A2 * [cfg.B; cfg.B];

    I = diag(strcmpi(cfg.info.sampling,'flow')); I = I(:,strcmpi(cfg.info.sampling,'flow') & cfg.info.used_in_estimation);
    if cfg.mixed_sampling
        tmp = zeros(ny,nx); 
        tmp(strcmpi(cfg.info.sampling,'stock'),:) = C(strcmpi(cfg.info.sampling,'stock'),:);
    else
        tmp = zeros(ny,nx);
%         I = eye(ny);
    end
    C = [tmp I];

end

function [A, C, Sigma] = add_differences_in_transition_eq(A,C,Sigma,cfg)
% ADD_DIFFERENCES_IN_TRANSITION_EQ add differenced measurements in the
% transition equation. This is used when the growth rates of the
% measurements are used in the observation equations (instead of the
% detrended measurements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_latent_states = cfg.n_latent_states; % # of latent states
    nx = size(A,1);
    which_measurement_take_difference = ... take differences of which measurements
        cfg.info.used_in_estimation.*cfg.info.take_differences;
    n_add_states = ... # of states to add to the transition eqn
        sum(which_measurement_take_difference); 

    tmp     = diag(cfg.info.take_differences);
    tmp_flo = tmp(find(which_measurement_take_difference & strcmpi(cfg.info.sampling,'flow')), ...
                          cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'flow'));
    tmp_sto = tmp(find(which_measurement_take_difference & strcmpi(cfg.info.sampling,'stock')), ...
        cfg.info.used_in_estimation & strcmpi(cfg.info.sampling,'stock'));
%     if cfg.flows
    A = blkdiag(A,zeros(n_add_states)); % expand dimensions of A
    A(nx+1:nx+size(tmp_flo,1),n_latent_states+1:n_latent_states+size(tmp_flo,2)) = tmp_flo;
    A(end-size(tmp_sto,1)+1:end,1:n_latent_states) =  ...
        C(which_measurement_take_difference & strcmpi(cfg.info.sampling,'stock'),1:n_latent_states);
%     else
%         A = blkdiag(A,zeros(n_add_states)); % expand dimensions of A
%         A(nx+1:nx+size(tmp_flo,1),n_latent_states+1:n_latent_states+size(tmp_flo,2)) = tmp_flo;
%         A(end-size(tmp_sto,1)+1:end,1:n_latent_states) =  ...
%             C(which_measurement_take_difference & strcmpi(cfg.info.sampling,'stock'),1:n_latent_states);
%     end

    tmpC = diag(which_measurement_take_difference);
    C = [C -tmpC(:,logical(which_measurement_take_difference))];
    Sigma = blkdiag(Sigma,zeros(n_add_states));
end

function [sys] = getABCD(a,b,c,Sig_eta,Sig_eps,cfg)
% GETABCD construct the ABCD system given the system matrices
% construct:
%   x(t+1) = A * x(t) + B * e(t+1)
%   y(t+1) = C * x(t) + D * e(t+1)
%  with E[e(t)*e(t)'] = SIGMA
% 
% starting from:
%   x(t+1) = a * x(t) + eta(t+1)
%   y(t)   = c * x(t) + eps(t)
%  with E[eta(t)*eta(t)'] = Sig_eta
%  and  E[eps(t)*eps(t)'] = Sig_eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=cfg.h;
[ny,nx] = size(c);
j=1;
for i=1:nx
    if any(Sig_eta(i,:))
    S(j,:) = zeros(1,nx);
    S(j,i) = 1;
    j=j+1;
    end
end
j=1;
U=[];
for i=1:ny
    if any(Sig_eps(i,:))
    U(j,:) = zeros(1,ny);
    U(j,i) = 1;
    j=j+1;
    end
end
num_eta = size(S*Sig_eta*S',1);
if isempty(U)
    num_error=0;
    USig_epsU = [];
else
    num_error = size(U*Sig_eps*U',1);
    USig_epsU = U*Sig_eps*U';
end
num_stoch = num_eta + num_error;
SIGMA = blkdiag(S*Sig_eta*S', USig_epsU);


I=eye(nx);
B1 = [];
for i=1:nx;
    if any(Sig_eta(i,:))
        B1 = [B1 I(:,i)];
    end
end

I=eye(ny);
D=[];
for i=1:ny
    if any(Sig_eps(i,:))
        D = [D I(:,i)]; 
    end
end
D = [c*B1 D];

sys.A = a;
sys.B = [B1 zeros(nx,num_error)];
sys.C = c * a;
sys.D = D; 
sys.Sig_eta = sys.B*SIGMA*sys.B';
sys.Sig_eps = sys.D*SIGMA*sys.D';
sys.SIGMA = SIGMA;
end


function [H] = HComputeIntegral(A,B,h)
% HCOMPUTEINTEGRAL compute the H matrix from the integral
% H(theta,h) := int_0^h exp(A*s) B ds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nx,nw]=size(B);

Onx = zeros(nx); Onxnw = zeros(nx,nw); I = eye(nx);
Q = full(B*B');
C = [-A'    I   Onx Onxnw;
    Onx    -A' Q   Onxnw;
    Onx    Onx A   full(B);
    Onx    Onx Onx Onxnw];
eCh = expm(C*h);
H = eCh(2*nx+1:3*nx,3*nx+1:end);
end



function  [minABCD] = get_structural_minimal_state_representation(A,C,Sig_eta,Sig_eps,cfg,mixed_frequency);
h = cfg.h;
mx=size(A,1);
measurements = cfg.info.used_in_estimation;
ix_flow  = strcmpi(cfg.info.sampling,'flow');
ix_stock = strcmpi(cfg.info.sampling,'stock');
Syf = measurements & ix_flow;
Sys = measurements & ix_stock;
Sxf = eye(mx);
errors = measurements & cfg.info.with_meas_error;
ny = size([C(Sys,:);C(Syf,:)],1);

C = cfg.scale*C;
Ah = expm(A*h);
minABCD.A = Ah;

minABCD.B = [eye(mx) zeros(mx,sum(Syf)) zeros(mx,mixed_frequency*2*sum(cfg.slow_freq_flow)) zeros(mx,sum(errors))];
minABCD.C = [
C(Sys,:)*Ah
C(Syf,:) * Sxf *( A\ (Ah-eye(mx)))/h];

% if mixed_frequency
%     n_slow = sum(cfg.slow_freq_flow); slow_var = cfg.slow_freq_flow;
%     cslow = C(Syf,:) * Sxf *( A\ (Ah-eye(mx)))/h;
%     tmp = zeros(ny,2*n_slow);
%     tmpflow = tmp(Syf,:);
%     tmpflow(n_slow,:) = 
%     cslow = cslow(slow_var,:); 
%     minABCD.C = [minABCD.C ];
%     Ah = [Ah;
%           zeros(n_slow,mx)
% end


c = [C(Sys,:); C(Syf,:)*Sxf];
I=eye(ny);
D=[];
for i=1:ny
    if any(Sig_eps(i,:))
        D = [D I(:,i)]; 
    end
end

minABCD.D = [blkdiag(C(Sys,:)*eye(mx), C(Syf,Syf)*eye(sum(Syf))) , D];
minABCD.SIGMA = blkdiag(Sig_eta,Sig_eps);
minABCD.Sig_eta = Sig_eta;
minABCD.Sig_eps = Sig_eps;

end