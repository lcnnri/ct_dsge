function [ID,OC] = isIdentified(theta,cfg,y,varargin)
%ISDENTIFIED Function that checks identification conditions of the model
%presented in 
% Christensen, Neri, Parra-Alvarez (2024). 
% It applies Proposition 2 of that paper.
%   Dependencies: 
%       KalmanFilter_mixed_freq_ABCD
%       prepareFilter_loglin
%   
%   Inputs: 
%       theta       (array)     vector of structural parameters
%       cfg         (struct)    structure of configurations settings
%       y           (array)     nYxT matrix of data
%       varargin    (char)      Choose which part of the proposition 2 of paper
%                               Christensen, Neri, Parra-Alvarez (2024):   
%                               allowed char. are 'a', 'b' or 'c'
%   Outputs:
%       ID          (bool)      true: theta is identified; false: otherwise
%       OC          (bool)      true: the order condition is met; false:
%                                       otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

derivs_theta = get_dtheta(theta,cfg,y);

if isempty(varargin)
    prop = 'a';
else
    prop = varargin{1};
    if ~ischar(prop)
        warning("identification matrix Delta can be of type 'a', 'b', or 'c'");
        prop = 'a';
    end
end

[info, A, B, C, Sig_ups, Ah, Bh, Ch, Dh, Kh, Sig_epsilon, Sig_nu] = get_matrices(theta,cfg,y);
if ~check_minimality(Ah,Kh,Ch)
    warning('The AKC system is not minimal')
end

nYs = info.nys;
nY  = info.ny;
nX  = info.nx;
dimtheta=numel(theta);
switch prop %Proposition 2 
    case 'a'  % (a)
        n_eq = nYs*nX + 2*nX*nY + nY*(nY+1)/2;
        rank_cond = dimtheta + nX^2;

        DeltaTheta = ...
            [derivs_theta.Cs;
            derivs_theta.Ah;
            derivs_theta.Kh;
            derivs_theta.Ch;
            derivs_theta.Sig_nu];

        DeltaT = ...
            [-kron(eye(nX), C(info.idx_stocks,:));
            kron(Ah',eye(nX)) - kron(eye(nX),Ah);
            kron(Kh',eye(nX));
            - kron(eye(nX),Ch);
            zeros(nY*(nY + 1)/2, nX^2)];

        Delta_plus = ...
            [DeltaTheta, DeltaT];
        Delta = Delta_plus;
    case 'b' % (b)
        rank_cond = dimtheta + nX^2;
        n_eq = 2*nX*nY + nY*(nY+1)/2;

        DeltaTheta = ...
            [derivs_theta.A;
            derivs_theta.C;
            derivs_theta.Kh;
            derivs_theta.Sig_nu];

        DeltaT = ...
            [kron(A',eye(nX))- kron(eye(nX),A);
            -kron(eye(nX),C);
            kron(Kh',eye(nX));
            zeros(nY*(nY+1)/2,nX^2)];

        Delta_a = ...
            [DeltaTheta,  DeltaT];
        Delta = Delta_a;

    case 'c' % (c)
        if rank(C(info.idx_stocks,:)) < size(C,2)
            error('The rank assumption for Cs is not satisfied')
        end
        mW = size(B,2);
        % sig_ups is known by assumption. adjust theta
        num_msmt_err_to_fix=numel(find(diag(Sig_ups)));
        dimtheta=dimtheta-num_msmt_err_to_fix;
        theta=theta(1:dimtheta);
        cfg.par_ident=cfg.par_ident(1:dimtheta);
        %
        rank_cond = dimtheta + nX^2 + mW^2;
        n_eq = nX*nY + (nX-mW)*mW + nY*(nY+1)/2;

        DeltaTheta = ...
            [derivs_theta.A(:,1:dimtheta);
            derivs_theta.B(:,1:dimtheta);
            derivs_theta.C(:,1:dimtheta)];

        DeltaT = ...
            [kron(A',eye(nX))-kron(eye(nX),A);
            kron(full(B)',eye(nX));
            -kron(eye(nX),C)];

        DeltaU = ...
            [zeros(nX^2,mW^2);
            kron(eye(mW), full(B));
            zeros(nX*nY, mW^2)];

        Delta_c = ...
            [DeltaTheta, DeltaT, DeltaU];
        Delta = Delta_c;
end

%% Include Calibrated parameters
if any(~logical(cfg.par_ident))
    r = ~logical(cfg.par_ident);
    nR = sum(r);
    ix = find(r);
    Restr = zeros(nR,numel(theta));
    for i=1:nR
        Restr(i,ix(i)) = 1;
    end

    Delr  = [Restr, zeros(nR,size(Delta,2)-dimtheta)];
    Delta = [Delr; Delta];
    n_eq = n_eq + nR;
end
%% Rank and order restrictions
rdel = rank(full(Delta));
rankerror=false;
if rdel <  rank_cond % rank condition
    warning('\theta is not locally identified, rank condition is not met')
    rankerror = true;
end

ordererror = false;
if dimtheta > n_eq % order condition
    warning('\theta is not locally identified, order condition is not met')
    ordererror = true;
end

id = ~any(rankerror+ordererror); % (id =1) system is identified, not identified otherwise

ID=id; OC=~ordererror;
end % eof




function [info, A, B, C, Sig_ups, Ah, Bh, Ch, Dh, Kh, Sig_epsilon, Sig_nu] = ...
                                                get_matrices(theta,cfg,y,varargin)
%GET_MATRICES gets the matrices of the system, both CT and DT
% Inputs:
%       theta (doubl) vector of parameters
%       cfg (doub) configuration settings)
%       y (doub) data
%       vectorize (bool) get vectorized system matrices or not

if ~isempty(varargin)
    vectorize = varargin{1};
else
    vectorize = false;
end

[A, b, C, d, B, Sig_ups, eigviol] = ... get the system matrices
    prepareFilter_loglin(theta,[],cfg); % continuous-time system matrices
C = C(cfg.measurement,:);
Sig_ups = Sig_ups(cfg.measurement,cfg.measurement);

% EDSSM
[ny,nx] = size(C);
h  = cfg.h;
Ah = expm(A*h);
if strcmpi(cfg.model,'EM-SSR')
    Ah = eye(nx) + A*h;
end

if strcmpi(cfg.model,'S-SSR')
    Wh = Ah;
    stock_measurements = true(ny,1);
    flow_measurements = [];
elseif strcmpi(cfg.model,'EM-SSR')
    Ah = eye(nx) + A*h;
    Wh = Ah;
    stock_measurements = true(ny,1);
    flow_measurements = [];
else
    Wh = A\(Ah-eye(nx))/h;
    stock_measurements = strcmpi(cfg.info.sampling(cfg.measurement), 'stock');
    flow_measurements  = strcmpi(cfg.info.sampling(cfg.measurement), 'flows');
end
Ch = C;
Ch(flow_measurements,:)  = C(flow_measurements,:)*Wh;
Ch(stock_measurements,:) = C(stock_measurements,:);

% C matrix for stocks
Chs = Ch(stock_measurements,:);
nys = size(Chs,1);

if strcmpi(cfg.model,'EM-SSR')
    Sig_eta = full(B*B')*h;
else
%  Covariance matrices
Om = full(B*B');
O  = zeros(nx);       % empty matrix
I  = eye(nx);         % identity matrix

PP = [     -A       I       O       O;
    O      -A       Om      O;
    O       O       A'      I;
    O       O       O       O];

P = expm(PP * h);
F3 = P(2*nx+1:3*nx, 2*nx+1:3*nx);
H2 = P(nx+1:2*nx, 3*nx+1:4*nx);
G2 = P(nx+1:2*nx, 2*nx+1:3*nx);
G3 = P(2*nx+1:3*nx,3*nx+1:4*nx);
K1 = P(1:nx, 3*nx+1:4*nx);
Wh = G3';

Sig_eta  = F3'*G2;
end

Sig_eps  = blkdiag(Sig_eta,Sig_ups);


% ABCD
Ah  = Ah;
Bh  = [eye(nx) zeros(nx,ny)];
Ch  = Ch;
Dh  = [Ch, eye(ny)];

Sig_epsilon = blkdiag(Sig_eta,Sig_ups);

x_0 = zeros(size(Ah,1),1);
P_0 = reshape((eye(nx^2) - kron(Ah,Ah))\Sig_eta(:),nx,nx);

%% Kalman filter iteration
[nll, ~, ~,~, ~, K, invF] = ...
    KalmanFilter_mixed_freq_ABCD(y, Ah, Ch, Sig_eta, Dh*Sig_epsilon*Dh', x_0, P_0,cfg);

Cs = C(stock_measurements,:);

if vectorize
    A = vec(A);
    B = vec(full(B));
    Cs = vec(Cs);
    C = vec(C);
    Sig_ups = vech(Sig_ups);
    Sig_epsilon= vech(Sig_epsilon);
    Ah = vec(Ah);
    Bh = vec(Bh);
    Ch = vec(Ch);
    Dh = vec(Dh);
    Kh = vec(K(:,:,end-1));
    Sig_nu = vech(inv(invF));
else
    Kh = K(:,:,end-1);
    Sig_nu = inv(invF);
end

info.nll = nll;
info.nx = nx;
info.ny = ny;
info.nys = nys;
info.idx_stocks = stock_measurements;
info.Cs = Cs;
end


function [dtheta] = get_dtheta(theta,cfg,y)
%GET_DTHETA gets the derivatives of the system matrices
%   inputs
%       theta (doub)    vector of parameters
%       cfg (struct)    configurartion settings
%       y (doub)        data
%   outputs
%       derivative wrt theta
[info, A, B, C, Sig_ups, Ah, Bh, Ch, Dh, Kh, Sig_epsilon, Sig_nu] = get_matrices(theta,cfg,y,true);
Cs = info.Cs;
cfg.par_ident=true(numel(theta),1);
dtheta= [];
tol = sqrt(eps);

dtheta.A = []; dtheta.B = []; dtheta.C = [];
dtheta.Sig_ups = []; dtheta.Sig_epsilon = []; dtheta.Sig_nu = [];
dtheta.Ah = []; dtheta.Bh = []; dtheta.Ch = [];  dtheta.Dh = [];
dtheta.Kh = []; dtheta.Cs = [];

for i=1:numel(theta)
    itheta = theta;
    itheta(i) = itheta(i) + tol;
    [info, iA, iB, iC, iSig_ups, iAh, iBh, iCh, iDh, iKh, iSig_epsilon,iSig_nu] = get_matrices(itheta,cfg,y,true);
    iCs = info.Cs;

    dtheta.A = [dtheta.A, (A-iA)/tol];
    dtheta.B = [dtheta.B, (B-iB)/tol];
    dtheta.C = [dtheta.C, (C-iC)/tol];
    dtheta.Cs = [dtheta.Cs, (Cs-iCs)/tol];

    dtheta.Ah = [dtheta.Ah, (Ah-iAh)/tol];
    dtheta.Bh = [dtheta.Bh, (Bh-iBh)/tol];
    dtheta.Ch = [dtheta.Ch, (Ch-iCh)/tol];
    dtheta.Dh = [dtheta.Dh, (Dh-iDh)/tol];
    dtheta.Kh = [dtheta.Kh, (Kh-iKh)/tol];

    dtheta.Sig_ups = [dtheta.Sig_ups, (Sig_ups-iSig_ups)/tol];
    dtheta.Sig_epsilon = [dtheta.Sig_epsilon, (Sig_epsilon-iSig_epsilon)/tol];
    dtheta.Sig_nu = [dtheta.Sig_nu, (Sig_nu-iSig_nu)/tol];
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