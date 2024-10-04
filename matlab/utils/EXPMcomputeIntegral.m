function [SIG_ETA, Wh] = EXPMcomputeIntegral(A, Omega, C, h, flowmeas, varargin)
% EXPMCOMPUTEINTEGRAL computes integrals involving the exponential matrix
% Computation of integrals involving the exponential matrix applied to a
% state space with stocks and flows.
%-------------------------------------------------------------------------%
% Inputs:
%   A       : (nxn)     drift matrix
%   Omega   : (nxn)     diffusion matrix (B*B')
%   C       : (mxn)     loading matrix
%   h       : (scalar)  time-step (e.g., 1/12, 1/4, ... )
%   flowmeas: (logical) if the state-space includes flow variables (true)
%
% References:   Van Loan (1978); 
%               Harvey, Stock (1993); 
%               Christensen, Neri, Parra-Alvarez (2022)
%
%=========================================================================%
% The continuous-time state space is:
%   dx(t) = A x(t) dt + B dW(t);        %transition equation
%   y(t)  = C x(t);                     %measurement equation
%
%-------------------------------------------------------------------------%
% flowmeas = true; =======================================================%
% With y(t) being FLOW VARIABLES, the EXACT DISCRETE MODEL of the
% transition equation, for time-step h fixed for all t>0 is:
% x(t)  = exp(A * h) * x(t-h) + etas(t);
% yf(t) = C * W(h)   * x(t-h) + etaf(t);
%
% with measurement equation
% y(t) = [0 I] [x(t)' yf(t)']';
%
% The variance-covariance matrix of the reduced form innovations
% SIG_ETA = [   Sig_etas        Sig_etas_etaf;
%               Sig_etas_etaf'  Sig_etaf];
%
% Autocovariances of order greater than zero are 0;
%
%-------------------------------------------------------------------------%
% flowmeas = false; ======================================================%
% With only STOCK VARIABLES in the model, the EXACT DISCRETE MODEL of the
% transition equation, for fixed delta is:
% x(t) = exp(A * h) * x(t-h) + etas(t)
%
% with measurement equation
% y(t) = C * x(t);
% The variance covariance matrix of the reduced form innvoations
% SIG_ETA = [Sig_etas];
% Autocovariances of order greater than zero are 0;
%
%=========================================================================%
%-------------------------------------------------------------------------%

if ~isempty(varargin) % is it for a ABCD minimal representation?
    ABCD = logical(varargin{1});
else
    ABCD = false;
end

Om = full(Omega); % rename variable for pretty
% [m, n] = size(C);

[n, ~] = size(A);
O = zeros(n);       % empty matrix
I = eye(n);         % identity matrix

PP = [     -A       I       O       O;
            O      -A       Om      O;
            O       O       A'      I;
            O       O       O       O];


% [     F1      G1       H1      K1;
%        O      F2       G2      H2;
%        O       O       F3      G3;
%        O       O       O       F4];

P = expm(PP * h);

F3 = P(2*n+1:3*n, 2*n+1:3*n);
H2 = P(n+1:2*n, 3*n+1:4*n);
G2 = P(n+1:2*n, 2*n+1:3*n);
G3 = P(2*n+1:3*n,3*n+1:4*n);
K1 = P(1:n, 3*n+1:4*n);

% int_{0}^h e^(As) ds
Wh = G3';

% int_0^h e^(A(s)) * Om * e^(A'(s)) ds
Sig_etas  = F3'*G2;
% G2 = int_0^h e^(-A(h-s)) Om e^(A'(s)) ds 
% F3 = e^(A(h))' 

%flowmeas=false-----------------------------------------------------------%
SIG_ETA = Sig_etas;
%-------------------------------------------------------------------------%
if flowmeas
    % int_0^h int_0^s e^(Ar) * Om * e^(A'r) dr ds
    Sig_etaf = (F3'*K1) + (F3'*K1)';
    

    
    % int_0^h int_0^s e^(A(s-r)) * Om * e^(A'r) dr ds
    Sig_etas_etaf = F3'*H2;
    % H_2 = int_0^h e^(-A(h-s)) * Om * exp(A'(s)
    
    if ABCD
        Sig_etaf        = Sig_etaf; 
        Sig_etas_etaf   = Sig_etas_etaf;
    else
        Sig_etaf        = Sig_etaf/(h^2);
        Sig_etas_etaf   = Sig_etas_etaf/h;
    end
%flowmeas=true------------------------------------------------------------%    
    SIG_ETA =   [   Sig_etas        Sig_etas_etaf;
                    Sig_etas_etaf'  Sig_etaf];
%-------------------------------------------------------------------------%
end

end %eof