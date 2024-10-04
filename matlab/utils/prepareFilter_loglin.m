% PREPAREFILTER construct system matrices of the state space
% representation.
% Inputs:
%   theta: (array) vector of structural parameters
%   y    : (array) vector of data
%   cfg  : (struct)cofiguration settings
%                   containing information cfg.measErr regarding
%                   measurement errors
%
% References: Christensen, Neri, Parra-Alvarez (2020)
%=========================================================================%
function [A, b, C, D, B, Sig_eps,eigviol,cfg]= prepareFilter_loglin(theta,y,cfg)

gamma = theta(1);
alpha = theta(2);
delta = theta(3);
rho   = theta(4);
sigz  = theta(5);
rhoz  = theta(6);
eta   = theta(7);
sigk  = theta(8);

try
    [G1,F,impact,st_st,~, A, B, C] = model_solution_loglin(theta); 
catch
    return
end
% check
smallno = 1e-15;
A(abs(A) < smallno) = 0;
B(abs(B) < smallno) = 0;
C(abs(C) < smallno) = 0;
eigA = eig(A);
eigviol = 0;
if max(real(eigA)) >= 0
    eigviol = 1;
end

cfg.ss = st_st;


css = cfg.ss(1);
kss = cfg.ss(2);
zss = cfg.ss(3);
nss = cfg.ss(4);
yss = cfg.ss(5);
xss = cfg.ss(6);
rss = cfg.ss(7);

%% The model

D = [css; nss; yss; rss]; % steady state

b = []; % model the demean transition
if ~isfield(cfg, 'info')
    cfg.info.used_in_estimation = logical([1 1 0 0]');
end
Sig_eps = eps*diag(cfg.info.used_in_estimation); %eps*diag(cfg.info.used_in_estimation);
if ~isfield(cfg, 'measErr')
    cfg.measErr = false;
end
if cfg.measErr
    sig = theta(8+1:end).^2; % variance of measurement error
%     if cfg.flows
%         sig = cfg.h * sig;
%     end
    append_sig_eps = double(diag(cfg.info.used_in_estimation));
    append_sig_eps([find(append_sig_eps)]) = sig;
    Sig_eps = Sig_eps + append_sig_eps;
    
else
    Sig_eps = Sig_eps + zeros(numel(cfg.info.used_in_estimation)); % matrix of zeros
end



 
 end
