function [G1, F, impact, st_st, eu, A, B, C] = model_solution_loglin(theta)

%% =======================================================================%
%                          Structural parameters                          %
%-------------------------------------------------------------------------%
global gamma alpha delta rho sigz rhoz eta sigk n_v n_g n_p nEErrors nVars;

gamma    = theta(1);
alpha    = theta(2);
delta    = theta(3);
rho      = theta(4);
sigz     = theta(5);
rhoz     = theta(6);
eta      = theta(7);
sigk     = theta(8);
n_v      = 1; % Number of controls
n_g      = 2; % Number of states (both endogenous and exogenous)
n_p      = 1; % Number of static equations
nEErrors = 1; % Number of expectational errors
nVars    = n_v + n_g + n_p;

%% =======================================================================%
%                       Deterministic steady state                        %
%-------------------------------------------------------------------------%
global varsSS


zss   = 0;
nss   = (1-alpha)/(gamma*(1-alpha*(delta+eta)/(rho+delta+eta)));
kss   = ((alpha*exp(zss))/(rho+delta+eta))^(1/(1-alpha))*nss;
css   = exp(zss)*kss^alpha*nss^(1-alpha) - (delta+eta)*kss;
yss   = exp(zss)*kss^alpha*nss^(1-alpha);
rss   = alpha*exp(zss)*kss^(alpha-1)*nss^(1-alpha) - delta;
 
%-------------------------------------------------------------------------%
%                       Store steady state values                         %
%-------------------------------------------------------------------------%
varsSS      = zeros(nVars,1);
varsSS(1,1) = log(css);
varsSS(2,1) = log(kss);
varsSS(3,1) = zss;
varsSS(4,1) = log(nss);
st_st       = varsSS;
st_st(5,1)  = log(yss);
st_st(6,1)  = log(yss-css);
st_st(7,1)  = log(rss);

%% =======================================================================%
%                Linearization of Equilibrium Conditions                  %
%-------------------------------------------------------------------------%
%                 Prepare automatic differentiation                       %  
%-------------------------------------------------------------------------%
vars = zeros(2*nVars+nEErrors+2,1);
vars = myAD(vars);

%-------------------------------------------------------------------------%
%                      Evaluate derivatives                               %
%-------------------------------------------------------------------------%
derivativesIntermediate = equilibrium_conditions_loglin(vars);

%-------------------------------------------------------------------------%
%                   Extract out derivative values                         %
%-------------------------------------------------------------------------%
derivs = getderivs(derivativesIntermediate);

%-------------------------------------------------------------------------%
%                       Unpack derivatives                                %
%-------------------------------------------------------------------------%
mVarsDerivs    = derivs(:,1:nVars);
mVarsDotDerivs = derivs(:,nVars+1:2*nVars);
mEErrorsDerivs = derivs(:,2*nVars+1:2*nVars+nEErrors);
mShocksDerivs  = derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+2);

%-------------------------------------------------------------------------%
% Rename derivatives to match notation in paper: g0*E(dot X) = c + g1*X   %
%-------------------------------------------------------------------------%
g0  = mVarsDotDerivs;
g1  = -mVarsDerivs;
c   = sparse(nVars,1);
psi = -mShocksDerivs;
pi  = -mEErrorsDerivs;

%-------------------------------------------------------------------------%
%            Solve out static constraints of the linear model             %
%-------------------------------------------------------------------------%
[~,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi);
n_g_red = n_g;
from_spline = speye(n_g_red + n_v);
Bn = full(inv_state_red*from_spline);

%% =======================================================================%
%                       Solve the linear model                            %
%                       dx = G1*x*dt + impact*dZ                          %
%-------------------------------------------------------------------------%
% [G1,~,impact,q,a,b,z,eu]= gensysct(full(g0),full(g1),full(c),full(psi),full(pi));
% [~,~,~,~,F] = schur_solver(g0,g1,c,psi,pi,1,1,1);
[G1,~,impact,eu,F] = schur_solver(g0,g1,c,psi,pi,1,1,1);

if eu(1)==0 || eu(2) == 0
    fprint('Model has no solution or solution is not unique')
    fprint('Algorithm stops here! Check the model parameters')
    return;
end

%% =======================================================================%
%                     State-space representation                          %
%                    dx = A*x*dt + B*dZ, x = [k,z]                        %
%                        m = Cx, m = [c,n,y]                              %
%-------------------------------------------------------------------------%
A = [G1(2,1)*F(1,1) + G1(2,2), G1(2,1)*F(1,2) + G1(2,3);                   % k
     G1(3,1)*F(1,1) + G1(3,2), G1(3,1)*F(1,2) + G1(3,3)];                  % z
B = impact(2:end,1:end);

phi_nk = Bn(end,1)*F(1,1) + Bn(end,2);
phi_nz = Bn(end,1)*F(1,2) + Bn(end,3);
rr=(alpha*kss^(alpha-1)*nss^(1-alpha))/rss;
C = [F(1,1), F(1,2);                                                        % c
     phi_nk, phi_nz ;                                                       % n
     alpha + (1-alpha)*(phi_nk),1+(1-alpha)*(phi_nz);                       % y
     ((alpha-1)+(1-alpha)*phi_nk)*rr, (1 + (1-alpha)*phi_nz)*rr];           % r


