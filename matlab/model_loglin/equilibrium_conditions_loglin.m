function vResidual = equilibrium_conditions_loglin(vars)
%-------------------------------------------------------------------------%
% Inputs:  (1) vars: vector which contains the variables in the system,   %
%                    the time derivatives of those variables, the         %
%                    expectational errors, and the shocks.                %                    
%                                                                         % 
% Outputs: (1) vResduals: residuals of equilibrium conditions, evaluated  %
%                         at vars                                         %
%% =======================================================================%
%                             Housekeeping                                %
%-------------------------------------------------------------------------%
%                       Declare global variables                          %
%-------------------------------------------------------------------------%
global gamma alpha delta rho sigz rhoz eta sigk nEErrors nVars varsSS

%-------------------------------------------------------------------------%
%             Unpack vars (vars are in log-deviations from ss)            % 
%-------------------------------------------------------------------------%
% 1. Controls                                                             % 
%-------------------------------------------------------------------------%
lc = vars(1) + varsSS(1);
%-------------------------------------------------------------------------%
% 2. States (endo + exo)                                                  % 
%-------------------------------------------------------------------------%
lk = vars(2) + varsSS(2);
z  = vars(3) + varsSS(3);
%-------------------------------------------------------------------------%
% 3. Static variables                                                     % 
%-------------------------------------------------------------------------%
ln = vars(4) + varsSS(4);
%-------------------------------------------------------------------------%
% 4. Time derivatives (controls+states)                                   % 
%-------------------------------------------------------------------------%
lc_dot = vars(nVars+1);
lk_dot = vars(nVars+2);
z_dot  = vars(nVars+3);
%-------------------------------------------------------------------------%
% 5. Expectational errors                                                 % 
%-------------------------------------------------------------------------%
lc_EErrors = vars(2*nVars+1);
%-------------------------------------------------------------------------%
% 6. Shocks                                                               % 
%-------------------------------------------------------------------------%
lk_shock  = vars(2*nVars+nEErrors+1);
z_shock   = vars(2*nVars+nEErrors+2);

%% =======================================================================%
%                             Model                                       %
%-------------------------------------------------------------------------%
%                   Define static equations                               % 
%-------------------------------------------------------------------------%
n = (((1-alpha)*exp(z)*exp(lk)^alpha)/(gamma*exp(lc)))^(1/alpha); 

%-------------------------------------------------------------------------%
%                   Equilibrium dynamics                                  % 
%-------------------------------------------------------------------------%
% Euler Equation
lc_Residual = lc_dot + lc_EErrors - (alpha*exp(z)*exp(lk)^(alpha-1)*exp(ln)^(1-alpha) - rho -  delta - eta);

% Capital dynamics 
% lk_Residual = lk_dot - ((exp(z)*exp(lk)^(alpha)*exp(ln)^(1-alpha) - exp(lc))/exp(lk) - (delta+eta) - 0.5*sigk^2) - sigk*lk_shock;
lk_Residual = lk_dot - ((exp(z)*exp(lk)^(alpha)*exp(ln)^(1-alpha)  - exp(lc))/exp(lk) - (delta+eta) - 0.5*sigk^2) - sigk*lk_shock;

% Labor supply
ln_Residual = log(n) - ln;

% Law of motion for aggregate shocks
z_Residual = z_dot + rhoz*z - sigz*z_shock;

vResidual = [lc_Residual;lk_Residual;z_Residual;ln_Residual];