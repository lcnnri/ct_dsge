function [CheckCO,minns,minSYS] = get_minimal_state_representation(SYS, derivs_flag,verbose)
% Derives and checks the minimal state representation
% Let x = A*x(-1) + B*u  and  y = C*x(-1) + D*u be a linear state space
% system, then this function computes the following representation
%     xmin = minA*xmin(-1) + minB*u  and  and  y=minC*xmin(-1) + minD*u
%
% -------------------------------------------------------------------------
% INPUTS
% SYS           [structure]
%     with the following necessary fields:
%     A:        [nspred by nspred] in DR order
%                 Transition matrix for all state variables
%     B:        [nspred by exo_nbr] in DR order
%                 Transition matrix mapping shocks today to states today
%     C:        [varobs_nbr by nspred] in DR order
%                 Measurement matrix linking control/observable variables to states
%     D:        [varobs_nbr by exo_nbr] in DR order
%                 Measurement matrix mapping shocks today to controls/observables today
%     and optional fields:
%     dA:       [nspred by nspred by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix A
%     dB:       [nspred by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix B
%     dC:       [varobs_nbr by nspred by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix C
%     dD:       [varobs_nbr by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix D
% derivs_flag   [scalar]
%                 (optional) indicator whether to output parameter derivatives
% verbose       [boolean]
%                 (optional) indicator whether to print text on screen
% -------------------------------------------------------------------------
% OUTPUTS
%   CheckCO:    [scalar]
%                 equals to 1 if minimal state representation is found
%   minns:      [scalar]
%                 length of minimal state vector
%   SYS         [structure]
%               with the following fields:
%     minA:     [minns by minns] in DR-order
%                 transition matrix A for evolution of minimal state vector
%     minB:     [minns by exo_nbr] in DR-order
%                 transition matrix B for evolution of minimal state vector
%     minC:     [varobs_nbr by minns] in DR-order
%                 measurement matrix C for evolution of controls, depending on minimal state vector only
%     minD:     [varobs_nbr by minns] in DR-order
%                 measurement matrix D for evolution of controls, depending on minimal state vector only
%     dminA:    [minns by minns by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix minA
%     dminB:    [minns by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix minB
%     dminC:    [varobs_nbr by minns by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix minC
%     dminD:    [varobs_nbr by u_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix minD
% -------------------------------------------------------------------------
% This function is called by
%   * get_identification_jacobians.m (previously getJJ.m)
% -------------------------------------------------------------------------
% This function calls
%   * check_minimality (embedded)
%   * minrealold (embedded)
% =========================================================================
% Copyright (C) 2019-2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
% =========================================================================
if nargin == 1
    derivs_flag = 0;
end
if nargin == 2
    verbose = true;
end
realsmall = 1e-7;
[nspred,exo_nbr] = size(SYS.B);
varobs_nbr = size(SYS.C,1);

% Check controllability and observability conditions for full state vector
CheckCO = check_minimality(SYS.A,SYS.B,SYS.C);
if CheckCO == 1 % If model is already minimal, we are finished
    minns = nspred;
    minSYS = SYS;
else
    %Model is not minimal    
    try
        minreal_flag = 1;
        % In future we will use SLICOT TB01PD.f mex file [to do @wmutschl], currently use workaround
        [mSYS,U] = minrealold(SYS,realsmall);
        minns = size(mSYS.A,1);
        CheckCO = check_minimality(mSYS.A,mSYS.B,mSYS.C);
        if CheckCO
            minSYS.A = mSYS.A;
            minSYS.B = mSYS.B;
            minSYS.C = mSYS.C;
            minSYS.D = mSYS.D;
            if derivs_flag
                totparam_nbr = size(SYS.dA,3);
                minSYS.dA = zeros(minns,minns,totparam_nbr);
                minSYS.dB = zeros(minns,exo_nbr,totparam_nbr);
                minSYS.dC = zeros(varobs_nbr,minns,totparam_nbr);
                % Note that orthogonal matrix U is such that (U*dA*U',U*dB,dC*U') is a Kalman decomposition of (dA,dB,dC)                    % 
                for jp=1:totparam_nbr
                    dA_tmp = U*SYS.dA(:,:,jp)*U';
                    dB_tmp = U*SYS.dB(:,:,jp);
                    dC_tmp = SYS.dC(:,:,jp)*U';
                    minSYS.dA(:,:,jp) = dA_tmp(1:minns,1:minns);
                    minSYS.dB(:,:,jp) = dB_tmp(1:minns,:);
                    minSYS.dC(:,:,jp) = dC_tmp(:,1:minns);
                end
                minSYS.dD = SYS.dD;
            end
        else
            minSYS = []; 
            minns = [];
            return;
        end
    catch
        minreal_flag = 0; % if something went wrong use below procedure
    end

    if minreal_flag == 0
        if verbose
            fprintf('Use naive brute-force approach to find minimal state space system.\n These computations may be inaccurate/wrong as ''minreal'' is not available, please raise an issue on GitLab or the forum\n')
        end
        % create indices for unnecessary states
        exogstateindex = find(abs(sum(SYS.A,1))>realsmall);
        minns = length(exogstateindex);
        % remove unnecessary states from solution matrices
        A_2 = SYS.A(exogstateindex,exogstateindex);
        B_2 = SYS.B(exogstateindex,:);
        C_2 = SYS.C(:,exogstateindex);
        D   = SYS.D;
        ind_A2 = exogstateindex;
        % minimal representation
        minSYS.A = A_2;
        minSYS.B = B_2;
        minSYS.C = C_2;
        minSYS.D = D;

        % Check controllability and observability conditions
        CheckCO = check_minimality(minSYS.A,minSYS.B,minSYS.C);

        if CheckCO ~=1
            % do brute-force search
            j=1;
            while (CheckCO==0 && j<minns)
                combis = nchoosek(1:minns,j);
                i=1;
                while i <= size(combis,1)
                    ind_A2 = exogstateindex;
                    ind_A2(combis(j,:)) = [];
                    % remove unnecessary states from solution matrices
                    A_2 = SYS.A(ind_A2,ind_A2);
                    B_2 = SYS.B(ind_A2,:);
                    C_2 = SYS.C(:,ind_A2);
                    D = SYS.D;
                    % minimal representation
                    minSYS.A = A_2;
                    minSYS.B = B_2;
                    minSYS.C = C_2;
                    minSYS.D = D;
                    % Check controllability and observability conditions
                    CheckCO = check_minimality(minSYS.A,minSYS.B,minSYS.C);
                    if CheckCO == 1
                        minns = length(ind_A2);
                        break;
                    end
                    i=i+1;
                end
                j=j+1;
            end
        end
        if derivs_flag
            minSYS.dA = SYS.dA(ind_A2,ind_A2,:);
            minSYS.dB = SYS.dB(ind_A2,:,:);
            minSYS.dC = SYS.dC(:,ind_A2,:);
            minSYS.dD = SYS.dD;
        end
    end
end

function CheckCO = check_minimality(a,b,c)
%% This function computes the controllability and the observability matrices of the ABCD system and checks if the system is minimal
%
% Inputs: Solution matrices A,B,C of ABCD representation of a DSGE model
% Outputs: CheckCO: equals 1, if both rank conditions for observability and controllability are fullfilled
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
CheckCO = CheckC&CheckO;        % equals 1 if minimal state
end % check_minimality end

function [mSYS,U] = minrealold(SYS,tol)
    % This is a temporary replacement for minreal, will be replaced by a mex file from SLICOT TB01PD.f soon
    a = SYS.A;
    b = SYS.B;
    c = SYS.C;
    [ns,nu] = size(b);
    [am,bm,cm,U,k] = ControllabilityStaircaseRosenbrock(a,b,c,tol);
    kk = sum(k);
    nu = ns - kk;
    nn = nu;
    am = am(nu+1:ns,nu+1:ns);
    bm = bm(nu+1:ns,:);
    cm = cm(:,nu+1:ns);
    ns = ns - nu;
    if ns
        [am,bm,cm,t,k] = ObservabilityStaircaseRosenbrock(am,bm,cm,tol);
        kk = sum(k);
        nu = ns - kk;
        nn = nn + nu;
        am = am(nu+1:ns,nu+1:ns);
        bm = bm(nu+1:ns,:);
        cm = cm(:,nu+1:ns);
    end
    mSYS.A = am;
    mSYS.B = bm;
    mSYS.C = cm;
    mSYS.D = SYS.D;
end


function [abar,bbar,cbar,t,k] = ObservabilityStaircaseRosenbrock(a,b,c,tol)
    %Observability staircase form
    [aa,bb,cc,t,k] = ControllabilityStaircaseRosenbrock(a',c',b',tol);
    abar = aa'; bbar = cc'; cbar = bb';
end

function [abar,bbar,cbar,t,k] = ControllabilityStaircaseRosenbrock(a, b, c, tol)
    % Controllability staircase algorithm of Rosenbrock, 1968
    [ra,ca] = size(a);
    [rb,cb] = size(b);
    ptjn1 = eye(ra);
    ajn1 = a;
    bjn1 = b;
    rojn1 = cb;
    deltajn1 = 0;
    sigmajn1 = ra;
    k = zeros(1,ra);
    if nargin == 3
        tol = ra*norm(a,1)*eps;
    end
    for jj = 1 : ra
        [uj,sj,vj] = svd(bjn1);
        [rsj,csj] = size(sj);
        %p = flip(eye(rsj),2);
        p = eye(rsj);
        p = p(:,end:-1:1);
        p = permute(p,[2 1 3:ndims(eye(rsj))]);
        uj = uj*p;
        bb = uj'*bjn1;
        roj = rank(bb,tol);
        [rbb,cbb] = size(bb);
        sigmaj = rbb - roj;
        sigmajn1 = sigmaj;
        k(jj) = roj;
        if roj == 0, break, end
        if sigmaj == 0, break, end
        abxy = uj' * ajn1 * uj;
        aj   = abxy(1:sigmaj,1:sigmaj);
        bj   = abxy(1:sigmaj,sigmaj+1:sigmaj+roj);
        ajn1 = aj;
        bjn1 = bj;
        [ruj,cuj] = size(uj);
        ptj = ptjn1 * ...
              [uj zeros(ruj,deltajn1); ...
               zeros(deltajn1,cuj) eye(deltajn1)];
        ptjn1 = ptj;
        deltaj = deltajn1 + roj;
        deltajn1 = deltaj;
    end
    t = ptjn1';
    abar = t * a * t';
    bbar = t * b;
    cbar = c * t';
end

end % Main function end
