function [T,Pss] =SolveRiccati(tspan, P0,A,Z,H,Q)
[T Pss] = ode45(@(t,X)mRiccati(t, X, A, Z, H, Q), tspan, P0(:));
end

function dXdt = mRiccati(t, X, F, Z, H, Q)
X = reshape(X, size(F)); %Convert from "n^2"-by-1 to "n"-by-"n"
dXdt = F*X + X*F' + Q;
% dXdt = A*X*A' - (A*X*Z'/(Z*X*Z'+H))*Z*X*A' + Q; %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end
