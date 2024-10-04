function theta0 = load_calibrated_params()
    % Load calibrated parameters
    S = paramstructbase();
    theta0 = [S.gamma; S.alpha; S.delta; S.rho; S.sigz; S.rhoz; S.eta; S.sigk];
end