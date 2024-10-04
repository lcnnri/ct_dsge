%% Steady state
z_ss    = 0;
n_ss    = (1-S.alpha)/(S.gamma*(1-S.alpha*(S.delta + S.eta)/(S.rho + S.delta + S.eta)));
k_ss    = (S.alpha/(S.rho + S.delta + S.eta))^(1/(1-S.alpha))*n_ss;
c_ss    = k_ss^S.alpha * n_ss^(1-S.alpha) - (S.delta + S.eta)*k_ss;
y_ss    = exp(z_ss)*k_ss^S.alpha*n_ss^(1-S.alpha);
ybar    = [c_ss; n_ss;y_ss];
y = [data{1}.ct; data{1}.nt; data{1}.yt];
N       = length(data{1}.yt); % number of observations
[A, b, C, D, R, B, H] = prepareFilter(theta,y(logical(cfg.measurement),:),cfg); % system matrices
[ny, nx] = size(C);
%% add measurement error if there is any
if cfg.measErr 
    for n=1:reps
        rng(666)
        y = [data{n}.ct; data{n}.nt; data{n}.yt];
        yy = y;
        err = randn(numel(measurement_error),N); % simulate errors
        y = y-ybar;
        for k=1:N
            y(logical(cfg.measurement),k) = y(logical(cfg.measurement),k) + measurement_error.*err(:,k);
        end
            y = y + ybar;
    end
end