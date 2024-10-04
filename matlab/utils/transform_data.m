% TRANSFORM_DATA     transform the data to stationary
% Inputs:
%   Yobs    (double)    array, cointains data
%   cfg     (struct)    contains configuration settings.
%
% Outputs:
%   y     (double)    array, transformed data
%   cfg   (struct)    contains configuration settings 
%% =======================================================================%
function [y, cfg] = transform_data(Yobs,cfg)

h = cfg.h;
S = double(diag(cfg.info.used_in_estimation)); % select measurements to use in estimation
S = logical(diag(S));                                  %boolean: selection of measurements

scale = [100/h 0   0     0;    % c
         0     100 0     0;    % n
         0     0   100/h 0;    % y
         0     0   0     100]; % r


Y     = log(Yobs);
missing = isnan(Y);

% 

if cfg.growth_rates
    y = [Y(1,2:end) - Y(1,1:end-1); % Delta log C
        Y(2,2:end);                % log N
        Y(3,2:end) - Y(3,1:end-1); % Delta log Y
        Y(4,2:end)];               % log r
    missing = missing(:,2:end);
    % Schorfheide uses
    % Y(2,t) = C * X(:,t) + beta (parameter capturing scale to fit data) -
    % beta is calibrated in the empirical application
else
    order=1; % order=1, linear detrending; order=2, quadratic trend etc...
    y = nan(size(Y));
    y(1,~missing(1,:)) = detrend(Y(1,~missing(1,:)),order);
%   X = [ones(1,239); (1:239)];
%   beta =   (X*X')\X*Y(1,2:end)'
%   trend_par = beta(2); DeltaMethod= 1/0.25; 
%   sig = Y(1,2:end)' - X'*beta
    y(2,:) = Y(2,:);
    y(3,~missing(3,:)) = detrend(Y(3,~missing(3,:)),order);
    y(4,:) = Y(4,:);
end

for i=1:numel(y(:,1))
    y(i,:) = scale(i,i)*y(i,:);
    y(i,:) = y(i,:) - mean(y(i,~missing(i,:)),2);
end
y = y(S,:);

cfg.scale   = scale;
cfg.Sy      = S;
end
