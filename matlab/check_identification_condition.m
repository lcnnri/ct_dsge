    
clear all;
clc;
close all;

fname = ... input name of MLE results from output/workspaces
    'Table 4 - a.mat'

choose_model = 1;  
% 1 - F-SSR
% 2 - S-SSR
% 3 - EM-SSR

choose_condition = 'a';
% 'a' - Proposition 2 (a)
% 'b' - Proposition 2 (b)
% 'c' - Proposition 2 (c) 



%% Load calibrated parameters and empirical data
addpath('utils/', 'Linearization/AutoDiff/','Linearization/phact/', 'model_loglin/');
load(['output/workspaces/' fname])

cfg = res{choose_model}.cfg; % extract configurations
theta_hat = res{choose_model}.theta_hat; % extract parameter estimate
[y, cfg] = ... transform the data
    transform_data(DATA',cfg);

[idFlag, OCflag] = ... identification, cf. Proposition 2
        isIdentified(theta_hat,cfg,y,choose_condition); 


if idFlag
    disp(['The model in ' fname(1:end-4)   ' is identified!'])
else
    error('The model is NOT identified!')
end