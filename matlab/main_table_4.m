%% MAIN_TABLE_4.m ====================================================
% Replicating Table 4 from Christensen, Neri, Parra-Alvarez (2024)
% Estimation of continuous-time linear DSGE models from discrete measurements

% %% Housekeeping
% clear all;
% clc;
% close all;
% 
% %% Configuration for Table 4
serial_id = 'replica Table 4'
Table = '4';
Panel = 'a';
if exist('overwrite_panel','var')
    Panel = overwrite_panel;
end
disp(['Table ' Table ', panel '  Panel])

%%


%% Load calibrated parameters and empirical data
addpath('utils/', 'Linearization/AutoDiff/','Linearization/phact/', 'model_loglin/');

table_configurations()

[theta0] = load_calibrated_params();
configuration_settings();
load_empirical_data();

%% ML Estimation
perform_ml_estimation();

%% Save and Display Results
display_results();
