% MAIN SCRIPT:
% Script to replicate Tables 4, 5, and 7 of the empirical application in 
% Section 6 of “Estimation of Continuous-time Linear DSGE Models from 
% Discrete-time Measurements” by Bent Jesper Christensen, Luca Neri, and 
% Juan Carlos Parra-Alvarez (2024), Journal of Econometrics, forthcoming. 
%
% This version: September, 2024. 

%% run Table 4
clear all; close all; clc;
main_table_4;
% pause

clear all; close all; clc;
overwrite_panel = 'b';
main_table_4;
% pause

%% run Table 5
clear all; close all; clc;
main_table_5;
% pause

clear all; close all; clc;
overwrite_panel = 'b';
main_table_5;
% pause

%% run Table 7
clear all; close all; clc;
main_table_7;
% pause

clear all; close all; clc;
overwrite_panel = 'b';
main_table_7;
% pause

%% check that the identification holds
check_identification_condition;