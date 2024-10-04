if ~(exist('plot_data')==1)
    plot_data=false;
end
%% load data =============================================================%
data_path = 'data/dataset.xlsx';
t1 = datetime(1960,01,01); % first quarter 1960
t2 = datetime(2019,12,01); % last  quarter 2019

% monthly data -----------------------------------------------------------%
%-------------------------------------------------------------------------%
% Read the data using readtable
opts = readtable(data_path, 'Sheet', 'monthly_data');

% Extract the string months and convert them to datetime format
% Assuming the date is in the first column
opts.strmonths = opts.(opts.Properties.VariableNames{1}); % Access first column by name
opts.strmonths = datetime(opts.strmonths, 'InputFormat', 'dd/MM/yyyy');

% Generate the monthly date range
datesm = t1:calmonths(1):t2;

% Find the start and end indices within the date range
opts.logicstartm = (t1 == opts.strmonths);  % Logical index for start date
opts_tmp.ixstartm = find(opts.logicstartm, 1);  % Find the first occurrence of the start date

opts.logicendm = (t2 == opts.strmonths);    % Logical index for end date
opts_tmp.ixendm = find(opts.logicendm, 1);  % Find the first occurrence of the end date

% Handle cases where the date is not found
if isempty(opts_tmp.ixstartm)
    error('Start date not found in the dataset.');
end

if isempty(opts_tmp.ixendm)
    error('End date not found in the dataset.');
end

% Extract the data matrix within the date range
M = opts{opts_tmp.ixstartm:opts_tmp.ixendm, {'C', 'r'}};

% Add NaN columns - as not all variables are measured monthly
M = [M(:,1), nan(length(M),1), nan(length(M),1), M(:,2)];
clear opts

% quarterly data ---------------------------------------------------------%
%-------------------------------------------------------------------------%
% Read the data using readtable
opts = readtable(data_path, 'Sheet', 'quarterly_data');

% Extract the string quarters and convert them to datetime format
% Assuming the date is in the first column
opts.strquarters = opts.(opts.Properties.VariableNames{1}); % Access first column by name
opts.strquarters = datetime(opts.strquarters, 'InputFormat', 'dd/MM/yyyy');

% Generate the quarterly date range
datesq = t1:calmonths(3):t2; % Quarterly frequency (3 months)

% Find the start and end indices within the date range
opts.logicstartq = (t1 == opts.strquarters);  % Logical index for start date
opts_tmp.ixstartq = find(opts.logicstartq, 1);  % Find the first occurrence of the start date

opts.logicendq = (datesq(end) == opts.strquarters);  % Logical index for end date
opts_tmp.ixendq = find(opts.logicendq, 1);  % Find the first occurrence of the end date

% Handle cases where the date is not found
if isempty(opts_tmp.ixstartq)
    error('Start date not found in the dataset.');
end

if isempty(opts_tmp.ixendq)
    error('End date not found in the dataset.');
end

% Extract the data matrix within the date range
Q = opts{opts_tmp.ixstartq:opts_tmp.ixendq, {'C', 'N', 'Y', 'r'}}; 

Q(:,2) = ... Hours worked divided by the average number in a month
    Q(:,2)/(4.34524 * 7 * 24);  

% Shift dates by 2 months 
datesq = datesq + calmonths(2);

clear opts opts_tmp

%% 

iC=1; iN=2; iY=3; ir=4; % location of measurements
% plot data --------------------------------------------------------------%
if plot_data % plot data
    load_palette
    tiledlayout(3,1)
    f1=nexttile;
%     f1=figure(1) % plot Consumption growth
    plot(datesm(2:end),12 * 100 * (log(M(2:end,iC)./M(1:end-1,iC))-theta0(7)*(1/12)),'-.', ...
        'Color',colors(3,:),'LineWidth',1.5)
    hold on
    plot(datesq(2:end), 4 *100* (log(Q(2:end,iC)./Q(1:end-1,iC))-theta0(7)*0.25), ...
        'o', 'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
    grid on;
    recessionplot;
    hold off;
        legend({'Monthly', 'Quarterly'},'Interpreter','latex','Location','northeast')
    ylabel('Consumption','Interpreter','latex')
    setmyfig_out1(f1)
    f2=nexttile;
    plot(datesq(2:end), 100*Q(2:end,iN), ...
               'o', 'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
    recessionplot;
    ylim(100*[min(Q(2:end,iN)),max(Q(2:end,iN))])
    grid on;
    ylabel('Hours worked','Interpreter','latex')
    setmyfig_out1(f2)
    f3=nexttile;
    plot(datesq(2:end), 100*log(Q(2:end,ir)),'o', 'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
    recessionplot;
    ylim([min(100*log(Q(2:end,ir))),max(100*log(Q(2:end,ir)))])
    grid on;
    hold off;
    ylabel('Real interest rate','Interpreter','latex')
    xlabel('Month','Interpreter','latex')
    setmyfig_out1(f3)
print('output/figures/quarterly_measurements.eps', '-depsc')
end

% change the bootstrap type in case of 
if any(what_freq(cfg.measurement)=='m')
    % nonparametric bootstrap not implemented for mixed frequency measurements
    % therefore change it to 'Wild bootstrap'
    bs_type = 'wild'; 
    ix = zeros(length(M),1);
    for q=1:numel(Q(:,1))
        ix = ix + (datesm == datesq(q))';
    end
else
    ix = ones(length(Q),1);
end

ix=logical(ix); ixq=logical(ixq); ixm=logical(ixm);
DATA = nan(length(ix),size(which_meas,1));
if any(ixm) && any(what_freq(cfg.measurement)=='m')
    DATA(:,ixm) = M(:,ixm);
end

if ~any(what_freq(cfg.measurement)=='m')
    DATA(ix,ixm) = Q(:,ixm);
end
if any(ixq)
    DATA(ix,ixq) = Q(:,ixq);
end

