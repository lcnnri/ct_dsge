plot_data=true;
which_meas = ...  data analogues of which model measurements to use
    {'c' 1 ;  % consumption
     'n' 1 ;  % hours worked
     'y' 0 ;  % output
     'r' 0};  % interest rate
what_freq =  ['q';'q';'q'; 'q'];
for ny = 1:size(which_meas,1)
    cfg.measurement(ny,1) = logical(which_meas{ny,2});
    ixm(ny,1) = strcmpi(what_freq(ny),'m');
    ixq(ny,1) = strcmpi(what_freq(ny),'q');
end
load_empirical_data;

load('quart_freq_c_n_1000_bsse_alpha5_20230914.mat')
qX = res{1}.s_states;

load('mixed_freq_c_n_1000_bsse_alpha5_20230914.mat')
mX = res{1}.s_states;

if any(what_freq(cfg.measurement)=='m')
    dates_plot = datesm;
else
    dates_plot=datesq;
end
qt0 = 2;
mt0 = qt0*3;
figure(4)
tiledlayout(2,1)
nexttile;
plot(datesm(mt0:end), mX(mt0:end,1)*100,'Color',colblind(4,:),'LineWidth',1.5)
hold on;
plot(datesq(qt0:end), qX(qt0:end,1)*100,'diamond','MarkerFaceColor',colblind(4,:),'MarkerEdgeColor',colblind(4,:));
grid on
hold off
recessionplot;
legend({'Monthly','Quarterly'})
ylabel('Capital ($\%$ deviation from the steady state)','Interpreter', 'latex')

nexttile;
plot(datesm(mt0:end), mX(mt0:end,2)*100,'Color',colblind(1,:),'LineWidth',1.5)
hold on;
plot(datesq(qt0:end), qX(qt0:end,2)*100,'diamond','MarkerFaceColor',colblind(1,:),'MarkerEdgeColor',colblind(1,:))
grid on
hold off
recessionplot;
% legend()
ylabel('TFP','Interpreter', 'latex')
xlabel('Month', 'Interpreter','latex')
print('-depsc','smoothed_states.eps')


figure(5)
plot(datesm(mt0:end), mX(mt0:end,1)*100,'Color',colblind(4,:),'LineWidth',1.5)
hold on;
plot(datesq(qt0:end), qX(qt0:end,1)*100,'diamond','MarkerFaceColor',colblind(4,:),'MarkerEdgeColor',colblind(4,:));
grid on
plot(datesm(mt0:end), mX(mt0:end,2)*100,'Color',colblind(1,:),'LineWidth',1.5)
hold on;
plot(datesq(qt0:end), qX(qt0:end,2)*100,'diamond','MarkerFaceColor',colblind(1,:),'MarkerEdgeColor',colblind(1,:))
hold off
recessionplot;
ylabel('TFP','Interpreter', 'latex')
xlabel('Month', 'Interpreter','latex')
print('-depsc','smoothed_states_all.eps')
