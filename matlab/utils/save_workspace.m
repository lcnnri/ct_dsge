pars = {'psi','alpha','delta','rho','sigz','rhoz','eta','sigk','sig_eps1','sig_eps2','sig_eps3'};
xtrastr = [];
% for i=1:length(cfg.par_ident)
%     if cfg.par_ident(i)==0
%         if i~=8
%             xtrastr = [xtrastr '_fix_' pars{i}];
%         end
%     end
% end
if ~cfg.flows 
    xtrastr = [xtrastr '_STOCKTREATMENT'];
end
if  ~isfield(cfg,'em') 
    cfg.em = false;
end
if cfg.em 
    xtrastr = [xtrastr '_EULERAPPROX'];
end

if exist('wrapper') && exist('stock_data')
    if stock_data
       fname = [ num2str(reps) 'MC_stock_' cfg.strmeasurement xtrastr '.mat'];
    else
        fname = [ num2str(reps) 'MC_flows_' cfg.strmeasurement xtrastr '.mat'];
    end
elseif flowdata
    fname = [ num2str(reps) 'MC_flows_' cfg.strmeasurement xtrastr '.mat'];
else
    fname = [ num2str(reps) 'MC_stock_' cfg.strmeasurement xtrastr '.mat'];
end

if exist('wrapper') && exist('wrapper_data')
    if strcmpi(wrapper_data(end-10:end),'monthly_240')
        fname = [fname(1:end-4) '240_measurements_monthly.mat'];
    elseif strcmpi(wrapper_data(end-9:end),'yearly_240')
        fname = [fname(1:end-4) '240_measurements_yearly.mat'];
    end
    if strcmpi(wrapper_data(end-6:end),'monthly')
        if strcmpi(wrapper_data(end-15:end-8), 'prec_120')
            fname = [fname(1:end-4) '_prec_120_monthly.mat']
        else
            fname = [fname(1:end-4) '_monthly.mat']
        end
    elseif strcmpi(wrapper_data(end-5:end),'yearly')
        if strcmpi(wrapper_data(end-14:end-7), 'prec_120')
            fname = [fname(1:end-4) '_prec_120_yearly.mat']
        else
            fname = [fname(1:end-4) '_yearly.mat']
        end
    end
end
if isfield(cfg,'serial_id')
    fname = [cfg.serial_id '_' fname];
end
fname = ['../output/simulation_results/' fname]; 
save(fname, 'output', 'data', 'reps', 'S', 'cfg', 'xitflag')