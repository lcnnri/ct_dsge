% display
labs   = {'psi'; 'alpha'; 'rho'; 'delta'; 'sigma_z';'rho_z';'eta'; 'sigma_k'; 'sigma_c'; 'sigma_n'; 'sigma_y'; 'sigma_r'};
labs = labs([true(8,1); res{1}.cfg.measurement]);
for i=1:numel(res)
    
    idx = logical(res{i}.cfg.par_ident);
    labels = labs(idx);
    model = repmat(res{i}.model,sum(idx),1);
    par_est = res{i}.theta_hat(idx);
    if res{i}.cfg.bootstrap_inference
        if res{i}.cfg.compute_hessians
            ci_est = res{i}.BS_stude_ci;
        else
            ci_est = res{i}.BS_basic_ci;
        end
        se_bs  = res{i}.BS_squared_se;
    else
        ci_est = res{i}.stude_ci;
        se_bs  = res{i}.se;
    end
    length = ci_est(:,2) - ci_est(:,1);
    res{i}.summary = table(model,labels, par_est, se_bs, ci_est, length);
    disp(res{i}.summary)
end

% SAVE RESULTS for posterity % ===========================================%
fname = ['output/workspaces/' serial_id];

fname = [fname ' - ' Panel];

save(fname, 'res', 'DATA')

%% PRINT RESULTS TO SCREEN 
if res{1}.cfg.measurement(4);
    num_models=4;
else
    num_models = 3;
end
idx = logical(res{1}.cfg.par_ident);
par_est = []; ses = [];
for i=1:numel(res)
    par_est = [par_est res{i}.theta_hat(idx)];
    ses     = [ses res{i}.BS_squared_se];
end
string = [];
order_in_table = 1:6;
order_in_table = order_in_table(1:size(par_est,1));
order_in_table(2) = 1;order_in_table(1)=2;
if numel(order_in_table)==6
    order_in_table = order_in_table([1:3  6]); 
end
par_est =  par_est(order_in_table,:);
ses = ses(order_in_table,:);
pars_to_display = size(order_in_table,2);
for j = 1:sum(pars_to_display)
    string = [string; 
        {sprintf('&  %2.4f ',  [par_est(j,:)]);
         sprintf('& (%2.4f) ', [ses(j,:)])}];
end
string(:) % print to screen and copy-paste in TeX


fid = fopen([ fname '.txt'],'wt');
for row = 1:size(string,1)
    fprintf(fid, '%s%s\n', string{row,:}, '\\');
end
fclose(fid);