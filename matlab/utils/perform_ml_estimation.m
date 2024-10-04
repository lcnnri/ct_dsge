restrtheta0 = nan(size(theta0));
for i=1:numel(theta0)
    trtheta0(i,1) = trans(theta0(i), trspec(i));
    cfg.pars(i,1) = trans(theta0(i), trspec(i));
end
trspecest = trspec(cfg.par_ident==1);
pd = makedist('Normal');
z_ = abs(icdf(pd,alpha/2)); % z score (critical value of standard normal)



% estimation
fprintf('Empirical Estimation------\n')
for i_model=1:numel(models)
    cfg.model = ... assign representation
        models{i_model}; 
    
    % sanity checks 
    fprintf('Model %s------\n',cfg.model)
    if strcmpi(cfg.model,'MX-SSR') || strcmpi(cfg.model,'F-SSR')
        cfg.flows=true;
        cfg.mixed_sampling=false;
        if strcmpi(cfg.model,'MX-SSR')
            cfg.mixed_sampling=true;
        end
    else
        cfg.flows=false;
    end
    if strcmpi(cfg.model,'EM-SSR')
        cfg.flows=false; cfg.em=true;
        cfg.mixed_sampling=false;
    end
    if strcmpi(cfg.model,'S-SSR')
        cfg.flows=false; cfg.em=false;
        cfg.mixed_sampling=false;
    end
    
    % trasnform the data according to cfg
    [y, cfg] = transform_data(DATA',cfg);
    objective = @(theta) estimate_theta2(y,cfg,theta,minimal_system_repr);


    if knitro_toolbox % ML Estimation
        [trtheta_hat,fval,exitflag,output,lambda,grad,hess] = ...
            knitro_nlp(@(theta)estimate_theta2(y,cfg,theta,minimal_system_repr), ...
            trtheta0(logical(cfg.par_ident)), [],[],[],[],...
            [],[],[],[],options);
    else
        [trtheta_hat,fval,exitflag,output,grad,hess] = ...
            fmincon(@(theta)objective(theta), ...
            trtheta0(logical(cfg.par_ident)), ...
            [],[],[],[],[],[],[], ...
            options);
    end
    
    % transform estimated parameters back to constrained values
    trtheta = trtheta0;
    trtheta(logical(cfg.par_ident)) = trtheta_hat; %store results
    for i=1:numel(trtheta_hat)
        xhat(i) = invtrans(trtheta_hat(i),trspecest(i));
    end
    theta_hat = theta0;
    theta_hat(cfg.par_ident==1) = xhat;
    clear xhat;
    % get the estimated structural shocks
    [~, x, P,ll, struct_shocks,v,K,xs,invF,sys] =  estimate_theta2(y,cfg,trtheta_hat,0);
    T = length(y);

    %% calculate hessian and grads for asymptotic se
     disp('Computing derivatives...')
     calculating_ses=true;
    parfor (t=2:T,num_workers)
%         for t=2:T % debug
        if knitro_toolbox
            [~,~,~,~,~,gLL(t,:)] = ...
                knitro_nlp(@(theta)estimate_theta2(y(:,1:t),cfg,theta,minimal_system_repr,calculating_ses), ...
                theta_hat(cfg.par_ident==1), [],[],[],[],...
                theta_hat(cfg.par_ident==1),theta_hat(cfg.par_ident==1),[],[],options_no_verbose);
        else
            [~,~,~,~,~,gLL(t,:)] = ...
                fmincon(@(theta)estimate_theta2(y(:,1:t),cfg,theta,minimal_system_repr,calculating_ses), ...
                theta_hat(cfg.par_ident==1), [],[],[],[],...
                theta_hat(cfg.par_ident==1),theta_hat(cfg.par_ident==1),[],options_no_verbose);
        end
    end

    % derivatives for studentized confidence intervals
    cfg.pars_invtrans = theta_hat';
    T = length(gLL); GLL = zeros(sum(cfg.par_ident));
    for t=1:length(gLL)
        GLL = GLL + gLL(t,:)' * gLL(t,:)/T;
    end
    hessian_mat = -hessian('getLL', theta_hat(cfg.par_ident==1), y, cfg,true,minimal_system_repr);
    GGLL = zeros(sum(cfg.par_ident));
    for t=1:size(hessian_mat,1)
        GGLL = GGLL + reshape(hessian_mat(t,:)',sum(cfg.par_ident),sum(cfg.par_ident));
    end
    OMEGA_T = GGLL*inv(GLL)*GGLL;
    se      = sqrt(diag(OMEGA_T)/T);

    % store results
    res{i_model}.DATA           = y;
    res{i_model}.dates_month   = datesm;
    res{i_model}.dates_quart   = datesq;
    res{i_model}.theta_hat      = theta_hat;
    res{i_model}.model          = cfg.model;
    res{i_model}.x              = x';
    res{i_model}.LL             = -fval;
    res{i_model}.exitflag       = exitflag;
    res{i_model}.output         = output;
    res{i_model}.grad           = grad;
    res{i_model}.par_est        = cfg.par_ident;
    res{i_model}.trspec         = trspec;
    res{i_model}.trtheta_hat    = trtheta;
    res{i_model}.i_logL         = -ll;
    res{i_model}.struct_shocks  = struct_shocks;
    res{i_model}.pred_errors    = v;
    res{i_model}.Kgain          = K;
    res{i_model}.P              = P;
    res{i_model}.s_states       = xs';
    res{i_model}.struct_shocks  = struct_shocks;
    res{i_model}.invF           = invF;
    res{i_model}.Ah             = sys.A;
    res{i_model}.B              = sys.B;
    res{i_model}.C              = sys.C;
    res{i_model}.D              = sys.D;
    res{i_model}.Sig_eps        = sys.Sig_eps;
    res{i_model}.Sig_eta        = sys.Sig_eta;
    res{i_model}.cfg            = cfg;
    res{i_model}.se             = se;
    res{i_model}.Omega_T        = OMEGA_T;
    res{i_model}.stude_ci       = [theta_hat(cfg.par_ident==1) - z_*se, ...
        theta_hat(cfg.par_ident==1) + z_*se];

    % cfg.minimal_representation=true;
[~, ~, ~,~, ~,~,Kh,~ ,invF,sys] = estimate_theta2(y,cfg,theta_hat(logical(cfg.par_ident)), ...
    true, ...
    true, ...
    false);


   % [A, b, C, D, B, Sig_eps, eigviol] = ... get the system matrices
   %      prepareFilter_loglin(theta_hat,[],cfg); % continuous-time system matrices
   % 
   %  [idFlag, OCflag] = ... identification, cf. Proposition 2
   %      isIdentified(theta_hat,cfg,y,'a'); 
   %  % [idFlag, OCflag] = ...
   %  %   isIdentified(theta_hat,cfg,y,'b');
   % 
   %  res{i_model}.IS_IDENTIFIED = idFlag;
   %  res{i_model}.IS_ORDERCOND  = OCflag;
end
%
if 1==2
    labs   = {'psi'; 
              'alpha'; 
              'rho'; 
              'delta'; 
              'sigma_z';
              'rho_z';
              'eta'; 
              'sigma_k'; 
              'sigma_c'; 
              'sigma_n'; 
              'sigma_y'; 
              'sigma_r'};
    idx = logical([ones(8,1); cfg.measurement]);
    MX_SSR = res{1}.theta_hat;
    F_SSR = res{2}.theta_hat;
    S_SSR = res{3}.theta_hat;
    EM_SSR = res{4}.theta_hat;

    table(labs(idx), MX_SSR, F_SSR, S_SSR,EM_SSR )

    if any(what_freq(cfg.measurement)=='m')
        dates_plot = datesm;
    else
        dates_plot=datesq;
    end

    figure(4)
    tiledlayout(2,1)
    nexttile;
    plot(dates_plot(2:end), res{1}.x(:,1)*100,'Color',colblind(4,:),'LineWidth',1.5)
    hold on;
    plot(dates_plot(2:end), res{2}.x(:,1)*100,'Color',colblind(1,:),'LineWidth',1.5);
    grid on
    hold off
    recessionplot;
    legend(models(1:2))
    ylabel('Capital ($\%$ deviation from the steady state)','Interpreter', 'latex')

    nexttile;
    plot(dates_plot(2:end), res{1}.x(:,2),'Color',colblind(4,:),'LineWidth',1.5)
    hold on;
    plot(dates_plot(2:end), res{2}.x(:,2),'Color',colblind(1,:),'LineWidth',1.5);
    grid on
    hold off
    recessionplot;
    legend(models(1:2))
    ylabel('TFP','Interpreter', 'latex')
end

if bootstrap_inference % if bootstrap inference to run

    if knitro_toolbox % suppress output for speed
        options = knitro_options('outlev',0);
    else
        options.Display = 'off';
    end
    for i_model = 1:numel(models) % iterate through models
        model = models{i_model};
        
        % extract values and configurations
        trspec  = res{i_model}.trspec;
        cfg     = res{i_model}.cfg;

        theta_ML   = res{i_model}.theta_hat;
        trtheta_ML = res{i_model}.trtheta_hat;
        trspecest  = trspec(logical(cfg.par_ident));

        % initialize storage
        trtheta_hat_bs                      = nan(BS, numel(trtheta_ML));
        trtheta_hat_bs(:,~cfg.par_ident)    = repmat(trtheta0(~cfg.par_ident)',BS,1);
        theta_hat_bs                        = nan(BS, numel(theta_ML));
        theta_hat_bs(:,~cfg.par_ident)      = repmat(theta_ML(~cfg.par_ident)',BS,1);
        fval_bs                             = nan(BS,1);
        exitflag_bs                         = nan(BS,1);
        grad_bs                             = nan(sum(cfg.par_ident),1,BS);
        cfg.pars_invtrans                   = res{i_model}.theta_hat;
        idx_par_est                         = (cfg.par_ident==1);
        bad_draw                            = zeros(BS,1);
        %         w = randn(T,n,BS);
        parfor(b=1:BS, num_workers)
        % for b=1:BS % debug
            disp([cfg.model ': bootstrap no.' num2str(b)]);
            
            % bootstrap the state-space
            [YY,e,XX] = ...
                bootstrap_Y2(y, ... original data
                cfg, ... configuration settings
                trtheta_ML(idx_par_est)', ... MLE 
                minimal_system_repr, ... use the minimial state representation flag
                bs_type ... type of bootstrap (nonparametric / wild)
                );   

            % estimate structural parameters from bootstrapped data
            try
                if knitro_toolbox
                    [b_trtheta_hat,b_fval,b_exitflag,output,~,b_grad,b_hess] = ...
                        knitro_nlp(@(theta)estimate_theta2(YY,cfg,theta,minimal_system_repr), ...
                        trtheta_ML(logical(cfg.par_ident))', [],[],[],[],...
                        [],[],[],[],options);
                else
                    [b_trtheta_hat,b_fval,b_exitflag,output,~,b_grad,b_hess] = ...
                        fmincon(@(theta)estimate_theta2(YY,cfg,theta,minimal_system_repr), ...
                        trtheta_ML(logical(cfg.par_ident))', [],[],[],[],[],[],[], ...
                        options);
                    
                end
            catch % if b is a bad bootstrap sample -> discard
                bad_draw(b)   = 1;
                b_trtheta_hat = cfg.pars(logical(cfg.par_ident));
                b_fval        = Inf;
                b_exitflag    = 999;
                output        = {};
                b_grad        = zeros(sum(cfg.par_ident),1);
                b_hess        = [];
            end
            
            %store results
            b_xhat =zeros(numel(b_trtheta_hat),1);
            for i=1:numel(b_trtheta_hat)
                b_xhat(i) = invtrans(b_trtheta_hat(i),trspecest(i));
            end
            b_theta_hat                   = theta0;
            b_theta_hat(idx_par_est)      = b_xhat;
            trtheta_hat_bs(b,idx_par_est) = b_trtheta_hat;
            fval_bs(b,1)                  = b_fval;
            exitflag_bs(b,1)              = b_exitflag;
            grad_bs(:,:,b)                = b_grad;
            theta_hat_bs(b,:)             = b_theta_hat;

            if compute_hessians % bootstrap asymptotic se:
                % derivatives for studentized confidence intervals
                % ---THIS TAKES A WHILE---
                T = length(YY); gLL = zeros(T,sum(cfg.par_ident));
                for t=2:T
                    if knitro_toolbox
                        [~,~,~,~,~,gLL(t,:)] = ...
                            knitro_nlp(@(theta)estimate_theta2(YY(:,1:t),cfg,theta,minimal_system_repr,calculating_ses), ...
                            b_theta_hat(idx_par_est), [],[],[],[],...
                            b_theta_hat(idx_par_est), b_theta_hat(idx_par_est),[],[],options_no_verbose);
                    else
                        [~,~,~,~,~,gLL(t,:)] = ...
                            fmincon(@(theta)estimate_theta2(YY(:,1:t),cfg,theta,minimal_system_repr,calculating_ses), ...
                            b_theta_hat(idx_par_est) , [],[],[],[],...
                            b_theta_hat(idx_par_est) ,b_theta_hat(idx_par_est) ,[],options_no_verbose);
                    end
                end

                GLL = zeros(sum(cfg.par_ident));
                for t=1:length(gLL)
                    GLL = GLL + gLL(t,:)' * gLL(t,:)/T;  % expectation of inner product of gradients
                end
                
                % hessian 
                hessian_mat = -hessian('getLL', b_xhat,YY, cfg,1,1);
                GGLL = zeros(sum(cfg.par_ident));
                for t=1:size(hessian_mat,1)
                    GGLL = GGLL + reshape(hessian_mat(t,:)',sum(cfg.par_ident),sum(cfg.par_ident));
                end
                OMEGA_T         = GGLL*inv(GLL)*GGLL;    % covariance matrix of estimates (QML)
                se_bs(b,:)      = sqrt(diag(OMEGA_T)/T); % se
                tstats_bs(b,:)  = (b_xhat' - theta_ML(cfg.par_ident==1)') ./ (se_bs(b,:)); % t-stat
            else
                se_bs(b,:)=0;
                tstats_bs(b,:)=0;
            end

        end
        % store bootstrap results
        res{i_model}.BS_trtheta_hat = ... transformed bootstrap parameter estimates
            trtheta_hat_bs;  
        res{i_model}.BS_theta_hat   = ... bootstrap structural parameter estimates
            theta_hat_bs;    
        res{i_model}.BS_theta_hat_biasCorrect = ... bias-correction cf. Tang and Chen (2009, JoE)
            2*res{i_model}.theta_hat - mean(theta_hat_bs,1)';
        res{i_model}.BS_fval        = ... negative log-likelihoods
            fval_bs;  
        res{i_model}.BS_exitflag    = ... exitflags from estimation (nonnegs are good)
            exitflag_bs;  
        res{i_model}.BS_grad        = ... gradient from the optim routine
            grad_bs; 
        res{i_model}.BS_num         = ... number of bootstrap samples generated
            BS; 
        res{i_model}.BS_num_workers = ... number of workers in the parallelization
            num_workers;  
        res{i_model}.BS_se          = ... analytical standard errors of bootstrap samples; if ~compute_hessians -> se_bs=0;
            se_bs;  
        res{i_model}.BS_perce_ci    = ... bootstrap percentile confidence intervals
            [quantile(theta_hat_bs(:,cfg.par_ident==1),alpha/2,1)', ...
            quantile(theta_hat_bs(:,cfg.par_ident==1),1-alpha/2,1)'];
        res{i_model}.BS_basic_ci    = ... bootstrap basic confidence intervals
            [2*res{i_model}.theta_hat(cfg.par_ident==1) - quantile(res{i_model}.BS_theta_hat(:,cfg.par_ident==1),1-alpha/2,1)', ...
             2*res{i_model}.theta_hat(cfg.par_ident==1) - quantile(res{i_model}.BS_theta_hat(:,cfg.par_ident==1),alpha/2,1)'];
        res{i_model}.BS_stude_ci    = ... bootstrap studentized confidence intervals
            [res{i_model}.theta_hat(cfg.par_ident==1)  + quantile(tstats_bs,alpha/2,1)'.*res{i_model}.se, ...
             res{i_model}.theta_hat(cfg.par_ident==1)  + quantile(tstats_bs,1-alpha/2,1)'.*res{i_model}.se];
        res{i_model}.BS_squared_se  = ... bootstrap squared standard errors
            std(res{i_model}.BS_theta_hat(:,cfg.par_ident==1),[],1)';
        res{i_model}.BS_bad_draw    = ... collection of bad bootstrap draws
            bad_draw; 
    end

end