if indic.approach==1
        mu_budget       = 1;
        mu_target       = 1;
        mu_opt_final    = 1;
        mu_rc           = 1;
        mu_imp          = 1;
        mu_defH         = 1; 
        H               = 1;
        c               = 1;
        lhc             = 1;
        lhd             = 1;
        llc             = 1; 
        lld             = 1;
        yd              = 1;
        yc              = 1;
 
 elseif indic.approach==2
        tauul           = 0.7; 
        guessLF = 0; % for simulation in comp equilibrium
end

if ~isfile(sprintf('simulation_results/ControlsRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget))

    % initial guess optimal policy and Lagrange multi
   

    for t=time

        %-- read in model equations with numeric parameter values; 
        %           only variables: policy and lagrange multipliers of gov. problem
        %           emission target is dynamic => every period a new
        %           problem!

        if indic.approach==2

                 [model_dyn, model_param_dyn, varsModel_dyn, paramsModel_dyn]=model_eq_Ramsey...
                            (Obj_ram_dynamic, symms.optim, [symms.params, symms.targets, E], [params, targets_num, E_vec(t)], symms.optim,...
                            'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'));
   
        elseif indic.approach ==1

            
        end

        %-- solve for optimal policy

        guess= eval(varsModel);

        modFF = @(x)Ramsey_model(x);
        options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

        [opt_pol, fval] = fsolve(modFF, guess, options);
        fprintf('ramsey solved with %d ', fval )

        % update initial guess
        tauul=opt_pol(varsModel=='tauul'); %   could be vector or scalar; initial value for next round
        opt_pol_sim(:,t)=opt_pol(varsModel=='tauul');

        % find optimal allocation: pass optimal policy into comp. equilibrium
        % model
        % FASTER: USE ANALYTICALLY SOLVED MODEL
        [params, pols_num, model_pars]=params_bgp_rep_agent(symms.params, f, pol, indic, T, opt_pol_sim(:,t), zetaa_calib);
        [ybgp, xpbgp, solution]= simul_bgp(list, x, x_init, params, pols_num, model_pars, t, guessLF);

        % save results
        y_simRam(:,t)=ybgp';
        x_simRam(:,t)=x_init; % to save technology in correct period

        % update initial values and use as initial guess for solution
        x_init=xpbgp;
        guessLF=solution;

    end

    save(sprintf('simulation_results/PAControlsRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget),'y_simRam');
    save(sprintf('simulation_results/PAStatesRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget),'x_simRam');

else
    load(sprintf('simulation_results/ControlsRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget));
    load(sprintf('simulation_results/StatesRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget));
end
  