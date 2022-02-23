if ~isfile(sprintf('simulation_results/StaticOptPol_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ))

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
        muu_target      = 1;
        kt_lab          = 0;
       % guessLF         = 0; % for simulation in comp equilibrium
end

    % initial guess optimal policy and Lagrange multi
   

    for t=1:P

        %-- read in model equations with numeric parameter values; 
        %           only variables: policy and lagrange multipliers of gov. problem
        %           emission target is dynamic => every period a new
        %           problem!

        if indic.approach==2

            [indexx, model, model_param, varsModel, paramsModel]=model_eq_Ramsey...
                            (Obj_ram, symms.optim, [symms.params, symms.targets, E], [params, targets_num, E_vec(t)], symms.optim,...
                            'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'), list, taul_sym);
        elseif indic.approach ==1

            [indexx, model, model_param, varsModel, paramsModel]=model_eq_Ramsey...
                             (Obj_ramPA, symms.optimPA, [symms.params, symms.targets, E], [params, targets_num, E_vec(t)], symms.optimPA,...
                             'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'), list, taul_sym);
        end

        %-- solve for optimal policy

        guess= eval(varsModel);
        % transform guess to unbounded variables
        guess_trans=trans_guess(indexx, guess); 

        modFF = @(x)Ramsey_model(x);
        options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

     
        [opt_pol, ~, exitf] = fsolve(modFF, guess_trans, options);
        fprintf('ramsey solved with %d, eppsilon %.2f, withtarget %d, period%d', exitf,params(list.params=='eppsilon'), indic.withtarget, t );

        % update initial guess
        opt_trans=trans_allo_out(indexx, opt_pol);
        multis(:,t) = opt_trans(varsModel=='muu_target'|varsModel=='kt_lab'); 
        tauul=opt_trans(varsModel=='tauul'); %   could be vector or scalar; initial value for next round
        opt_pol_sim(:,t)=opt_trans(varsModel=='tauul');
        vc           = pols_num(list.pol=='vc');
        vdd          = pols_num(list.pol=='vdd');
        vars_tosolve = [y,xp];
        
        ybgp=zeros(size(list.y));
        xpbgp=zeros(size(list.xp));

        % find optimal allocation: pass optimal policy into comp. equilibrium
        % model
        % FASTER: USE ANALYTICALLY SOLVED MODEL
        %[params, pols_num, model_pars]=params_bgp_rep_agent(symms.params, f, pol, indic, T, opt_pol_sim(:,t), zetaa_calib);
        %[ybgp, xpbgp, solution]= simul_bgp(list, x, x_init, params, pols_num, model_pars, t, guessLF);
        
        % analytic solution
        levels=solution_eqbm(x_init, params, eval(pol), list, vars_tosolve );
                        
            % save results
            for i =list.y
                ybgp(list.y==i) =levels(vars_tosolve==i);
            end
        
            for i = list.xp
                xpbgp(list.xp==i)= levels(vars_tosolve==i);
            end

        % save results
        y_simRam(:,t)=transpose(ybgp);
        x_simRam(:,t)=transpose(x_init); % to save technology in correct period

        % update initial values and use as initial guess for solution
        x_init=xpbgp;
        %guessLF=solution;

    end

save(sprintf('simulation_results/StaticControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
save(sprintf('simulation_results/StaticStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
save(sprintf('simulation_results/StaticOptPol_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'opt_pol_sim');

elseif isfile(sprintf('simulation_results/StaticOptPol_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ))
  
    fprintf('ramsey exists')
    
    load(sprintf('simulation_results/StaticControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
    load(sprintf('simulation_results/StaticStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
    load(sprintf('simulation_results/StaticOptPol_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'opt_pol_sim');
end
  