if indic.approach==1
%         mu_budget       = 1;
%         mu_target       = 1;
%         mu_opt_final    = 1;
%         mu_rc           = 1;
%         mu_imp          = 1;
%         mu_defH         = 1; 
%         H               = 1;
%         c               = 1;
%         lhc             = 1;
%         lhd             = 1;
%         llc             = 1; 
%         lld             = 1;
%         yd              = 1;
%         yc              = 1;
 
 elseif indic.approach==2
        guess           = repmat(0.7, size(symms.optim_dynamic)); %subs(symms.optim_dynamic,zer )
end

    % initial guess optimal policy and Lagrange multi
  

        %-- read in model equations with numeric parameter values; 
        %           only variables: policy and lagrange multipliers of gov. problem
        %           emission target is dynamic => every period a new
        %           problem!

if indic.approach==2

     [model_dyn, model_param_dyn, varsModel_dyn, paramsModel_dyn]=model_eq_Ramsey...
                (Obj_ram_dynamic, symms.optim_dynamic, [symms.params, symms.targets, transpose(vecs(:,list.joint=='E'))], [params, targets_num, E_vec], symms.optim_dynamic,...
                'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'));

elseif indic.approach ==1

fprintf('primal approach dynamic not yet coded')
end

%-- solve for optimal policy

modFF = @(x)Ramsey_model(x);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[opt_pol_sim, ~, exitflag] = fsolve(modFF, guess, options);
fprintf('ramsey solved with exitflag %d ', exitflag );


% optimal allocation from analytic solution: static 
vc           = pols_num(list.pol=='vc');
vdd          = pols_num(list.pol=='vdd');
vars_tosolve = [y,xp];

for t=1:P
% find optimal allocation: pass optimal policy into comp. equilibrium
% x= states, params = numeric parameters, pols = numeric vector of
% solution
    tauul=opt_pol_sim(t);
   
    levels=solution_eqbm(x_init, params, eval(pol), list, vars_tosolve );
                
    % save results
    ybgp=zeros(size(list.y));
    xpbgp=zeros(size(list.xp));
    for i =list.y
        ybgp(list.y==i) =levels(vars_tosolve==i);
    end

    for i = list.xp
        xpbgp(list.xp==i)= levels(vars_tosolve==i);
    end

    y_simRam(:,t)=transpose(ybgp);
    x_simRam(:,t)=x_init; % to save technology in correct period

    % update initial values and use as initial guess for solution
    x_init=xpbgp;
end


save(sprintf('simulation_results/DynamicControlsRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget),'y_simRam');
save(sprintf('simulation_results/DynamicStatesRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget),'x_simRam');

  