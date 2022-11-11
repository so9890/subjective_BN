
    % this loop is to test the effect of tax progressivity on the output
    % share

%% calibrate model
[params, pols_num, model_pars]=params_bgp_rep_agent(symms.params, f, pol, indic, T, nan, zetaa_calib);

%% simulate model for T years

% initial guess => empty because is updated in code in initial iteration 
guess0 = 0;
% iterate over solving model taking previous period technology and hours
% worked as initial conditions

for t=time
    [ybgp, xpbgp, solution]= simul_bgp( list, x, x_init, params, pols_num, model_pars, t, guess0);

    % welfare
    if indic.util==0
    W_simLF(:,t)=log(ybgp(list.y=='c'))-(ybgp(list.y=='hl')+...
        params(list.params=='zetaa')*ybgp(list.y=='hh'))^(1+params(list.params=='sigmaa'))/(1+params(list.params=='sigmaa'));
    end
    % save results
    y_simLF(:,t)=ybgp';
    x_simLF(:,t)=x_init; % to save technology in correct period
    
    % update initial values and use as initial guess for solution
    x_init=xpbgp;
    guess0=solution;
end

%solution_LF(tt)=[y_simLF;x_simLF];
helper=[y_simLF;x_simLF;W_simLF];

if ~isnan(sol_mat)
    sol_mat(:,:,tt)= reshape(helper,...
        [size(helper,1),size(helper,2),1]);
end
% save simulations
% save(sprintf('simulation_results/states_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d_tauul%d_util%d', ...
%     indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10, ...
%     params(list.params=='thetad')*10, indic.het_growth, round(pols_num(list.pol=='tauul')*100), indic.util)...
%     ,'x_simLF')
% save(sprintf('simulation_results/controls_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d_tauul%d_util%d', ...
%     indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10,...
%     params(list.params=='thetad')*10, indic.het_growth, round(pols_num(list.pol=='tauul')*100), indic.util)...
%     , 'y_simLF')
