for tt=1:length(gri.tauul)
    indic.tauul_ex=gri.tauul(tt);
    
    % this loop is to test the effect of tax progressivity on the output
    % share

%% calibrate model
[params, pols_num, model_pars, E]=params_bgp_rep_agent(symsparams, f, pol, indic, T, nan);

%% simulate model for T years

% initial guess => empty because is updated in code in initial iteration 
guess0 = 0;
% iterate over solving model taking previous period technology and hours
% worked as initial conditions

for t=time
    [ybgp, xpbgp, solution]= simul_bgp( list, x, x_init, params, pols_num, model_pars, t, guess0);

    % save results
    y_simLF(:,t)=ybgp';
    x_simLF(:,t)=x_init; % to save technology in correct period
    
    % update initial values and use as initial guess for solution
    x_init=xpbgp;
    guess0=solution;
end

solution_LF(string(indic.tauul_ex))=[y_simLF;x_simLF];

% save simulations
save(sprintf('simulation_results/states_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d_tauul%d_util%d', ...
    indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10, ...
    params(list.params=='thetad')*10, indic.het_growth, round(indic.tauul_ex*100), indic.util)...
    ,'x_sim')
save(sprintf('simulation_results/controls_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d_tauul%d_util%d', ...
    indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10,...
    params(list.params=='thetad')*10, indic.het_growth, round(indic.tauul_ex*100), indic.util)...
    , 'y_sim')
end