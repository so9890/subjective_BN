%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour

% NEXT: iNCLUDE LAMBDAA AS A POLICY VARIABLE
%% include path to package
clc, clear
if isfile('/home/sonja/Documents/projects/Overconsumption/codding_model/Own/tools')
    cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')
else 
    cd('C:\Users\user\Documents\subjective_BN\codding_model\Own')
end

folder = '../tools';
addpath(genpath(folder))

% create folder to save figures
%if CHECK IF DOES OR DOES NOT EXIST 
mkdir('figures/Rep_agent');
mkdir('simulation_results');

%% indicators for model versions

indic.fullDisposal = 0; % == 0 if gov. revenues are fully consumed (baseline), ==1 if gov revenues are fully disposed of
indic.het_growth   = 1; % == 0 if there is equal growth across sectors, ==1 if growth in the sustainable sector is slower 
indic.util         = 0; % == 0 if uses CRRA with gammaa=1 (bgp compatible, hours do not react to wage changes); 
                        % == 1 if CRRA gamma!=1 KPR preferences (bgp),
                        % should also see that income and substitution
                        % effect cancel due to bgp compatibility
indic.withtarget   = 0; % ==1 if uses swf with target

indic.var          = string('zero');% 'zetaa'; % which parameter to change in simulations

%% read in model file, Objective fcn gov and variables

model_rep_agent;

%% initialise variables for parameterisation

% grids depending on variation
zetaa_calib=1.4;
gri.zetaa= [zetaa_calib];
gri.tauul=linspace(0,1.5,50);

%solution_LF=containers.Map;
%solution_Ramsey=containers.Map;

% number of periods for simulation 
T=101;
time=1:T; % vector of periods (1 is the initial period)

% initialise matrices to save results
y_simLF=zeros(length(y),T);
x_simLF=zeros(length(x),T);
W_simLF=zeros(1,T); % social welfare


y_simRam=zeros(length(y),T);
x_simRam=zeros(length(x),T);
opt_pol_sim=zeros(length(pol(ismember(pol, ['tauul']))),T);

% initial conditions technology
Ad=8;
Ac=4;

x_init = eval(x);

% read in parameter values
% calibrate model
[params, pols_num, ~]=params_bgp_rep_agent(symms.params, f, pol, indic, T, nan, zetaa_calib);

%% simulation: Laissez-faire
% read in parameters, write model as function to be solved numerically,
% solve for T periods

if ~isfile(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
     T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
     params(list.params=='thetad'), indic.het_growth, pols_num(list.pol=='tauul'), indic.util))
 
    sol_mat=zeros( length(y)+length(x)+1, T, length(gri.zetaa)); % for each value of zetaa have a matrix of variables (rows) over time (columns)
    for tt=1:length(gri.zetaa)
        indic.zetaa=gri.zetaa(tt);

        simulation; 
    end
    save(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
         T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
         params(list.params=='thetad'), indic.het_growth, pols_num(list.pol=='tauul'), indic.util)...
         ,'sol_mat')
else 
    fprintf('LF exists')
    load(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
     T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
     params(list.params=='thetad'), indic.het_growth, pols_num(list.pol=='tauul'), indic.util))
end
 
%% calibration emissions and emission targets
% use this to calibrate relation of production and output
[targets_num, E_vec]=calibration_emissions(squeeze(sol_mat(list.y=='yd',1,gri.zetaa==1.4)), symms.targets, T);

%% grid method to maximise welfare function 

% takes ages!

%% Primal approach read in model
[symms, Obj_ramPA]=primal_problem(y, x, list, symms, E);
[Obj_ramPA_dynamic ]=primal_problem_dynamic(y, x, list, symms, E);

%% simulation Ramsey

Ac0=x_init(list.x=='Ac');
Ad0=x_init(list.x=='Ad');

if ~isfile(sprintf('simulation_results/ControlsRamsey_hetgrowth%d_util%d_withtarget%d.mat',indic.het_growth, indic.util, indic.withtarget))

    % initial guess optimal policy and Lagrange multi
    mu_budget       = 1;
    mu_target       = 1;
    mu_opt_final    = 1;
    mu_rc           = 1;
    mu_imp          = 1;
    mu_defH         = 1; 
    tauul           = 0.7;
    guessLF         = 0; % for simulation in comp equilibrium

    
    
    for t=time

        %-- read in model equations with numeric parameter values; 
        %           only variables: policy and lagrange multipliers of gov. problem
        %           emission target is dynamic => every period a new
        %           problem!

        [model, model_param, varsModel, paramsModel]=model_eq_Ramsey...
                            (Obj_ram, symms.optim, [symms.params, symms.targets, E], [params, targets_num, E_vec(t)], symms.optim,...
                            'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'));
%        [model, model_param, varsModel, paramsModel]=model_eq_Ramsey...
%                             (Obj_ramPA, symms.optimPA, [symms.params, symms.targets, E], [params, targets_num, E_vec(t)], symms.optimPA,...
%                             'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'));

        %-- solve for optimal policy

        guess= eval(varsModel);

        modFF = @(x)Ramsey_model(x);
        options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

        [opt_pol, fval] = fsolve(modFF, guess, options);

        tauul=opt_pol(varsModel=='tauul'); %   could be vector or scalar; initial value for next round
        opt_pol_sim(:,t)=tauul;

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
    
    %% Plots

%% Ramsey versus laissez faire

% number of figures in row in subplot
nn=5;
% list of variables to be plotted
%list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x, "welfare"];
list.plot=[list.y list.x, "welfare"];

% combine variables into one matrix and list 

plottsLF=sol_mat(:,:,params(list.params=='zetaa')==zetaa_calib);
% welfare
welf_sim=log(y_simRam(list.y=='c',:))-(y_simRam(list.y=='hl',:)+...
        params(list.params=='zetaa').*y_simRam(list.y=='hh',:)).^(1+params(list.params=='sigmaa'))./(1+params(list.params=='sigmaa'));

plottsRam=[y_simRam; x_simRam; welf_sim ];



figure(3)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plottsLF(list.joint==list.plot(i),:), time, plottsRam(list.joint==list.plot(i),:), 'LineWidth', 1.6)
legend(sprintf('LF: %s', list.plot(i)), sprintf('Ramsey: %s', list.plot(i)), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
ytickformat('%.2f')
end

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
plot(time, plottsLF(list.joint=='yd',:)./plottsLF(list.joint=='yc',:), time, plottsRam(list.joint=='yd',:)./plottsRam(list.joint=='yc',:), 'LineWidth', 1.6)
legend(sprintf('LF: $y_d/y_c$'), sprintf('Ramsey: $y_d/y_c$'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
ytickformat('%.2f')
%set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+2)
plot(time, opt_pol_sim, 'LineWidth', 1.6)
legend(sprintf('Optimal $\\tau_l$'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
ytickformat('%.2f')

sgtitle('Laissez Faire versus Ramsey')
path=sprintf('figures/Rep_agent/Ram_LF_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad0,Ac0,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget);
saveas(gcf,path)

%% -- comparison params

% choose parameter of interest here
if indic.var=='zetaa'
    gri.params=gri.zetaa;
elseif indic.var=='tauul'
    gri.params=gri.tauul;
end

% number of figures in row in subplot
nn=3;
% list of variables to be plotted
list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x];
% combine variables into one matrix and list 

plotts1=solution_LF(string(gri.params(1)));
plotts2=solution_LF(string(gri.params(2)));

list.joint=[list.y, list.x];

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plotts1(list.joint==list.plot(i),:), time, plotts2(list.joint==list.plot(i),:))
title(sprintf('%s', indic.var))
legend(sprintf('%s, value%d', list.plot(i), round(gri.params(1)*10)), sprintf('%s, value%d', list.plot(i), round(gri.params(2)*10)))
end

%%

fig=figure(2);
plot(time, plotts1(list.joint=='yd',:)./plotts1(list.joint=='yc',:), time, plotts2(list.joint=='yd',:)./plotts2(list.joint=='yc',:), 'LineWidth', 1.6)

if indic.var=='zetaa'
    lgd=legend(sprintf('$y_d/y_c, zeta$ =%.2f', gri.params(1)), sprintf('$y_d/y_c, zeta$ =%.2f', gri.params(2)), 'Interpreter', 'latex');
elseif indic.var=='tauul'
    lgd=legend(sprintf('$y_d/y_c, tau_l$ %.2f', gri.params(1)), sprintf('$y_d/y_c tau_l$ %.2f', gri.params(2)), 'Interpreter', 'latex');
end

set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')

path=sprintf('figures/Rep_agent/Yd_Yc_ratio_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_fullDisp%d_HetGrowth%d_tauul%.3f_util%d.png', T-1, ...
 params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad,Ac,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.fullDisposal, indic.het_growth,pols_num(list.pol=='tauul'), indic.util);
saveas(fig, path)


%% -- evolution of economy over time

% number of figures in row in subplot
nn=3;
% list of variables to be plotted
list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x];
% combine variables into one matrix and list 

plotts=[y_simLF;x_simLF];
list.joint=[list.y, list.x];

figure(1)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plotts(list.joint==list.plot(i),:))
legend(sprintf('%s', list.plot(i)))
end

suptitle('BGPs')
path=sprintf('figures/Rep_agent/bgps_periods%d_eppsilon%d_zeta%d_Ad0%d_Ac0%d_thetac%d_thetad%d_fullDisp%d_HetGrowth%d_tauul%d_util%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa')*10, Ad,Ac,...
    params(list.params=='thetac')*10, round(params(list.params=='thetad')*10) , indic.fullDisposal, indic.het_growth, round(pols_num(list.pol=='tauul')*100), indic.util);
saveas(gcf,path)

%% aggregate price level

p_sim=(y_simLF(list.y=='pd',:).^(1-params(list.params=='eppsilon'))+y_simLF(list.y=='pc',:).^(1-params(list.params=='eppsilon'))).^(1/(1-params(list.params=='eppsilon')));
% check market clearing
if indic.fullDisposal==0
    demand_output= params(list.params=='psii').*(y_simLF(list.y=='xd',:)+y_simLF(list.y=='xc',:))+ y_simLF(list.y=='c',:)+ y_simLF(list.y=='G',:);
else
    demand_output= params(list.params=='psii').*(y_simLF(list.y=='xd',:)+y_simLF(list.y=='xc',:))+ y_simLF(list.y=='c',:);
end
figure(2)
plot(time, p_sim)

figure(3)
plot(time, demand_output, time, y_simLF(list.y=='Y',:))
legend('demand', 'supply')
