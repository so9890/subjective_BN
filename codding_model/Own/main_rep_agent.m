%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour



%% include path to package
clc, clear
folder='/home/sonja/Documents/projects/Overconsumption/codding_model/Own/tools';
addpath(genpath(folder))

cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')

% create folder to save figures
%if CHECK IF DOES OR DOES NOT EXIST 
mkdir('figures/Rep_agent');
mkdir('simulation_results');

%% indicators for model versions

indic.fullDisposal = 0; % == 0 if gov. revenues are fully consumed (baseline), ==1 if gov revenues are fully disposed of
indic.het_growth   = 0; % == 0 if there is equal growth across sectors, ==1 if growth in the sustainable sector is slower 
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
gri.zetaa= [1.4];
%gri.tauul=[0.181, 0.5];

solution_LF=containers.Map;
%solution_Ramsey=containers.Map;

% number of periods for simulation 
T=101;
time=1:T; % vector of periods (1 is the initial period)

% initialise matrices to save results
y_simLF=zeros(length(y),T);
x_simLF=zeros(length(x),T);

y_simRam=zeros(length(y),T);
x_simRam=zeros(length(x),T);
opt_pol_sim=zeros(length(pol(ismember(pol, ['tauul']))),T);

% initial conditions technology
Ad=8;
Ac=4;

x_init = eval(x);
%% simulation: Laissez-faire
% read in parameters, write model as function to be solved numerically,
% solve for T periods

% simulation with values of zetaa in gri.zetaa; tauul
for zz=1:length(gri.zetaa)
    indic.zetaa=gri.zetaa(zz);
    simulation;
end

% for tt=1:length(gri.tauul)
%     indic.tauul_ex=gri.tauul(tt);
%     simulation;
% end
%% simulation Ramsey

% calibrate model
[params, pols_num, ~, E_num]=params_bgp_rep_agent(symsparams, f, pol, indic, T, nan);

% initial guess optimal policy and Lagrange multi
mu_target= 0.2;
tauul = 0.1;
guessLF=0; % for simulation in comp equilibrium
 
for t=time

    %-- read in model equations with numeric parameter values; 
    %           only variables: policy and lagrange multipliers of gov. problem
    %           emission target is dynamic
    
    [model, model_param, varsModel, paramsModel]=model_eq_Ramsey...
                        (Obj_ram, symms.optim, [symsparams,symstargets], [params, E_num(t)], symms.optim,...
                        'Ramsey_model', x_init, x, pols_num(pol~='tauul'), pol(pol~='tauul'));
    
    %-- solve for optimal policy
    guess= eval(varsModel);
    modFF = @(x)Ramsey_model(x);
    options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-14);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

    [opt_pol, fval] = fsolve(modFF, guess, options);
    
    tauul=opt_pol(varsModel=='tauul'); %   could be vector or scalar; initial value for next round
    opt_pol_sim(:,t)=tauul;
    
    % find optimal allocation: pass optimal policy into comp. equilibrium
    % model
    [params, pols_num, model_pars, E]=params_bgp_rep_agent(symsparams, f, pol, indic, T, opt_pol_sim(:,t));
    [ybgp, xpbgp, solution]= simul_bgp(list, x, x_init, params, pols_num, model_pars, t, guessLF);

    % save results
    y_simRam(:,t)=ybgp';
    x_simRam(:,t)=x_init; % to save technology in correct period
    
    % update initial values and use as initial guess for solution
    x_init=xpbgp;
    guessLF=solution;
  
end

%% Plots

%% Ramsey versus laissez faire

% number of figures in row in subplot
nn=3;
% list of variables to be plotted
list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x];
% combine variables into one matrix and list 

plottsLF=[y_simLF; x_simLF];
plottsRam=[y_simRam; x_simRam];

list.joint=[list.y, list.x];

figure(3)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plottsLF(list.joint==list.plot(i),:), time, plottsRam(list.joint==list.plot(i),:), 'LineWidth', 1.6)
legend(sprintf('LF: %s', list.plot(i)), sprintf('Ramsey: %s', list.plot(i)), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
end

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
plot(time, plottsLF(list.joint=='yd',:)./plottsLF(list.joint=='yc',:), time, plottsRam(list.joint=='yd',:)./plottsRam(list.joint=='yc',:), 'LineWidth', 1.6)
legend(sprintf('LF: y_d/y_c'), sprintf('Ramsey: y_d/y_c'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
%set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+2)
plot(time, opt_pol_sim, 'LineWidth', 1.6)
legend(sprintf('Optimal $\tau_l$'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')


suptitle('Laissez Faire versus Ramsey')
path=sprintf('figures/Rep_agent/Ram_LF_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad,Ac,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util);
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
