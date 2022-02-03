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

indic.fullDisposal = 1; % == 0 if gov. revenues are fully consumed (baseline), ==1 if gov revenues are fully disposed of
indic.het_growth   = 0; % == 0 if there is equal growth across sectors, ==1 if growth in the sustainable sector is slower 

%% read in model file and variables
if indic.fullDisposal==0
    model_rep_agent;
else
    model_rep_agent_fullDisposal;
end

%% read in parameters and write model as function to be solved numerically
[params, pols_num, model_pars]=params_bgp_rep_agent(symsparams, f, pol, indic);

%% simulate model for 30 years (until 2050)

% number of periods for simulation 
T=201;
time=1:T; % vector of periods (1 is the initial period)

% initialise matrices to save results
y_sim=zeros(length(y),T);
x_sim=zeros(length(x),T);

% initial guesses for hours supplied

H=(1-pols_num(list.pol=='tauul'))^(1/(1+params(list.params=='sigmaa')));
hh=H/4; % one fourth in high skill hours
hl=H-params(list.params=='zetaa')*hh;

hours=[H, hh, hl];
list.hours=["H", "hh", "hl"];

% initial conditions technology
Ad=8;
Ac=4;

x_init=eval(x);

% iterate over solving model taking previous period technology and hours
% worked as initial conditions

for t=time
    [ybgp, xpbgp]= simul_bgp(hours, list, x, x_init, params, pols_num, model_pars, indic);

    % save results
    y_sim(:,t)=ybgp';
    x_sim(:,t)=x_init; % to save technology in correct period
    
    % update initial values
    x_init=xpbgp;
    hours=[ybgp(list.y=='H'), ybgp(list.y=='hh'), ybgp(list.y=='hl')];
end

% save simulations
save(sprintf('simulation_results/states_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d', ...
    indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10, ...
    params(list.params=='thetad')*10, indic.het_growth)...
    ,'x_sim')
save(sprintf('simulation_results/controls_fullDisposal%d_T%d_initialAd%dAc%d_eppsilon%d_zetaa%d_thetac%d_thetad%d_HetGrowt%d', ...
    indic.fullDisposal, T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa')*10, params(list.params=='thetac')*10,...
    params(list.params=='thetad')*10, indic.het_growth)...
    , 'y_sim')

%% Plots
%-- evolution of economy over time

% number of figures in row in subplot
nn=3;
% list of variables to be plotted
list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x];
% combine variables into one matrix and list 
plotts=[y_sim;x_sim];
list.joint=[list.y, list.x];

figure(1)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plotts(list.joint==list.plot(i),:))
legend(sprintf('%s', list.plot(i)))
end

suptitle('BGPs')
path=sprintf('figures/Rep_agent/bgps_periods%d_eppsilon%d_zeta%d_Ad0%d_Ac0%d_thetac%d_thetad%d_fullDisp%d_HetGrowth%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa')*10, Ad,Ac,...
    params(list.params=='thetac')*10, params(list.params=='thetad')*10 , indic.fullDisposal, indic.het_growth);
saveas(gcf,path)

%% aggregate price level

p_sim=(y_sim(list.y=='pd',:).^(1-params(list.params=='eppsilon'))+y_sim(list.y=='pc',:).^(1-params(list.params=='eppsilon'))).^(1/(1-params(list.params=='eppsilon')));
% check market clearing
if indic.fullDisposal==0
    demand_output= params(list.params=='psii').*(y_sim(list.y=='xd',:)+y_sim(list.y=='xc',:))+ y_sim(list.y=='c',:)+ y_sim(list.y=='G',:);
else
    demand_output= params(list.params=='psii').*(y_sim(list.y=='xd',:)+y_sim(list.y=='xc',:))+ y_sim(list.y=='c',:);
end
figure(2)
plot(time, p_sim)

figure(3)
plot(time, demand_output, time, y_sim(list.y=='Y',:))
legend('demand', 'supply')
