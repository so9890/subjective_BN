%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour



%% include path to package
clc, clear
folder='/home/sonja/Documents/projects/Overconsumption/codding_model/Own/tools'
addpath(genpath(folder))

cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')

% create folder to save figures
%if CHECK IF DOES OR DOES NOT EXIST 
mkdir('figures/Rep_agent');

%% read in model file and variables
model_rep_agent;

%% read in parameters and write model as function to be solved numerically
[params, pols_num, model_pars]=params_bgp_rep_agent(symsparams, f, pol, list);

%% simulate model for 30 years (until 2050)

% number of periods for simulation 
T=151;
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
    [ybgp, xpbgp]= simul_bgp(hours, list, x, x_init, params, pols_num, model_pars, symsparams, pol);

    % save results
    y_sim(:,t)=ybgp';
    x_sim(:,t)=x_init; % to save technology in correct period
    
    % update initial values
    x_init=xpbgp;
    hours=[ybgp(list.y=='H'), ybgp(list.y=='hh'), ybgp(list.y=='hl')];
end

%save('x_sim')
%save('y_sim')
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
subplot(length(list.plot)/nn,nn,i)
plot(time, plotts(list.joint==list.plot(i),:))
legend(sprintf('%s', list.plot(i)))
end

suptitle('BGPs')
path=sprintf('figures/Rep_agent/bgps_periods%d_eppsilon%d_zeta%d.png', T-1, ...
    params(list.params=='eppsilon'),round(params(list.params=='zetaa')*10) );
saveas(gcf,path)
