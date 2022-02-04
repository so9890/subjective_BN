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
indic.het_growth   = 1; % == 0 if there is equal growth across sectors, ==1 if growth in the sustainable sector is slower 
indic.util         = 0; % == 0 if uses CRRA with gammaa=1 (bgp compatible, hours do not react to wage changes); 
                        % == 1 if CRRA gamma!=1 KPR preferences (bgp),
                        % should also see that income and substitution
                        % effect cancel due to bgp compatibility


%% read in model file and variables
%if indic.fullDisposal==0
    model_rep_agent;
%else
%    model_rep_agent_fullDisposal;
%end

%% read in parameters and write model as function to be solved numerically
gri.tauul=[0.181, 0.5];
%keyys= cellfun(@num2str,num2cell(gri.tauul(:)),'uniformoutput',false);
mapp=containers.Map;

for tt=1:length(gri.tauul)
    indic.tauul_ex=gri.tauul(tt);
    
    % this loop is to test the effect of tax progressivity on the output
    % share
    
[params, pols_num, model_pars]=params_bgp_rep_agent(symsparams, f, pol, indic);

%% simulate model for 30 years (until 2050)

% number of periods for simulation 
T=11;
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

x_init = eval(x);

% initial guess => empty because is updated in code in initial iteration 
guess0 = 0;
% iterate over solving model taking previous period technology and hours
% worked as initial conditions

for t=time
    [ybgp, xpbgp, solution]= simul_bgp(hours, list, x, x_init, params, pols_num, model_pars, indic, t, guess0);

    % save results
    y_sim(:,t)=ybgp';
    x_sim(:,t)=x_init; % to save technology in correct period
    
    % update initial values and use as initial guess for solution
    x_init=xpbgp;
    hours=[ybgp(list.y=='H'), ybgp(list.y=='hh'), ybgp(list.y=='hl')];
    guess0=solution;
end

mapp(string(indic.tauul_ex))=[y_sim;x_sim];

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

%% Plots

%% -- comparison tauul
% number of figures in row in subplot
nn=3;
% list of variables to be plotted
list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd" "lld"])), list.x];
% combine variables into one matrix and list 

plotts1=mapp(string(gri.tauul(1)));
plotts2=mapp(string(gri.tauul(2)));

list.joint=[list.y, list.x];

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plotts1(list.joint==list.plot(i),:), time, plotts2(list.joint==list.plot(i),:))
legend(sprintf('%s_tauul%d', list.plot(i), gri.tauul(1)), sprintf('%s_tauul%d', list.plot(i), gri.tauul(2)))
end
%%
for i=1:length(list.plot)
%subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plotts1(list.joint=='yd',:)./plotts1(list.joint=='yc',:), time, plotts2(list.joint=='yd',:)./plotts2(list.joint=='yc',:))
legend(sprintf('tauul%d', gri.tauul(1)), sprintf('tauul%d', gri.tauul(2)))
end


%% -- evolution of economy over time

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
path=sprintf('figures/Rep_agent/bgps_periods%d_eppsilon%d_zeta%d_Ad0%d_Ac0%d_thetac%d_thetad%d_fullDisp%d_HetGrowth%d_tauul%d_util%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa')*10, Ad,Ac,...
    params(list.params=='thetac')*10, round(params(list.params=='thetad')*10) , indic.fullDisposal, indic.het_growth, round(pols_num(list.pol=='tauul')*100), indic.util);
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
