%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour

%% include path to package
clc, clear
if isfile('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')
    cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')
else 
    cd('C:\Users\user\Documents\subjective_BN\codding_model\Own')
end
%%
folder = '../tools';
addpath(genpath(folder))

% create folder to save figures
mkdir('figures/Rep_agent');
mkdir('simulation_results');

%% indicators for model versions

%for ss=0:1
%for ttt=0:1
indic.subst        = 1; %ss % == 0 if complements, ==1 if substitutes
indic.epps         = [0.4, 4]; % vector of values for eppsilon, depending on indic.subst
indic.fullDisposal = 0; % == 0 if gov. revenues are fully consumed (baseline), ==1 if gov revenues are fully disposed of
indic.het_growth   = 1; % == 0 if there is equal growth across sectors, ==1 if growth in the sustainable sector is slower 
indic.util         = 2; % == 0 if uses CRRA with gammaa=1 (bgp compatible, hours do not react to wage changes); 
                        % == 1 if CRRA gamma!=1 KPR preferences (bgp),
                        % should also see that income and substitution
                        % effect cancel due to bgp compatibility
                        % ==2 if uses MaCurdy preferences (following Boppart and Krusell) 

indic.withtarget   = 0;  %ttt; % ==1 if uses swf with target; ==0 if no target
indic.approach     = 2; % ==1 if uses primal approach, ==2 if uses dual approach (maxmise over optimal policy (tauul, lambdaa) directly
indic.var          = string('zero');% 'zetaa'; % which parameter to change in simulations

%% read in model file, Objective fcn gov and variables

model_rep_agent;

%% initialise variables for parameterisation

% grids depending on variation
zetaa_calib=1.4;
tauul_calib=0.181;
gri.zetaa= [zetaa_calib];
gri.tauul=[-0.2, 0, tauul_calib, 0.7];

%solution_LF=containers.Map;
%solution_Ramsey=containers.Map;

% number of periods for simulation 
T=101; % all time periods
P=60;  % periods for which to solve ramsey problem explicitly
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

%% - if want to compare different tauul; uncomment the following
for taut=1:length(gri.tauul) % loop over values for tauul
indic.tauul_ex=gri.tauul(taut);
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
end
%% calibration emissions and emission targets
% use this to calibrate relation of production and output
[targets_num, E_vec]=calibration_emissions(squeeze(sol_mat(list.y=='yd',1,gri.zetaa==1.4)), symms.targets, P);

%% grid method to maximise welfare function 

% takes ages!

%% primal and dynamic models
if indic.approach==1
    %- read in static problem
    [symms, list, Obj_ramPA]=primal_problem(y, x, list, symms, E);
    
    %- dynamic primal problems
    [Obj_ramPA_dynamic, symms.optim_dynamic, vecs, list ]=problem_dynamic(y, x, list, symms, E, Obj_ramPA, indic, pol, P);
else
    %- dynamic dual approach
     [Obj_ram_dynamic, symms.optim_dynamic, vecs, list]=problem_dynamic(y, x, list, symms, E, Obj_ram, indic, pol, P);
end
%% simulation Ramsey: static

%- initial values (guesses for ramsey problem updated in loop)
Ac1     = x_init(list.x=='Ac');
Ad1     = x_init(list.x=='Ad');

ramsey_solve_static;
%end
%end
%% dynamic problem
% static gives same results but easier to solve!
ramsey_solve_dynamic;

%% plots laissez-faire
 taul1=load(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
         T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
         params(list.params=='thetad'), indic.het_growth, gri.tauul(1), indic.util)...
         ,'sol_mat');
 taul2=load(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
         T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
         params(list.params=='thetad'), indic.het_growth, gri.tauul(2), indic.util)...
         ,'sol_mat');
taul3=load(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
         T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
         params(list.params=='thetad'), indic.het_growth, gri.tauul(3), indic.util)...
         ,'sol_mat');
taul4=load(sprintf('simulation_results/fullSimLF_T%d_initialAd%dAc%d_eppsilon%.2f_zetaa%.2f_thetac%.2f_thetad%.2f_HetGrowt%d_tauul%.3f_util%d.mat', ...
         T, Ad, Ac,params(list.params=='eppsilon'), params(list.params=='zetaa'), params(list.params=='thetac'), ...
         params(list.params=='thetad'), indic.het_growth, gri.tauul(4), indic.util)...
         ,'sol_mat');

list.plot=[list.y, list.x, "welfare"];
list.plot_mat=list.plot;
plottsLF1=taul1.sol_mat(:,1:T,params(list.params=='zetaa')==zetaa_calib);
plottsLF2=taul2.sol_mat(:,1:T,params(list.params=='zetaa')==zetaa_calib);
plottsLF3=taul3.sol_mat(:,1:T,params(list.params=='zetaa')==zetaa_calib);
plottsLF4=taul4.sol_mat(:,1:T,params(list.params=='zetaa')==zetaa_calib);

nn=5;

figure(1) %gcf=figure('Visible','off');
%         
        for i=1:length(list.plot)
        subplot(floor(length(list.plot)/nn)+1,nn,i)
        plot(time, plottsLF1(list.plot_mat==list.plot(i),:), time, plottsLF2(list.plot_mat==list.plot(i),:),...
             time, plottsLF3(list.plot_mat==list.plot(i),:), time, plottsLF4(list.plot_mat==list.plot(i),:), 'LineWidth', 1.3)
        legend(sprintf('tauul %.3f', gri.tauul(1)), sprintf('tauul %.3f', gri.tauul(2)), sprintf('tauul %.3f', gri.tauul(3)), sprintf('tauul %.3f', gri.tauul(4)), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
        ytickformat('%.2f')
        title(sprintf('%s', list.plot(i)), 'Interpreter', 'latex')
        end
        
%         subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
        plot(time, plottsLF1(list.plot_mat=='yd',:)./plottsLF1(list.plot_mat=='yc',:),...
            time, plottsLF2(list.plot_mat=='yd',:)./plottsLF2(list.plot_mat=='yc',:), ...
            time, plottsLF3(list.plot_mat=='yd',:)./plottsLF3(list.plot_mat=='yc',:), ...
            time, plottsLF4(list.plot_mat=='yd',:)./plottsLF4(list.plot_mat=='yc',:), 'LineWidth', 1.6)
        legend(sprintf('tauul %.3f', gri.tauul(1)), sprintf('tauul %.3f', gri.tauul(2)), sprintf('tauul %.3f', gri.tauul(3)), sprintf('tauul %.3f', gri.tauul(4)), ...
            'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
        ytickformat('%.2f')
            path=sprintf('figures/Rep_agent/comparison_tauul_ydyc_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_util%d_withtarget%d.png', T-1, ...
            params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad,Ac,...
            params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, indic.util, indic.withtarget);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
      %  close gcf
