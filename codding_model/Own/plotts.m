%% Plots
%% Ramsey versus laissez faire

for prob=["static"]% "dynamic"]
for ss=0:1
 
    indic.subst= ss;
  
    if ss==0
        params(list.params=='eppsilon')=epps(1);
    else
       params(list.params=='eppsilon')=epps(2);
    end

for ttt=0:1
  indic.withtarget=ttt;

if prob == "static"
% choose whether to use static or dynamic problem
    load(sprintf('simulation_results/StaticControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
    load(sprintf('simulation_results/StaticStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
else
     load(sprintf('simulation_results/DynamicControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
     load(sprintf('simulation_results/DynamicStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
end

% if dynamic
if prob == "dynamic"
    P=31;
else
    P=60;
end
time=1:P;
T=P;

% number of figures in row in subplot
nn=5;
% list of variables to be plotted
%list.plot=[list.y(~ismember(list.y, ["pcL" "pdL" "lhc" "llc" "lhd"
%"lld"])), list.x, "welfare"];
list.plot=[list.y list.x, "welfare"];

% combine variables into one matrix and list 

plottsLF=sol_mat(:,1:T,params(list.params=='zetaa')==zetaa_calib);

% welfare
welf_sim=log(y_simRam(list.y=='c',:))-(y_simRam(list.y=='hl',:)+...
        params(list.params=='zetaa').*y_simRam(list.y=='hh',:)).^(1+params(list.params=='sigmaa'))./(1+params(list.params=='sigmaa'));

plottsRam=[y_simRam(:,1:T); x_simRam(:,1:T); welf_sim(:,1:T) ];
% optimal policy
opt_pol_sim= 1-y_simRam(list.y=='H',1:T).^(1+params(list.params=='sigmaa'));

% list as in matrices of results
list.plot_mat=[list.y, list.x, "welfare"];

%% figures
%- all in one figure
gcf=figure('Visible','off');

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plottsLF(list.plot_mat==list.plot(i),:), time, plottsRam(list.plot_mat==list.plot(i),:), 'LineWidth', 1.3)
%legend(sprintf('LF'), sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
ytickformat('%.2f')
title(sprintf('%s', list.plot(i)), 'Interpreter', 'latex')
end

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
plot(time, plottsLF(list.plot_mat=='yd',:)./plottsLF(list.plot_mat=='yc',:), time, plottsRam(list.plot_mat=='yd',:)./plottsRam(list.plot_mat=='yc',:), 'LineWidth', 1.6)
%legend(sprintf('LF'), sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
ytickformat('%.2f')
title('$y_d/y_c$', 'Interpreter', 'latex')
%set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+2)
plot(time, opt_pol_sim, 'LineWidth', 1.6)
title(sprintf('Optimal $\\tau_l$'), 'Interpreter', 'latex')%, 'box', 'off', 'Location', 'best')
ytickformat('%.2f')

sgtitle('Laissez Faire versus Ramsey')
path=sprintf('figures/Rep_agent/%sRam_LF_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d_lgd0.png', prob, T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget);
exportgraphics(gcf,path,'Resolution', 400)
% saveas(gcf,path)
close gcf

%- separate figures 
for lgdd=0:1
for i=1:length(list.plot)
gcf=figure('Visible','off');
pp=plot(time, plottsLF(list.plot_mat==list.plot(i),:), time, plottsRam(list.plot_mat==list.plot(i),:), 'LineWidth', 1.6);
set(pp, {'LineStyle'}, {'--'; '-'}, {'color'}, {'k'; 'k'})  
if lgdd==1
    legend(sprintf('LF'), sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
end
ytickformat('%.2f')

path=sprintf('figures/Rep_agent/%sRam_LF_separate_%s_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d_lgd%d.png', prob, list.plot(i), T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget, lgdd);
exportgraphics(gcf,path,'Resolution', 400)
end


gcf=figure('Visible','off');
plot(time, plottsLF(list.plot_mat=='yd',:)./plottsLF(list.plot_mat=='yc',:), time, plottsRam(list.plot_mat=='yd',:)./plottsRam(list.plot_mat=='yc',:), 'LineWidth', 1.6)
set(pp, {'LineStyle'}, {'--'; '-'}, {'color'}, {'k'; 'k'})  
if lgdd==1
    legend(sprintf('LF'), sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
end
ytickformat('%.2f')
title('$y_d/y_c$', 'Interpreter', 'latex')
%set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
path=sprintf('figures/Rep_agent/%sRam_LF_separate_ydyc_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d_lgd%d.png', prob, T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget, lgdd);
exportgraphics(gcf,path,'Resolution', 400)


gcf=figure('Visible','off');
pp=plot(time, opt_pol_sim, 'LineWidth', 1.6);
set(pp, {'LineStyle'}, {'-'}, {'color'}, {'k'})  

if lgdd==1
    legend( sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
end
ytickformat('%.2f')

path=sprintf('figures/Rep_agent/%sRam_LF_separate_tauul_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d_lgd%d.png', prob, T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget, lgdd);
exportgraphics(gcf,path,'Resolution', 400)
close gcf

end
end
end
end
%% test: static versus 
for ss=0:1
    for ttt=0:1
indic.subst=ss;

indic.withtarget=ttt;
dync  = load(sprintf('simulation_results/DynamicControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
dyns  = load(sprintf('simulation_results/DynamicStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
statc = load(sprintf('simulation_results/StaticControlsRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'y_simRam');
stats =  load(sprintf('simulation_results/StaticStatesRamsey_hetgrowth%d_util%d_withtarget%d_eppsilon%.2f_dual%d_zetaa%.2f_thetac%.2f_thetad%.2f_initialAd%dAc%d.mat',...
        indic.het_growth, indic.util, indic.withtarget, params(list.params=="eppsilon"), indic.approach, params(list.params=='zetaa'), params(list.params=='thetac'), ...
        params(list.params=='thetad'), Ad, Ac ),'x_simRam');
% if dynamic
P=31;
time=1:P;
T=P;
nn=5;
%welf_sim=log(y_simRam(list.y=='c',:))-(y_simRam(list.y=='hl',:)+...
 %       params(list.params=='zetaa').*y_simRam(list.y=='hh',:)).^(1+params(list.params=='sigmaa'))./(1+params(list.params=='sigmaa'));

plottsRam_dyn=[dync.y_simRam(:,1:T); dyns.x_simRam(:,1:T) ];
plottsRam_stat=[statc.y_simRam(:,1:T); stats.x_simRam(:,1:T) ];

% list as in matrices of results
list.plot_mat=[list.y, list.x];
% list of variables to be plotted
list.plot=[list.y, list.x];

figure(1)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plottsRam_stat(list.plot_mat==list.plot(i),:), time, plottsRam_dyn(list.plot_mat==list.plot(i),:), 'LineWidth', 1.6)
legend(sprintf('Static: %s', list.plot(i)), sprintf('Dynamic: %s', list.plot(i)), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
ytickformat('%.2f')
end

subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
plot(time, plottsRam_stat(list.plot_mat=='yd',:)./plottsRam_stat(list.plot_mat=='yc',:), time, plottsRam_dyn(list.plot_mat=='yd',:)./plottsRam_dyn(list.plot_mat=='yc',:), 'LineWidth', 1.6)
legend(sprintf('Static: $y_d/y_c$'), sprintf('Dynamic: $y_d/y_c$'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
ytickformat('%.2f')
%set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')


sgtitle('Laissez Faire versus Ramsey')
path=sprintf('figures/Rep_agent/dynamic_vs_static_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_tauul%.3f_util%d_withtarget%d.png', T-1, ...
    params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
    params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, pols_num(list.pol=='tauul'), indic.util, indic.withtarget);
saveas(gcf,path)

close gcf 
    end
end
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

plottss=[y_simLF;x_simLF];
list.joint=[list.y, list.x];

figure(1)

for i=1:length(list.plot)
subplot(floor(length(list.plot)/nn)+1,nn,i)
plot(time, plottss(list.joint==list.plot(i),:))
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
