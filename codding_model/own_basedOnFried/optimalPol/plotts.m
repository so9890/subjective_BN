if ~isfile('figures/testfig')
    mkdir('figures');
end
% this script plots results

syms hh hl Y F E Emnet G pg pn pf pee tauf taul taus wh wl ws lambdaa C Lg Lf Ln xn xg xf sn sff sg real
symms.plotsvars =[hh hl Y F E G  C xn xg xf sn sff sg Emnet]; % LG Lf LN missing in allvars! 
symms.plotspol = [tauf taul taus lambdaa]; 
symms.plotsprices = [pg pn pf pee wh wl ws];

list.plotsvars=string(symms.plotsvars);
list.plotspol=string(symms.plotspol);
list.plotsprices=string(symms.plotsprices);

% read in results
BAU=load('LF_BAU.mat');
FB_LF=load('FB_LF_SIM_NOTARGET.mat');
SP_T=load('SP_target.mat');
SP_NOT=load('SP_notarget.mat');

% x indices
Year =transpose(year(['2025'; '2030';'2035'; '2040';'2045'; '2050';'2055'; '2060'; '2065';'2070';'2075';'2080'],'yyyy'));
time = 1:T;
% number of figures in row in subplot
nn=5;
% plot
%allvars=BAU.LF_SIM'; % transpose in LF case so that same order as SP : each column a different period
allvars= SP_T.sp_all';
lisst=list.plotspol;
% plot(time, allvars(find(list.allvars=='sff'),:))
        %%% with subplots
          gcf=figure('Visible','off');
        
        for i=1:length(lisst)
            varr=string(lisst(i));
        subplot(floor(length(lisst)/nn)+1,nn,i)
        plot(time,allvars(find(list.allvars==varr),:),  'LineWidth', 1.3)
        ytickformat('%.2f')
        xticklabels(Year)
        title(sprintf('%s', varr), 'Interpreter', 'latex')
        end
        
%         subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+1)
%         plot(time, plottsLF(list.plot_mat=='yd',:)./plottsLF(list.plot_mat=='yc',:), time, plottsRam(list.plot_mat=='yd',:)./plottsRam(list.plot_mat=='yc',:), 'LineWidth', 1.6)
%         %legend(sprintf('LF'), sprintf('Ramsey'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best');
%         ytickformat('%.2f')
%         xticklabels(Year)
%         title('$y_d/y_c$', 'Interpreter', 'latex')
%         %set(lgd, 'Interpreter', 'latex', 'box', 'off', 'Location', 'best')
%         
%         subplot(floor(length(list.plot)/nn)+1,nn,length(list.plot)+2)
%         plot(time, opt_pol_simRam, 'LineWidth', 1.6)
%         title(sprintf('Optimal $\\tau_l$'), 'Interpreter', 'latex')%, 'box', 'off', 'Location', 'best')
%         ytickformat('%.2f')
%         xticklabels(Year)
        sgtitle('Social Planner Allocation')
        path=sprintf('figures/sp_subplots_variables.png');
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
        
%% separate plots        

for i= 1:length(lisst)
    ss = lisst(i);
    
        gcf=figure('Visible','off');
        pp=plot(time, plottsLF(list.plot_mat==list.plot(i),:), 'LineWidth', 1.6);
        set(pp, {'LineStyle'}, {'-'}, {'color'}, { 'k'}) 
        
%         if lgdd==1
%             legend(sprintf('with target' ), sprintf('without target'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best','FontSize', 20)
%         end
        ytickformat('%.2f')
        xticklabels(Year)
            ax=gca;
            ax.FontSize=13;
            xlabel('Periods', 'Fontsize', 20)
        path=sprintf('figures/Rep_agent/%s_onlyBAU_%s_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_util%d_lgd%d.png', prob, list.plot(i), T-1, ...
            params(list.params=='eppsilon'), params(list.params=='zetaa'), Ad1,Ac1,...
            params(list.params=='thetac'), params(list.params=='thetad') , indic.het_growth, indic.util,  lgdd);
        exportgraphics(gcf,path,'Resolution', 400)
end
  

       