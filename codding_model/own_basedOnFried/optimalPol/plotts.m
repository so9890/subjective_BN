function []=plotts(list, symms, T, params)

if ~isfile('figures/testfig')
    mkdir('figures');
end
% this script plots results

syms hh hl Y F E Emnet G pg pn pf pee tauf taul taus wh wl ws lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An real
symms.plotsvars =[hh hl Y F E G  C xn xg xf sn sff sg Emnet Ln Lg Lf Af Ag An SWF];  
symms.plotspol = [tauf taul taus lambdaa]; 
symms.plotsprices = [pg pn pf pee wh wl ws];

list.plotsvars=string(symms.plotsvars);
list.plotspol=string(symms.plotspol);
list.plotsprices=string(symms.plotsprices);

lisst = containers.Map({'VARS', 'POL', 'PRICES'}, {list.plotsvars, list.plotspol, list.plotsprices});
 
% read in results
helper=load('LF_BAU.mat');
bau=helper.LF_SIM;
helper=load('FB_LF_SIM_NOTARGET.mat');
fb_lf=helper.LF_SIM;
helper=load('SP_target.mat');
sp_t=helper.sp_all';
helper=load('SP_notarget.mat');
sp_not=helper.sp_all';
helper=load('OPT_notarget.mat');
opt_not=helper.opt_all';

RES = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT', 'OPT_NOT'}, {bau, fb_lf, sp_t, sp_not, opt_not});

% SWF comparison
betaa=params(list.params=='betaa');
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;
 SWF_PV= zeros(length(keys(RES)),1);
 
% x indices
Year =transpose(year(['2025'; '2030';'2035'; '2040';'2045'; '2050';'2055'; '2060'; '2065';'2070';'2075';'2080'],'yyyy'));
time = 1:T;

%% plot
for i =keys(RES)
    ii=string(i);
    fprintf('plotting %s',ii );
    allvars= RES(ii);
    % SEF calculation 
    SWF_PV(keys(RES)==ii)=vec_discount*allvars(find(list.allvars=='SWF'),:)';

%% 
for l =keys(lisst)
    ll=string(l);
    plotvars=lisst(ll);
    % number of figures in row in subplot
    if ll~='VARS'
        nn=2;
    else 
        nn=4;
    end
    %%% with subplots
    gcf=figure; %('Visible','off');
        
        for v=1:length(plotvars)
            varr=string(plotvars(v));
            subplot(floor(length(plotvars)/nn)+1,nn,v)
            plot(time,allvars(find(list.allvars==varr),:),  'LineWidth', 1.1)
            ytickformat('%.2f')
            xticklabels(Year)
            title(sprintf('%s', varr), 'Interpreter', 'latex')
        end
 
%         sgtitle('Social Planner Allocation')
        path=sprintf('figures/%s_subplots_%s_new.png', ii, ll);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
  end
end      
%% separate plots        

for i= 1:length(list)
    ss = list(i);
    
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
  
end
       