function []=plottsSP(list, symms, T, params)

if ~isfile('figures/testfig')
    mkdir(sprintf('figures/all_0505'));
end
% this script plots results

syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An real
symms.plotsvarsProd =[Y N E G F];  
symms.plotsvarsHH =[hh hl C SWF Emnet];  
symms.plotsvarsRes =[sn sff sg  Af Ag An];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsPri =[pg pf pee pn wh wl ws];  

% symms.plotspol = [tauf taul taus lambdaa]; 
% symms.plotsprices = [pg pn pf pee wh wl ws];

listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri});
 
% read in results
helper=load(sprintf('LF_BAU_spillovers%d.mat', indic.spillovers));
bau=helper.LF_SIM;
helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers));
fb_lf=helper.LF_SIM;
helper=load(sprintf('SP_target_active_set_0505_spillover%d.mat', indic.spillovers));
sp_t=helper.sp_all';
helper=load(sprintf('SP_notarget_active_set_0505_spillover%d.mat', indic.spillovers));
sp_not=helper.sp_all';
helper=load(sprintf('OPT_notarget_active_set_0505_spillover%d.mat', indic.spillovers));
opt_not=helper.opt_all';
helper=load(sprintf('OPT_target_active_set_0505_spillover%d.mat', indic.spillovers));
opt_t=helper.opt_all';


RES = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT' , 'OPT_T', 'OPT_NOT'},...
                        {bau, fb_lf, sp_t, sp_not, opt_t, opt_not});

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
%     if ll~='VARS'
%         nn=2;
%     else 
    nn=3;
%     end
    %%% with subplots
    gcf=figure('Visible','off');
        
        for v=1:length(plotvars)
            varr=string(plotvars(v));
            subplot(floor(length(plotvars)/nn)+1,nn,v)
            plot(time,allvars(find(list.allvars==varr),:),  'LineWidth', 1.1)
            ytickformat('%.2f')
            xticklabels(Year)
            title(sprintf('%s', varr), 'Interpreter', 'latex')
        end
 
%         sgtitle('Social Planner Allocation')
        path=sprintf('figures/all_0505/%s_subplots_%s_spillover%d.png', ii, ll, indic.spillovers);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
  end
end      
%% separate plots        

for i= 1:length(listt)
    ss = listt(i);
    
        gcf=figure('Visible','off');
        pp=plot(time, plottsLF(listt.plot_mat==listt.plot(i),:), 'LineWidth', 1.6);
        set(pp, {'LineStyle'}, {'-'}, {'color'}, { 'k'}) 
        
%         if lgdd==1
%             legend(sprintf('with target' ), sprintf('without target'), 'Interpreter', 'latex', 'box', 'off', 'Location', 'best','FontSize', 20)
%         end
        ytickformat('%.2f')
        xticklabels(Year)
            ax=gca;
            ax.FontSize=13;
            xlabel('Periods', 'Fontsize', 20)
        path=sprintf('figures/Rep_agent/%s_onlyBAU_%s_periods%d_eppsilon%.2f_zeta%.2f_Ad0%d_Ac0%d_thetac%.2f_thetad%.2f_HetGrowth%d_util%d_lgd%d.png', prob, listt.plot(i), T-1, ...
            params(listt.params=='eppsilon'), params(listt.params=='zetaa'), Ad1,Ac1,...
            params(listt.params=='thetac'), params(listt.params=='thetad') , indic.het_growth, indic.util,  lgdd);
        exportgraphics(gcf,path,'Resolution', 400)
end
  
end
       