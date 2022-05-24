function []=plottsSP(list, T, etaa, indic, params)

if ~isfile('figures/testfig')
    mkdir(sprintf('figures/all_1705'));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];

if indic.sep==1
    varlist=list.sepallvars;
else
    varlist=list.allvars;
end
% this script plots results

syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws wsg wsn wsf lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An real
symms.plotsvarsProd =[Y N E G F];  
symms.plotsvarsHH =[hh hl C SWF Emnet];  
symms.plotsvarsRes =[sn sff sg  Af Ag An];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
if indic.sep==0
    symms.plotsvarsPri =[pg pf pee pn wh wl ws];  
else
    symms.plotsvarsPri =[pg pf pee pn wh wl wsg wsn wsf];  
end

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
 
%% read in results
helper=load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
bau=helper.LF_BAU';
% helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
% fb_lf=helper.LF_SIM';
helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,   etaa));
sp_t=helper.sp_all';
helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,   etaa));
sp_not=helper.sp_all';
helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
opt_not_notaus=helper.opt_all';
helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
opt_t_notaus=helper.opt_all';

RES = containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                        {bau,  sp_t, sp_not, opt_t_notaus, opt_not_notaus});

%- results without taul
helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
opt_not_notaus_notaul=helper.opt_all';
helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
opt_t_notaus_notaul=helper.opt_all';


RES_polcomp_full   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus, opt_not_notaus});
RES_polcomp_notaul = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul});
                    
% %-- results with spillover 
% 
% % read in results
% helper=load(sprintf('LF_BAU_spillovers1_noskill0_sep0_etaa1.20.mat'));
% bau=helper.LF_BAU';
% % helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
% % fb_lf=helper.LF_SIM';
% helper=load(sprintf('SP_target_active_set_1705_spillover1_noskill0_sep0_etaa1.20.mat'));
% sp_t=helper.sp_all';
% helper=load(sprintf('SP_notarget_active_set_1705_spillover1_noskill0_sep0_etaa1.20.mat'));
% sp_not=helper.sp_all';
% helper=load(sprintf('OPT_notarget_active_set_1905_spillover1_taus0_noskill0_notaul0_sep0_etaa1.20.mat'));
% opt_not_notaus=helper.opt_all';
% helper=load(sprintf('OPT_target_active_set_0505_spillover1_taus0_noskill0_notaul0_sep0_etaa1.20.mat'));
% opt_t_notaus=helper.opt_all';
% 
% RES_spill= containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                         {bau,  sp_t, sp_not, opt_t_notaus, opt_not_notaus});
%% results without skill heterogeneity
% helper=load(sprintf('LF_BAU_spillovers%d_noskill1.mat', indic.spillovers));
% bau_ns=helper.LF_BAU';
% helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill1.mat', indic.spillovers));
% fb_lf_ns=helper.LF_SIM';
% helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill1.mat', indic.spillovers));
% sp_t_ns=helper.sp_all';
% helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill1.mat', indic.spillovers));
% sp_not_ns=helper.sp_all';
% % helper=load(sprintf('OPT_notarget_active_set_0505_spillover%d.mat', indic.spillovers));
% % opt_not=helper.opt_all';
% % helper=load(sprintf('OPT_target_active_set_0505_spillover%d.mat', indic.spillovers));
% % opt_t=helper.opt_all';
% % results without taus
% helper=load(sprintf('OPT_notarget_active_set_0505_spillover%d_taus%d_noskill1_notaul0_alt.mat', indic.spillovers, indic.taus));
% opt_not_notaus_ns=helper.opt_all';
% helper=load(sprintf('OPT_target_active_set_0505_spillover%d_taus%d_noskill1_notaul0.mat', indic.spillovers, indic.taus));
% opt_t_notaus_ns=helper.opt_all';
% 
% RES_ns = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                         {bau_ns, fb_lf_ns, sp_t_ns, sp_not_ns, opt_t_notaus_ns, opt_not_notaus_ns});

%% optimal policy scenarios

% RES_polcomp_full_ns   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_ns, opt_not_notaus_ns});
% RES_polcomp_notaul_ns = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul_ns, opt_not_notaus_notaul_ns});

%% SWF comparison
betaa=params(list.params=='betaa');
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;
 %SWF_PV= zeros(length(keys(RES)),1);
 TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1));
%  TableSWF_PV_ns=table(keys(RES_ns)',zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1));
%  TableSWF_PV_ns.Properties.VariableNames={'Allocation','FullModel', 'NoSkillHet', 'NoTaul', 'NoTaulNoSkillHet'};
 TableSWF_PV.Properties.VariableNames={'Allocation','FullModel', 'NoTaul'};
 TableSWF_PV_spill=TableSWF_PV;

% x indices
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%% Tables
for i =keys(RES)
    ii=string(i);
    allvars= RES(ii);
    % SWF calculation 
    TableSWF_PV.FullModel(TableSWF_PV.Allocation==ii)=vec_discount*allvars(find(varlist=='SWF'),:)';
end

%% plot: Subplots
for i =keys(RES)
    
    %- loop
    ii=string(i);
    allvars= RES(ii);

%% 
fprintf('plotting %s',ii );
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
            plot(time,allvars(find(varlist==varr),:),  'LineWidth', 1.1)
            ytickformat('%.2f')
            xticks(txx)
            xticklabels(Year10)
            title(sprintf('%s', varr), 'Interpreter', 'latex')
        end
 
%         sgtitle('Social Planner Allocation')
        path=sprintf('figures/all_1705/%s_subplots_%s_spillover%d_etaa%.2f_sep%d_BN%d_ineq%d.png', ii, ll, indic.spillovers, etaa, indic.sep, indic.BN, indic.ineq);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
  end
end     
%% All figures single
% RES = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                        {bau, fb_lf, sp_t, sp_not, opt_t_notaus, opt_not_notaus});

for i =keys(RES)
    ii=string(i);

    %- loop
    ii=string(i);
    allvars= RES(ii);
%% 
fprintf('plotting %s',ii );
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);

    for v=1:length(plotvars)
    gcf=figure('Visible','off');
        

        varr=string(plotvars(v));
        %subplot(floor(length(plotvars)/nn)+1,nn,v)
        main=plot(time,allvars(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
    path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_sep%d_BN%d_ineq%d_etaa%.2f.png',  ii,varr, indic.spillovers, indic.sep, indic.BN, indic.ineq, etaa);
    exportgraphics(gcf,path,'Resolution', 400)
    % saveas(gcf,path)
%    close gcf
%    pause
    end
end
end

%% comparison without skill heterogeneity
% for i =keys(RES)
%     ii=string(i);
%     allvars= RES(ii);
%     allvarsns=RES_ns(ii);
%     % SEF calculation 
%     TableSWF_PV_ns.FullModel(TableSWF_PV_ns.Allocation==ii)=vec_discount*allvars(find(varlist=='SWF'),:)';
%     TableSWF_PV_ns.NoSkillHet(TableSWF_PV_ns.Allocation==ii)=vec_discount*allvarsns(find(varlist=='SWF'),:)';
% % end
% %%
% fprintf('plotting %s',ii );
% for l =keys(lisst)
%     ll=string(l);
%     plotvars=lisst(ll);
%     % number of figures in row in subplot
% %     if ll~='VARS'
% %         nn=2;
% %     else 
%     nn=3;
% %     end
%     %%% with subplots
%     gcf=figure('Visible','off');
%         
%         for v=1:length(plotvars)
%             varr=string(plotvars(v));
%             subplot(floor(length(plotvars)/nn)+1,nn,v)
%             main=plot(time,allvars(find(varlist==varr),:), time,allvarsns(find(varlist==varr),:), 'LineWidth', 1.1);
%             
%            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
%            xticks(time)
%             ax=gca;
%             ax.FontSize=13;
% %             xlim([0.1, 0.9])
%             ytickformat('%.2f')
%             xticklabels(Year)
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
%         end
%  
% %         sgtitle('Social Planner Allocation')
%         path=sprintf('figures/all_0505/%s_subplots_%s_spillover%d_comp_noskill.png', ii, ll, indic.spillovers);
%         exportgraphics(gcf,path,'Resolution', 400)
%         % saveas(gcf,path)
%         close gcf
%   end
% end      


%% comparison with and without taul

for i =keys(RES_polcomp_full)
    
%     if indic.sep==0
%         lisst('Pri')= [pg pf pee pn wh wl ws];  
%     else
%         lisst('Pri')= [pg pf pee pn wh wl wsg wsn wsf];  
%     end
    
    
    ii=string(i);
%     if indic.noskill==0
         allvars= RES_polcomp_full(ii);
         allvarsnt=RES_polcomp_notaul(ii); 
         TableSWF_PV.NoTaul(TableSWF_PV.Allocation==ii)=vec_discount*allvarsnt(find(varlist=='SWF'),:)';

%     else
%         allvars= RES_polcomp_full_ns(ii);
%         allvarsnt=RES_polcomp_notaul_ns(ii);
%         TableSWF_PV_ns.NoTaulNoSkillHet(TableSWF_PV_ns.Allocation==ii)=vec_discount*allvarsnt(find(varlist=='SWF'),:)';
%     end

%% 
fprintf('plotting %s',ii );
for lgdind=1
for l =keys(lisst)
    ll=string(l);
    plotvars=lisst(ll);
    % number of figures in row in subplot
%     if ll~='VARS'
%         nn=2;
%     else 
%     nn=3;
%     end
    %%% with subplots
       for v=1:length(plotvars)
       gcf=figure('Visible','off');
        
     
            varr=string(plotvars(v));
%             subplot(floor(length(plotvars)/nn)+1,nn,v)
            main=plot(time,allvars(find(varlist==varr),:), time,allvarsnt(find(varlist==varr),:), 'LineWidth', 1.1);
            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
      
        if lgdind==1
          
        lgd=legend('full model', 'no income tax', 'Interpreter', 'latex');
           
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        end
%         sgtitle('Social Planner Allocation')
        path=sprintf('figures/all_1705/comp_notaul_%s_%s_spillover%d_sep%d_BN%d_ineq%d_etaa%.2f_lgd%d.png', ii, varr, indic.spillovers, indic.sep, indic.BN, indic.ineq, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
       end % variables in group
end % variable group
end % legend
end      
 % calculate contributes to SWF 
% tt=load('Table_SWF_etaa1.mat', 'TableSWF_PV');
TableSWF_PV.ContributionPERC=abs(TableSWF_PV.NoTaul-TableSWF_PV.FullModel)./abs(TableSWF_PV.FullModel)*100; 
%save(sprintf('Table_SWF_sep%d_etaa%.2f_BN%d_ineq%d.mat', indic.sep, etaa, indic.BN, indic.ineq), 'TableSWF_PV');

%% comparison to BAU
% RES = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                        {bau, fb_lf, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
for j=0
    indic.noskill=j;
    if indic.noskill==0
        bau=RES('BAU');
    else
        bau=RES_ns('BAU');
    end
for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
    ii=string(i);

    if indic.noskill==0
         allvars= RES(ii);
    else
         allvars= RES_ns(ii);
    end

%% 
fprintf('plotting %s',ii );
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
    gcf=figure('Visible','off');
        

        varr=string(plotvars(v));
        %subplot(floor(length(plotvars)/nn)+1,nn,v)
        main=plot(time,allvars(find(varlist==varr),:), time,bau(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
        title(sprintf('%s', varr), 'Interpreter', 'latex')
       if lgdind==1
           if contains(ii, 'SP')
              lgd=legend('social planner', 'bau', 'Interpreter', 'latex');
           else
              lgd=legend('ramsey planner', 'bau', 'Interpreter', 'latex');
           end
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
       end

    path=sprintf('figures/all_1705/%s_BAUComp%s_spillover%d_sep%d_BN%d_ineq%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.sep, indic.BN, indic.ineq, etaa, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
    % saveas(gcf,path)
%    close gcf
%     pause
    end
    end
  end
end      
end
%% comparison with and without target    
%- string to loop over 
ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
for j=0
    indic.noskill=j;
for i =[1,3]
    ii=ssr(i);
    %- read in data
    t=string(ssr(i));
    nt=string(ssr(i+1));
    
    if indic.noskill==0
         allvars= RES(t);
         allvarsnot=RES(nt); 
    else
         allvars= RES_ns(t);
         allvarsnot=RES_ns(nt);
    end

for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
        gcf=figure('Visible','off');
        varr=string(plotvars(v));
        main=plot(time,allvars(find(varlist==varr),:), time,allvarsnot(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
       if lgdind==1
          lgd=legend('wih emission target', 'no emission target', 'Interpreter', 'latex');
          set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
       end
    path=sprintf('figures/all_1705/%s_TargetComp%s_spillover%d_sep%d_BN%d_ineq%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.sep,indic.BN, indic.ineq,etaa, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
    % saveas(gcf,path)
%    close gcf
%     pause
    end
    end
  end
end      
end

%% comparison social planner and optimal policy (with and without labour tax)
eff= string({'SP_T', 'SP_NOT'});
opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

for withtaul=0:1
for i =[1,2]
    
    ie=eff(i);
    io=opt(i);
    
     if withtaul==1
        allvars= RES(io);
     else
         allvars =RES_polcomp_notaul(io);
     end
     allvarseff=RES(ie); 
    
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
        gcf=figure('Visible','off');
        varr=string(plotvars(v));
        main=plot(time,allvars(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
       if lgdind==1
          lgd=legend('optimal allocation', 'first best allocation', 'Interpreter', 'latex');
          set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
       end
    path=sprintf('figures/all_1705/%s_CompEff%s_spillover%d_sep%d_BN%d_ineq%d_etaa%.2f_withtaul%d_lgd%d.png', varr, ii, indic.spillovers, indic.sep,indic.BN, indic.ineq,etaa, withtaul, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
    % saveas(gcf,path)
%    close gcf
%     pause
    end
    end
  end

end

end
end
       