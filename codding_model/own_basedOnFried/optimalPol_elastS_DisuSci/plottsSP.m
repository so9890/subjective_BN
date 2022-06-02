function []=plottsSP(list, T, etaa, indic, params, Ems, plotts)

if ~isfile('figures/all_1705')
    mkdir(sprintf('figures/all_1705'));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.8 0.8 0.8];

if indic.ineq==0
    if indic.sep==1
        varlist=list.sepallvars;
    else
        varlist=list.allvars;
    end
else
    if indic.sep==1
        varlist=list.sepallvars_ineq;
    else
        varlist=list.allvars_ineq;
    end
end
% this script plots results

syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws wsg wsn wsf lambdaa Ch Cl C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An S real
%- additional vars
syms AgAf sgsff GFF EY CY hhhl whwl LgLf gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
if indic.ineq==0
    symms.plotsvarsHH =[hh hl C SWF Emnet]; 
else
    symms.plotsvarsHH =[hh hl Ch Cl  SWF Emnet];  
end
symms.plotsvarsRes =[sn sff sg  S Af Ag An];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [AgAf sgsff GFF EY CY hhhl whwl LgLf gAagg gAg gAf gAn Utilcon Utillab Utilsci ];

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
listt.plotsvarsAdd=string(symms.plotsvarsAdd);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%- variables which to plot in one graph plus legend
lissComp = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {string([hh hl]), string([sff sg sn S]), string([gAf gAg  gAn gAagg]), string([Lf Lg])});
legg = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {["high skill", "low skill"], ["fossil", "green", "non-energy", "total"], ["fossil", "green", "non-energy", "aggregate"], ["fossil", "green"]});
%% read in results
%- baseline results without reduction
helper=load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
bau=helper.LF_BAU';
% helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
% fb_lf=helper.LF_SIM';
helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,   etaa));
sp_t=helper.sp_all';
helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,  indic.BN_red,  etaa));
sp_not=helper.sp_all';
helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
opt_not_notaus=helper.opt_all';
helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
opt_t_notaus=helper.opt_all';

RES = containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                        {bau,  sp_t, sp_not, opt_t_notaus, opt_not_notaus});
%- add additional variables
[RES]=add_vars(RES, list, params, indic, varlist, symms);
varlist=[varlist, string(symms.plotsvarsAdd)];

%- results without taul
helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
opt_not_notaus_notaul=helper.opt_all';
helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
opt_t_notaus_notaul=helper.opt_all';


RES_polcomp_full   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus, opt_not_notaus});
RES_polcomp_notaul = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul});
%- add additional variables
RES_polcomp_full   =add_vars(RES_polcomp_full, list, params, indic, varlist, symms);
RES_polcomp_notaul =add_vars(RES_polcomp_notaul, list, params, indic, varlist, symms);

%- results with reduction 
% helper=load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red1_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
% bau_red=helper.LF_BAU';
% % helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
% % fb_lf=helper.LF_SIM';
% helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red1_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,   etaa));
% sp_t_red=helper.sp_all';
% helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red1_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,   etaa));
% sp_not_red=helper.sp_all';
% helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red1_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
% opt_not_notaus_red=helper.opt_all';
% helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red1_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
% opt_t_notaus_red=helper.opt_all';
% 
% RES_red = containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                         {bau_red,  sp_t_red, sp_not_red, opt_t_notaus_red, opt_not_notaus_red});
% 
% %- results without taul
% helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red1_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
% opt_not_notaus_notaul_red=helper.opt_all';
% helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red1_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, etaa));
% opt_t_notaus_notaul_red=helper.opt_all';
% 
% RES_polcomp_full_red   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_red, opt_not_notaus_red});
% RES_polcomp_notaul_red = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul_red, opt_not_notaus_notaul_red});
%            
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

% x indices
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));
%% Tables
betaa=params(list.params=='betaa');
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;
 %SWF_PV= zeros(length(keys(RES)),1);
 TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1));
%  TableSWF_PV_ns=table(keys(RES_ns)',zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1),zeros(length(keys(RES_ns)),1));
%  TableSWF_PV_ns.Properties.VariableNames={'Allocation','FullModel', 'NoSkillHet', 'NoTaul', 'NoTaulNoSkillHet'};
 TableSWF_PV.Properties.VariableNames={'Allocation','FullModel', 'NoTaul'};

%- all results
for i =keys(RES)
    ii=string(i);
    allvars= RES(ii);
    % SWF calculation 
    TableSWF_PV.FullModel(TableSWF_PV.Allocation==ii)=vec_discount*allvars(find(varlist=='SWF'),:)';
end
%- without taul
for i =keys(RES_polcomp_full)

     ii=string(i);
     allvars= RES_polcomp_full(ii);
     allvarsnt=RES_polcomp_notaul(ii); 
     TableSWF_PV.NoTaul(TableSWF_PV.Allocation==ii)=vec_discount*allvarsnt(find(varlist=='SWF'),:)';
end

TableSWF_PV.ContributionPERC=abs(TableSWF_PV.NoTaul-TableSWF_PV.FullModel)./abs(TableSWF_PV.FullModel)*100; 
save(sprintf('Table_SWF_sep%d_noskill%d_etaa%.2f_BN%d_ineq%d_red%d.mat', indic.sep, indic.noskill, etaa, indic.BN, indic.ineq, indic.BN_red), 'TableSWF_PV');

%% Plots
%% Subplots
% for i =keys(RES)
%     
%     %- loop
%     ii=string(i);
%     allvars= RES(ii);
% 
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
%             plot(time,allvars(find(varlist==varr),:),  'LineWidth', 1.1)
%             ytickformat('%.2f')
%             xticks(txx)
%             xticklabels(Year10)
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
%         end
%  
% %         sgtitle('Social Planner Allocation')
%         path=sprintf('figures/all_1705/%s_subplots_%s_spillover%d_etaa%.2f_sep%d_BN%d_ineq%d_red%d.png', ii, ll, indic.spillovers, etaa, indic.sep, indic.BN, indic.ineq, indic.BN_red);
%         exportgraphics(gcf,path,'Resolution', 400)
%         % saveas(gcf,path)
%         close gcf
%   end
% end     
%% All figures single
if plotts.single==1
    fprintf('plotting single graphs')
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
            if ll=="HH" && varr=="Emnet"
                main=plot(time,allvars(find(varlist==varr),:),time(3:end),Ems, 'LineWidth', 1.1);  
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

            else
                main=plot(time,allvars(find(varlist==varr),:), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

            end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
        path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f.png',  ii,varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
    %    pause
        end
    end
    end
end
    %% figures single overlayed
if plotts.singov==1
    fprintf('plotting single overlayed graphs')
    for lgdind=0:1
    for i =keys(RES)
        ii=string(i);

        %- loop
        ii=string(i);
        allvars= RES(ii);
    %% 
    fprintf('plotting %s',ii );
    for l =keys(lissComp) % loop over variable groups
        ll=string(l);
        plotvars=lissComp(ll); % here plotvars is a group of variable names which are to be plotted in the same graph

        gcf=figure('Visible','off');

        if length(plotvars)==2
            main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:), 'LineWidth', 1.1);    % plot vectors!        
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} ) 
    %    elseif ll=="LabourInp"
    %           main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),time,allvars(find(varlist==plotvars(3)),:), 'LineWidth', 1.1);    % plot vectors!        
    %           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; 'k'} )   
       elseif ll=="Growth"
              main=plot(time(1:end-1),allvars(find(varlist==plotvars(1)),1:end-1), time(1:end-1),allvars(find(varlist==plotvars(2)),1:end-1),...
              time(1:end-1),allvars(find(varlist==plotvars(3)),1:end-1),time(1:end-1),allvars(find(varlist==plotvars(4)),1:end-1),'LineWidth', 1.1);    % plot vectors!        
              set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   

        elseif ll=="Science"
              main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),...
                  time,allvars(find(varlist==plotvars(3)),:), time,allvars(find(varlist==plotvars(4)),:),'LineWidth', 1.1);    % plot vectors!        
               set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   
        end
           xticks(txx)
           if ll=="Growth"
              xlim([1, time(end-1)])
           else
              xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
         if lgdind==1
            pp=legg(ll);
            if length(pp)==2
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),  'Interpreter', 'latex');
            elseif length(pp)==3
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)),  'Interpreter', 'latex');

            elseif length(pp)==4
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)), sprintf('%s',pp(4)),  'Interpreter', 'latex');
            end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        end
        path=sprintf('figures/all_1705/SingleJointTOT_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png',  ii,ll, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
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
if plotts.notaul==1
    fprintf('plotting comparison no taul graphs')    
    for i =keys(RES_polcomp_full)

        ii=string(i);
             allvars= RES_polcomp_full(ii);
             allvarsnt=RES_polcomp_notaul(ii); 
%              TableSWF_PV.NoTaul(TableSWF_PV.Allocation==ii)=vec_discount*allvarsnt(find(varlist=='SWF'),:)';
%     % end
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
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
            path=sprintf('figures/all_1705/comp_notaul_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', ii, varr, indic.spillovers,indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            % saveas(gcf,path)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end      
end

%% comparison to BAU
if plotts.bau==1
    fprintf('plotting bau graphs')    
    bau=RES('BAU');

    for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
        ii=string(i);
        allvars= RES(ii);

    %% 
    fprintf('plotting %s',ii );
    for l = keys(lisst) % loop over variable groups
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

        path=sprintf('figures/all_1705/%s_BAUComp%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
    %     pause
        end
        end
      end
    end      
end
%% comparison with and without target    
%- string to loop over 
if plotts.comptarg==1
    fprintf('plotting comparison target graphs')    
    ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,3]
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));

%         if indic.noskill==0
             allvars= RES(t);
             allvarsnot=RES(nt); 
%         else
%              allvars= RES_ns(t);
%              allvarsnot=RES_ns(nt);
%         end

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
        path=sprintf('figures/all_1705/%s_TargetComp%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red,etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end
%% comparison social planner and optimal policy (with and without labour tax)
if plotts.compeff==1
    fprintf('plotting comparison efficient-optimal graphs')    
    eff= string({'SP_T', 'SP_NOT'});
    opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        allvars= RES(io);
        allvarsnotaul =RES_polcomp_notaul(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time,allvarseff(find(varlist==varr),:), time,allvars(find(varlist==varr),:), time,allvarsnotaul(find(varlist==varr),:), 'LineWidth', 1.2);            
           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'with income tax', ' no income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_1705/%s_CompEff%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', varr, io, indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red,etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end

    end
end

end

       