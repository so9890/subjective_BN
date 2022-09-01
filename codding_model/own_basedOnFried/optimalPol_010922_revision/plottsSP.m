function []=plottsSP(list, T, etaa, weightext,indic, params, Ems, plotts)

date="Aout22";
if ~isfile(sprintf('figures/all_%s', date ))
    mkdir(sprintf('figures/all_%s', date));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.6 0.6 0.6];

if indic.sep==1
    varlist=[list.sepallvars];
else
    varlist=[list.allvars];
end
% this script plots results

syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws wsg wsn wsf lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An A S real
%- additional vars
syms analyTaul PV CEVv CEVvPV CEVvDy AgAf sgsff GFF EY CY hhhl whwl LgLf gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg  S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul PV AgAf sgsff GFF EY CY hhhl whwl LgLf gAagg gAg gAf gAn Utilcon Utillab Utilsci];
% already exists: symms.addgov
symms.comp=[ CEVv CEVvDy CEVvPV ]; % for comparison of policy interventions, 

if indic.sep==0
    symms.plotsvarsPri =[pg pf pee pn wh wl ws];  
else
    symms.plotsvarsPri =[pg pf pee pn wh wl wsg wsn wsf];  
end

listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);
listt.plotsvarsAdd=string(symms.plotsvarsAdd);
list.comp=string(symms.comp);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%- variables which to plot in one graph plus legend
lissComp = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {string([hh hl]), string([sff sg sn S]), string([gAf gAg  gAn gAagg]), string([Lf Lg])});
legg = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {["high skill", "low skill"], ["fossil", "green", "non-energy", "total"], ["fossil", "green", "non-energy", "aggregate"], ["fossil", "green"]});

%% read in results
%- baseline results without reduction
if indic.xgrowth==0
    helper=load(sprintf('BAU_spillovers%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, etaa));
    bau=helper.LF_BAU';
    helper=load(sprintf('LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill,indic.sep, etaa));
    LF=helper.LF_SIM';
else
    helper=load(sprintf('BAU_xgrowth_spillovers%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, etaa));
    bau=helper.LF_SIM;
    helper=load(sprintf('LF_xgrowth_spillovers%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, etaa));
    LF=helper.LF_SIM;
end
helper=load(sprintf('SP_target_0508_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,  indic.xgrowth, indic.PV, etaa));
sp_t=helper.sp_all';
helper=load(sprintf('SP_notarget_0508_spillover%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,  indic.PV,etaa));
sp_not=helper.sp_all';
helper=load(sprintf('OPT_notarget_0308_spillover%d_taus0_noskill%d_notaul0_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,indic.PV, etaa));
opt_not_notaus=helper.opt_all';
helper=load(sprintf('OPT_target_0308_spillover%d_taus0_noskill%d_notaul0_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,indic.PV, etaa));
opt_t_notaus=helper.opt_all';

RES = containers.Map({'BAU','LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                        {bau,  LF, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
%- add additional variables
[RES]=add_vars(RES, list, params, indic, varlist, symms);
varlist=[varlist, string(symms.plotsvarsAdd)];

%-results with counterfactual productivity gap
if indic.noskill==0
    if indic.xgrowth==0 
        helper=load(sprintf('BAU_countec_spillovers%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,etaa));
        bau=helper.LF_SIM;
        helper=load(sprintf('LF_countec_spillovers%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, etaa));
        LF=helper.LF_SIM;
    end

%     helper=load(sprintf('SP_target_active_set_1705_countec_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,  indic.xgrowth, etaa));
%     sp_t=helper.sp_all';
%     helper=load(sprintf('SP_notarget_active_set_1705_countec_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern0_xgrowth%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,  indic.BN_red, indic.xgrowth,  etaa));
%     sp_not=helper.sp_all';
%     helper=load(sprintf('OPT_notarget_active_set_1905_countec_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_extern0_xgrowth%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, etaa));
%     opt_not_notaus=helper.opt_all';
%     helper=load(sprintf('OPT_target_active_set_1905_countec_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, etaa));
%     opt_t_notaus=helper.opt_all';
% 
%     RES_countec = containers.Map({'BAU','LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                             {bau,  LF, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
% 
%     [RES_countec]=add_vars(RES_countec, list, params, indic, varlist, symms);
end

%- version without endogenous growth
%- for comparison 
helper=load(sprintf('SP_target_0508_spillover%d_noskill%d_sep%d_xgrowth1_PV%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.PV, etaa));
sp_t=helper.sp_all';
helper=load(sprintf('SP_notarget_0508_spillover%d_noskill%d_sep%d_extern0_xgrowth1_PV%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep,indic.PV, etaa));
sp_not=helper.sp_all';
helper=load(sprintf('OPT_notarget_0308_spillover%d_taus0_noskill%d_notaul0_sep%d_extern0_xgrowth1_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep,indic.PV, etaa));
opt_not_notaus=helper.opt_all';
helper=load(sprintf('OPT_target_0308_spillover%d_taus0_noskill%d_notaul0_sep%d_xgrowth1_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep,indic.PV, etaa));
opt_t_notaus=helper.opt_all';

RES_xgr = containers.Map({ 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                        { sp_t, sp_not, opt_t_notaus, opt_not_notaus});

[RES_xgr]=add_vars(RES_xgr, list, params, indic, varlist, symms);

%%
%- results with other policy specifications
OTHERPOL={}; % cell of containers over which to loop
% add other government variables here!
% read in again because variables have changed!

% varlist: first separate  variables, then additional gov variables, then
% new variables from function 
if indic.sep==1
    varlist_polcomp=[list.sepallvars, list.addgov, string(symms.plotsvarsAdd)];
else
    varlist_polcomp=[list.allvars, list.addgov, string(symms.plotsvarsAdd)];
end

for i=0:5 % loop over policy versions
    if indic.sep==0
        error('have to change variable list in function add_vars below')
    end
    helper=load(sprintf('OPT_notarget_0308_spillover%d_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, i,indic.sep,indic.xgrowth,indic.PV, etaa));
    opt_not_notaus_notaul=[helper.opt_all';helper.addGov'];
    helper=load(sprintf('OPT_target_0308_spillover%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill,i, indic.sep, indic.xgrowth, indic.PV, etaa));
    opt_t_notaus_notaul=[helper.opt_all';helper.addGov'];
        
    if i==0
        RES_polcomp_notaul0 = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{opt_t_notaus_notaul, opt_not_notaus_notaul});
        if indic.sep==1
            RES_polcomp_notaul0 = add_vars(RES_polcomp_notaul0, list, params, indic, list.sepallvars, symms);
        else
           RES_polcomp_notaul0 = add_vars(RES_polcomp_notaul0, list, params, indic, list.allvars, symms);
        end
    else
        % add additional variables and save container to cell
       if indic.sep==1
            OTHERPOL{i} = add_vars(containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul}), list, params, indic, list.sepallvars, symms);
       else
            OTHERPOL{i} = add_vars(containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul}), list, params, indic, list.allvars, symms);

       end
       end
end

%%
% - results with externality
% helper=load(sprintf('SP_notarget_spillover%d_noskill%d_sep%d_extern1_weightext0.01_xgrowth%d_etaa%.2f.mat',  indic.spillovers, 0, indic.sep,indic.xgrowth,  etaa));
% sp_not=helper.sp_all';
% helper=load(sprintf('OPT_notarget_spillover%d_taus0_noskill%d_notaul%d_sep%d_extern1_weightext0.01_xgrowth%d_etaa%.2f.mat',indic.spillovers, 0, indic.notaul, indic.sep,indic.xgrowth,  etaa));
% opt_not_wt=helper.opt_all';
% 
% helper=load(sprintf('OPT_notarget_spillover%d_taus0_noskill%d_notaul_sep%d_extern1_weightext0.01_xgrowth%d_etaa%.2f.mat',indic.spillovers, 0, indic.sep, indic.xgrowth, etaa));
% opt_not_nt=helper.opt_all';
% 
% RES_ext = containers.Map({ 'SP' , 'OPT', 'OPT_notaul'},{sp_not, opt_not_wt, opt_not_nt});
% RES_ext = add_vars(RES_ext, list, params, indic, varlist, symms);

%-version for comparison with full model: robustness
% 
% RES_robext = containers.Map({ 'SP_T', 'OPT_T_NoTaus'},{sp_not, opt_not_wt}); % names to compare with baseline with target
% RES_robext = add_vars(RES_robext, list, params, indic, varlist, symms);


% %- results with counterfactual policy
% if plotts.countcomp==1|| plotts.countcomp2==1|| plotts.countcomp3==1
% if indic.xgrowth==1
%     helper=load(sprintf('COMPEqu_SIM_taufopt1_taulopt0_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, params(list.params=='etaa')));
%     count_taufopt=helper.LF_COUNT';
%     helper=load(sprintf('COMPEqu_SIM_taufopt0_taulopt1_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, params(list.params=='etaa')));
%     count_taulopt=helper.LF_COUNT';
% else
%     helper=load(sprintf('COMPEqu_SIM_taufopt1_taulopt0_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
%     count_taufopt=helper.LF_COUNT';
%     helper=load(sprintf('COMPEqu_SIM_taufopt0_taulopt1_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red,  params(list.params=='etaa')));
%     count_taulopt=helper.LF_COUNT';
% end
% RES_count = containers.Map({ 'taufOpt' , 'taulOpt'},{count_taufopt, count_taulopt});
% RES_count = add_vars(RES_count, list, params, indic, varlist, symms);
% end

%% Tables
if plotts.table==1
 %- discount vector
 betaa=params(list.params=='betaa');

% version with T=12
  disc=repmat(betaa, 1,T);
  expp=0:T-1;
 vec_discount= disc.^expp;

 %- Table
 TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1));
 TableSWF_PV.Properties.VariableNames={'Allocation','Integrated', 'Scenario 1: as 0 but no taul', 'Scenario 2: Gov consu no taul', 'Scenario 3: Gov consu with taul', 'Scenario 4: Lump Sum with taul', 'Scenario 5: Lump Sum no taul'};

%- all results
for i =keys(RES)
    ii=string(i);
    allvars= RES(ii);
%     fprintf('%s', ii)
    % SWF calculation 
    TableSWF_PV.Integrated(TableSWF_PV.Allocation==ii)=vec_discount*allvars(find(varlist=='SWF'),:)'+indic.PV*allvars(find(varlist=='PV'),1);
end
%- Other policies
for i =keys(RES_polcomp_notaul0)
     ii=string(i);
     count=0; % to keep track of which container is used
    for ccc=OTHERPOL
        pp=ccc{1};
        count=count+1;
        allvarsnt=pp(ii);
%         hh=vec_discount*allvarsnt(find(varlist_polcomp=='SWF'),:)'+indic.PV*allvarsnt(find(varlist_polcomp=='PV'),1);
%         fprintf('with round %d the value is %.6f', count, hh);
        TableSWF_PV{TableSWF_PV.Allocation==ii,count+2}=vec_discount*allvarsnt(find(varlist_polcomp=='SWF'),:)'+indic.PV*allvarsnt(find(varlist_polcomp=='PV'),1);
    end
end
%%
save(sprintf('Table_SWF_%s_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d.mat',date, indic.sep, indic.noskill, etaa, indic.xgrowth, indic.PV, indic.extern), 'TableSWF_PV');
end

%% table CEV
if plotts.cev==1
    %- calculate CEV for a pair of policy regimes each
    h1= OTHERPOL{3}; % taul can be used
    h2= OTHERPOL{2}; % taul cannot be used
    COMP = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist_polcomp, varlist_polcomp, symms, list, params, T, indic);
end

%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%% plot analy taul
if plotts.analyta==1
gcf=figure('Visible','off');
allvars=RES('OPT_T_NoTaus');
main=plot(time,allvars(varlist=='analyTaul',:), time,allvars(varlist=='taul',:)  );   
set(main,{'LineWidth'}, { 1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, { 'k'; orrange} )   

xticks(txx)
% ylim([-0.2, 3])
xlim([time(1), time(end)])

ax=gca;
ax.FontSize=13;
ytickformat('%.2f')
xticklabels(Year10)
lgd=legend('analytical taul', 'tau_{\iota}', 'Interpreter', 'latex');    
set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
path=sprintf('figures/all_July22/analyT_xgrowth%d_noskill%d.png', indic.xgrowth, indic.noskill);
exportgraphics(gcf,path,'Resolution', 400)
close gcf
end
%% plot emission limit
if plotts.limit==1
gcf=figure(); %('Visible','off');

main=plot(time(3:end), Ems, time(1:end), zeros(size(time)));   
set(main,{'LineWidth'}, { 1.3; 1}, {'LineStyle'},{'-'; '--'}, {'color'}, { orrange; grrey} )   
hold on 
plot(time(1:3), [Ems(1), Ems(1), Ems(1)], 'LineWidth', 0.8, 'LineStyle', '--', 'color', orrange);
hold off

xticks(txx)
ylim([-0.2, 3])
xlim([time(1), time(end)])

ax=gca;
ax.FontSize=13;
ytickformat('%.2f')
xticklabels(Year10)
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
ylabel('GtCO2-eq') 
path=sprintf('figures/all_%s/emission_lim.png', date);
exportgraphics(gcf,path,'Resolution', 400)
close gcf
end
%% comparison laissez-faire
if plotts.lf==1
    fprintf('plotting LF comp graphs')    
    lff=RES('LF');

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
            main=plot(time,allvars(find(varlist==varr),:), time,lff(find(varlist==varr),:), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('social planner', 'laissez-faire', 'Interpreter', 'latex');
               else
                  lgd=legend('ramsey planner', 'laissez-faire', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_1705/%s_LFComp%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,indic.xgrowth,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
        
    end
    end      
end

%% robustness comparison to full model
if plotts.robust==1
    fprintf('plotting robustness graphs')
    for lgdind=0:1
    for ccc=1:2
        if ccc==0
            whh='extern';
            RES_rob=RES_robext;
        elseif ccc==1
            whh='countec';
            RES_rob =RES_countec;
        elseif ccc==2
            whh='xgrowth';
            RES_rob=RES_xgr;
        end
    
    for i =keys(RES_rob)
        ii=string(i);
        allvars= RES(ii);
        allvars_rob=RES_rob(ii);
    %% 
    fprintf('plotting %s',ii );
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            main=plot(time,allvars_rob(find(varlist==varr),:),time,allvars(find(varlist==varr),:), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   

           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if lgdind==1

                lgd=legend('counterfactual', 'baseline', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            path=sprintf('figures/all_1705/SingleROB_%s_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', whh, ii,varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
            
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end 
    end
    end
end
%% counterfactual comparison to full optimal allocation
if plotts.countcomp==1
    allvars=RES('OPT_T_NoTaus');
    for ccc=0:1
        indic.tauf=ccc;
        if indic.tauf==1 % version with tauf optimal
            allvars_count=RES_count('taufOpt');
        else
            allvars_count=RES_count('taulOpt');
        end
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            
                main=plot(time,allvars(find(varlist==varr),:),time,allvars_count(find(varlist==varr),:), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k';  orrange} )   

           xticks(txx)
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
           end
            if lgdind==1
                if indic.tauf==1
                    lgd=legend('full model', '$\tau_{f}$ optimal, $\tau_l=0$', 'Interpreter', 'latex');
                else
                    lgd=legend('full model', '$\tau_{l}$ optimal, $\tau_f=0$', 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
        path=sprintf('figures/all_1705/CompCounterFac_taufopt%d_taulopt%d_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth0_etaa%.2f_lgd%d.png', indic.tauf, (1-indic.tauf), varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
    %% counterfactual comparison to full optimal allocation
    %- with LF
if plotts.countcomp2==1
    allvars=RES('OPT_T_NoTaus');
    allvarsLF=RES('LF');
    for ccc=0:1
        indic.tauf=ccc;
        if indic.tauf==1 % version with tauf optimal
            allvars_count=RES_count('taufOpt');
        else
            allvars_count=RES_count('taulOpt');
        end
    for lgdind=0:1
    for l ="HH"%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            main=plot(time,allvarsLF(find(varlist==varr),:), time,allvars_count(find(varlist==varr),:),time,allvars(find(varlist==varr),:), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k';  orrange; 'b'} )   

           xticks(txx)
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end
            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
            if lgdind==1
                if indic.tauf==1
                    lgd=legend( 'laissez-faire', '$\tau_{f}$ optimal, $\tau_l=0$', 'optimal policy','Interpreter', 'latex');
                else
                    lgd=legend('laissez-faire',  '$\tau_{l}$ optimal, $\tau_f=0$', 'optimal policy','Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
        path=sprintf('figures/all_1705/CompCounterFac_withLF_taufopt%d_taulopt%d_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_lgd%d.png', indic.tauf, (1-indic.tauf), varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
   %% counterfactual comparison to full optimal allocation
    %- with LF no optimal policy
if plotts.countcomp3==1
    allvarsLF=RES('LF');
    for ccc=0:1
        indic.tauf=ccc;
        if indic.tauf==1 % version with tauf optimal
            allvars_count=RES_count('taufOpt');
        else
            allvars_count=RES_count('taulOpt');
        end
    for lgdind=0:1
    for l = "HH"%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            main=plot(time,allvarsLF(find(varlist==varr),:), time,allvars_count(find(varlist==varr),:), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k';  orrange} )   

           xticks(txx)
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end
            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
            if lgdind==1
                if indic.tauf==1
                    lgd=legend( 'laissez-faire', '$\tau_{f}$ optimal, $\tau_l=0$', 'Interpreter', 'latex');
                else
                    lgd=legend('laissez-faire',  '$\tau_{l}$ optimal, $\tau_f=0$', 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
        path=sprintf('figures/all_1705/CompCounterFac_withLF_noopt_taufopt%d_taulopt%d_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_lgd%d.png', indic.tauf, (1-indic.tauf), varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% No exogenous target but externality
%- comparison across efficient with taul no taul
if plotts.extern==1
   fprintf('plotting externality graphs')
    allvarseff= RES_ext('SP');
    allvarsopt= RES_ext('OPT');
    allvarsopt_nt= RES_ext('OPT_notaul');

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            main=plot(time,allvarseff(find(varlist==varr),:),time,allvarsopt(find(varlist==varr),:), time,allvarsopt_nt(find(varlist==varr),:), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   

           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if lgdind==1

            lgd=legend('efficient', 'with income tax', 'no income tax', 'Interpreter', 'latex');
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
        path=sprintf('figures/all_1705/Extern_CompEff_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_weightext%.2f_xgrowth%d_etaa%.2f_lgd%d.png',  varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,weightext, indic.xgrowth, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% All figures single
if plotts.single==1
    fprintf('plotting single graphs')
    if indic.extern==1
        RESS=RES_ext;
    else
        RESS=RES;
    end
    for i ="LF"%keys(RESS)
        ii=string(i);
        allvars= RESS(ii);
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
            if indic.count_techgap==0
                path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_extern%d_etaa%.2f.png',  ii,varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth,indic.extern,  etaa);
            else
                path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_extern%d_countec_etaa%.2f.png',  ii,varr, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, indic.extern, etaa);
            end
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
        if plotts.regime_gov==0
            nt=0;
            helpp= RES;
            varl=varlist;
        else
            nt =3;
            helpp =OTHERPOL{nt};
            varl= varlist_polcomp;
        end
    for i =keys(helpp)
        %- loop
        ii=string(i);
   
        allvars =helpp(ii);
    %% 
    fprintf('plotting %s',ii );
    for l =keys(lissComp) % loop over variable groups
        ll=string(l);
        plotvars=lissComp(ll); % here plotvars is a group of variable names which are to be plotted in the same graph

        gcf=figure('Visible','off');

        if length(plotvars)==2
            main=plot(time,allvars(find(varl==plotvars(1)),:), time,allvars(find(varl==plotvars(2)),:), 'LineWidth', 1.1);    % plot vectors!        
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} ) 
    %    elseif ll=="LabourInp"
    %           main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),time,allvars(find(varlist==plotvars(3)),:), 'LineWidth', 1.1);    % plot vectors!        
    %           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; 'k'} )   
       elseif ll=="Growth"
              main=plot(time(1:end-1),allvars(find(varl==plotvars(1)),1:end-1), time(1:end-1),allvars(find(varl==plotvars(2)),1:end-1),...
              time(1:end-1),allvars(find(varl==plotvars(3)),1:end-1),time(1:end-1),allvars(find(varl==plotvars(4)),1:end-1),'LineWidth', 1.1);    % plot vectors!        
              set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   

        elseif ll=="Science"
              main=plot(time,allvars(find(varl==plotvars(1)),:), time,allvars(find(varl==plotvars(2)),:),...
                  time,allvars(find(varl==plotvars(3)),:), time,allvars(find(varl==plotvars(4)),:),'LineWidth', 1.1);    % plot vectors!        
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
        if indic.count_techgap==0
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, nt, ii,ll, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,indic.extern, indic.PV,  etaa, lgdind);
        else
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_countec_extern%d_PV%d_etaa%.2f_lgd%d.png',date,  nt, ii,ll, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.extern, indic.PV, etaa, lgdind);
        end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
    end
    end
end

 %% comparison POLICY scenarios
if plotts.notaul==1
    fprintf('plotting comparison across policies') 
    for nt =  1:length(OTHERPOL)
        pp = OTHERPOL{nt};
        count = nt;
        
    for i =keys(pp)

     ii=string(i);
     
     if plotts.regime_gov ==0
         bb=0;
        allvars= RES_polcomp_notaul0(ii);
        varl=varlist; 
     else
         bb=3;
         helpp=OTHERPOL{bb};
         allvars= helpp(ii);
         varl=varlist_polcomp;
     end
     allvarsnt=pp(ii); 
     
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
           for v=1:length(plotvars)
           gcf=figure('Visible','off');

               varr=string(plotvars(v));
               main=plot(time,allvars(find(varl==varr),:), time,allvarsnt(find(varlist_polcomp==varr),:), 'LineWidth', 1.1);

               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
               xticks(txx)
               xlim([1, time(end)])
               ax=gca;
               ax.FontSize=13;
               ytickformat('%.2f')
               xticklabels(Year10)

            if lgdind==1
                if count==1 % version with only no income tax, but integretaed transfers 
                    lgd=legend('benchmark policy', 'no income tax', 'Interpreter', 'latex');
                elseif count == 2 %
                     lgd=legend('benchmark policy', ['no redistribution,' newline 'no income tax'], 'Interpreter', 'latex');
                elseif count == 3 %
                     lgd=legend('benchmark policy', ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex');
                elseif count == 4 %
                     lgd=legend('benchmark policy', ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                elseif count == 5 %
                     lgd=legend('benchmark policy', ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
    %         sgtitle('Social Planner Allocation')
            path=sprintf('figures/all_July22/comp_bb%d_notaul%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', bb, count, ii, varr, indic.spillovers,indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            % saveas(gcf,path)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end
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
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'bau', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'bau', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_1705/%s_BAUComp%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,indic.xgrowth,  etaa, lgdind);
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
        path=sprintf('figures/all_July22/%s_TargetComp%s_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end
%% single alternative regimes with income tax
if plotts.regimes==1
    fprintf('plotting  other regimes with and without taul graphs')   
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt = [3] % 3 is version with gov consumes env revenues but income tax is available 
        
        % read in version to version without income tax
        RES_help=OTHERPOL{nt};
    for i =1:length(opt)
        io=opt(i);
        allvars =RES_help(io);

    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                
           main=plot(time,allvars(find(varlist_polcomp==varr),:));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   

           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

        path=sprintf('figures/all_July22/%s_SingleAltPol%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end


%% comparison social planner and optimal policy (with and without labour tax)
% graph incorporates with and without laissez faire allocation 
% ensure to use different variable list for policy scenario alternations
if plotts.compeff==1
    for withlff=0
    lff=RES('LF');
    fprintf('plotting comparison efficient-optimal graphs')   
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt =  1:length(OTHERPOL) % loop over policy scenarios
        RES_help=OTHERPOL{nt};
    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        if plotts.regime_gov ==0
            bb=0;
           allvars= RES(io);
           varl=varlist;
        else
            bb =3;
            helpp = OTHERPOL{bb};
            allvars = helpp(io);
            varl=varlist_polcomp;
        end
        allvarsnotaul =RES_help(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varlist==varr),1:T), time,allvarseff(find(varlist==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varlist_polcomp==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot(time,allvarseff(find(varlist==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varlist_polcomp==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
           end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

%             if varr=="C"
%                 ylim([0.56, 0.72]);
%             elseif varr=="hh"
%                 ylim([0.46, 0.51]);
%             elseif varr=="hl"
%                 ylim([0.31, 0.335]);                
%             end
           if lgdind==1
               if withlff==1

                   if nt ==1
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  'no income tax', 'Interpreter', 'latex');
                   elseif nt ==2
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['no redistribution,' newline 'no income tax'],  'Interpreter', 'latex');
                   elseif nt ==3
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex'); 
                   elseif nt==4
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                   elseif nt==5
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if nt ==1
                        lgd=legend( 'efficient', 'benchmark policy',  'no income tax', 'Interpreter', 'latex');
                   elseif nt ==2
                        lgd=legend( 'efficient', 'benchmark policy',  ['no redistribution,' newline 'no income tax'], 'Interpreter', 'latex');
                   elseif nt ==3
                        lgd=legend( 'efficient', 'benchmark policy',  ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex'); 
                   elseif nt==4
                        lgd=legend( 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                   elseif nt==5
                        lgd=legend( 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');                        
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_July22/%s_CompEff%s_bb%d_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png', varr, io, bb, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end

%% comparison alternative regimes with and without income tax
if plotts.compeff_dd==1
    for withlff=0
    lff=RES('LF');
    fprintf('plotting comparison efficient- other regimes with and without taul graphs')   
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt = [3,4] % 3 is version with gov consumes env revenues but income tax is available 
        
        % read in version to version without income tax
        if nt==3
            RES_dd=OTHERPOL{2}; % comparison is separate policy regime
        elseif nt ==4
            RES_dd=OTHERPOL{5}; % lump sum trans but no taul
        end
        RES_help=OTHERPOL{nt};
    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        allvars= RES_dd(io);
        allvarswithtaul =RES_help(io);
        allvarseff=RES(ie); 

    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ismember(varr, ["tauf", "taul", "GovCon", "Tls"])
                 
               if withlff==1
                   main=plot(time, lff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, {1; 1.2; 1},  {'LineStyle'},{'--'; '--'; ':'}, {'color'}, {grrey; orrange; 'b'} )   
               else
                   main=plot( time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, { 1.2; 1},  {'LineStyle'},{ '--'; ':'}, {'color'}, { orrange; 'b'} )   
               end

           else
                
               if withlff==1
                   main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
               else
                   main=plot(time,allvarseff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
               end
            end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
            if ismember(varr, ["tauf","taul", "GovCon", "Tls"]) % then no efficient graph
            
               if withlff==1
                  if nt ==3
                    lgd=legend('laissez-faire', 'no income tax', 'with income tax', 'Interpreter', 'latex'); 
                  elseif nt ==4
                     lgd=legend('laissez-faire', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if nt ==3
                        lgd=legend( 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   elseif nt ==4
                        lgd=legend('no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
            else    % with efficient graph            
               if withlff==1
                   if nt==3
                      lgd=legend('laissez-faire', 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   elseif nt ==4
                      lgd=legend('laissez-faire', 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   end
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                  if nt ==3
                        lgd=legend( 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  elseif nt ==4
                        lgd=legend('efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
            end
           end
        path=sprintf('figures/all_July22/%s_DDCompEff%s_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d_lff%d.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
    end
    end
end

%% only social planner
if plotts.compeff1==1
%     lff=RES('LF');
    fprintf('plotting comparison efficient')   
    eff= string({'SP_T', 'SP_NOT'});   
    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           
           main=plot( time,allvarseff(find(varlist==varr),:));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
              lgd=legend('efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_July22/%s_CompEff%s_onlyeff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d.png', varr, ie, indic.spillovers, indic.noskill, indic.sep,indic.xgrowth, indic.count_techgap, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
    for withlff=0
         lff=RES('LF');
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        
        if plotts.regime_gov==0
            nt=0;
            allvars= RES(io);
            varl=varlist;
        else
            nt=3;
            helpp=OTHERPOL{nt};
            allvars=helpp(io);
            varl=varlist_polcomp;
        end
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvars(find(varl==varr),:));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),:), time,allvars(find(varl==varr),:));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', 'optimal policy', 'Interpreter', 'latex');
               else
%                    if varr =="tauf"
%                       lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
%                    else
                     lgd=legend( 'efficient', ' optimal policy', 'Interpreter', 'latex');
%                    end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_July22/%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d_lff%d.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep,indic.xgrowth, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end

%% comparison social planner and optimal policy (with and without labour tax)
if plotts.compeff2==1
    %- only efficient and no income tax
    fprintf('plotting comparison efficient-optimal graphs')   
    for withlff=0
         lff=RES('LF');
    
    if indic.zero==0
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    else
        eff= string({'SP_NOT'});
         opt=string({'OPT_NOT_NoTaus'});
    end
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
      if withlff==1
            main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvarsnotaul(find(varlist==varr),:));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),:), time,allvarsnotaul(find(varlist==varr),:));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', ' no income tax', 'Interpreter', 'latex');
               else
                   if varr =="tauf"
                      lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
                   else
                     lgd=legend( 'efficient', ' no income tax', 'Interpreter', 'latex');
                   end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_1705/%s_CompEff%s_noopt_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_countec%d_etaa%.2f_lgd%d_lff%d.png', varr, io, indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red,indic.xgrowth, indic.zero, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end

end

       