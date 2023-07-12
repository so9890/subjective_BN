function []=plotts_extT(list, T, etaa, weightext,indic, params, Ems, plotts, percon, MOM, StatsEms)

% this script plots results

date="13Sept22_Tplus30";
if ~isfile(sprintf('figures/all_%s', date ))
    mkdir(sprintf('figures/all_%s', date));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.6 0.6 0.6];

%- variables
syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl wlf wlg wln ws wsg wsn wsf lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An A S real
syms analyTaul Hagg PV CEVv CEVvPV CEVvDy Tauf dTaulS dTaulHh dTaulHl dTaulAv dTaulAvS AgAf sffsg sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul Hagg PV Tauf dTaulAvS dTaulS dTaulHh dTaulHl dTaulAv AgAf sffsg sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAagg gAg gAf gAn Utilcon Utillab Utilsci];
% already exists: symms.addgov
symms.comp=[ CEVv CEVvDy CEVvPV ]; % for comparison of policy interventions, 

if indic.sep==0
    symms.plotsvarsPri =[pg pf pee pn wh wl ws];  
else
    symms.plotsvarsPri =[pg pf pee pn wh wl wlf wlg wln wsg wsn wsf];  
end

listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);
listt.plotsvarsAdd=string(symms.plotsvarsAdd);
list.comp=string(symms.comp);
list.growthrates =string([gAg gAf gAn gAagg]);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%- variables which to plot in one graph plus legend
lissComp = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {string([hh hl]), string([sff sg sn S]), string([gAf gAg  gAn gAagg]), string([Lf Lg])});
legg = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {["high skill", "low skill"], ["fossil", "green", "non-energy", "total"], ["fossil", "green", "non-energy", "aggregate"], ["fossil", "green"]});

%- new list 
% varlist_polcomp=[list.allvars, list.addgov, string(symms.plotsvarsAdd)];

varlist=[list.allvars, string(symms.plotsvarsAdd)];
%% read in results

%-initialise cells to store matrices: one for each combination of xgr and nsk
OTHERPOL={}; % cell to save containers of different policy versions
OTHERPOL_xgr={}; 
OTHERPOL_nsk={};
OTHERPOL_xgr_nsk={}; 

% baseline results 
for xgr=0:1
for nsk =0:1
    indic.xgrowth=xgr;
    indic.noskill=nsk;
    %- sp solution independent of policy
%     helper=load(sprintf('SP_notarget_1008_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', ...
%         indic.spillovers, indic.noknow_spill, indic.noskill, 1, indic.xgrowth,  indic.PV, indic.sizeequ, etaa));
%     sp_not=helper.sp_all';

%- other results
    for i=[4,5] % loop over policy versions
        
       helper = load(sprintf('BAU_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ0_sep%d_notaul0_countec%d_GovRev%d_etaa%.2f.mat', ...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, indic.sep, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
        LF = helper.COMP';
        
        helper=load(sprintf('OPT_notarget_plus30_0509_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',... 
            indic.spillovers,indic.noknow_spill, indic.noskill,i,  indic.sep, indic.xgrowth,indic.PV,  plotts.sizeequ, plotts.GOV, etaa));
        opt_not_notaus=helper.opt_all_all(1:T,:)';
        helper=load(sprintf('OPT_target_plus30_0509_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', ... 
            indic.spillovers,indic.noknow_spill, indic.noskill, i, indic.sep, indic.xgrowth,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
        opt_t_notaus=helper.opt_all_all(1:T,:)';
        if indic.noskill~=1

                helper=load(sprintf('SP_target_plus30_2609_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
                indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, etaa));
                sp_t=helper.sp_all_all(1:T,:)';
        else
               sp_t=opt_t_notaus;
        end
        %- with externality 
        %- RES without
        if plotts.nsk==1 && (indic.extern==1 || plotts.extern==1)
           error('do not have results with externality and homo skills')
        else
            if indic.noskill==1
                     opt_t_ext=helper.opt_all_all(1:T,:)';
            else
                helper= load(sprintf('OPT_notarget_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern1_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
                    indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, i,indic.sep,weightext, indic.xgrowth,indic.PV,indic.sizeequ, plotts.GOV, etaa));
                opt_t_ext=helper.opt_all(1:12,:)';
            end
        end

        RES = containers.Map({'LF', 'OPT_T_NoTaus', 'OPT_NOT_NoTaus', 'SP_T', 'Ext'},...
                                { LF, opt_t_notaus, opt_not_notaus, sp_t, opt_t_ext});
        %- add additional variables
        if xgr==0 && nsk==0
            OTHERPOL{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
        elseif xgr==0 && nsk==1
            OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
        elseif xgr==1 && nsk==0
            OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
        elseif xgr==1 && nsk==1
            OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
        end
    end
end
end

%- social planner only with target and nsk 0 xgr 0 
    helper=load(sprintf('SP_target_plus30_2609_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
         indic.spillovers,indic.noknow_spill, 0, indic.sep, 0, indic.PV,indic.sizeequ, etaa));
    sp_t=helper.sp_all_all(1:T,:)';
    RES_SP=containers.Map({'SP_T',}, {sp_t});
    RES_SP=add_vars(RES_SP, list, params, indic, list.allvars, symms, MOM);
  
%- no knowledge spillovers 
     helper=load(sprintf('OPT_notarget_plus30_0509_spillover%d_knspil1_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
         indic.spillovers, plotts.nsk,plotts.regime_gov,  indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
     opt_not_notaus=helper.opt_all_all(1:T,:)';
     helper=load(sprintf('OPT_target_plus30_0509_spillover%d_knspil1_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
         indic.spillovers, plotts.nsk, plotts.regime_gov, indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
     opt_t_notaus=helper.opt_all_all(1:T,:)';
     if plotts.nsk==1 && plotts.xgr==0
         sp_t= opt_t_notaus;
     else
         helper=load(sprintf('SP_target_plus30_2609_spillover%d_knspil1_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
            indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,indic.sizeequ, etaa));
         sp_t=helper.sp_all_all(1:T,:)';
     end    
        if plotts.nsk==1 && (indic.extern==1 || plotts.extern==1)
           error('do not have results with externality and homo skills')
        else
            if plotts.nsk==1
               opt_ext=opt_t_notaus;
            else     
                helper= load(sprintf('OPT_notarget_0509_spillover%d_knspil1_taus%d_noskill%d_notaul%d_sep%d_extern1_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
                    indic.spillovers, indic.taus, plotts.nsk, plotts.regime_gov,indic.sep,weightext, plotts.xgr,indic.PV,indic.sizeequ, plotts.GOV, etaa));
                opt_ext=helper.opt_all(1:12,:)';    
            end
        end
     RES_noknspil=containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus', 'SP_T', 'Ext'},...
                                 {opt_t_notaus, opt_not_notaus, sp_t, opt_ext});
                             
     RES_noknspil=add_vars(RES_noknspil, list, params, indic, list.allvars, symms, MOM);
    % no knowledge spillovers in benchmark
        helper=load(sprintf('OPT_notarget_plus30_0509_spillover%d_knspil1_taus0_noskill0_notaul%d_sep%d_extern0_xgrowth0_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            indic.spillovers, plotts.regime_gov,  indic.sep,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
     opt_not_notaus=helper.opt_all_all(1:T,:)';
     helper=load(sprintf('OPT_target_plus30_0509_spillover%d_knspil1_taus0_noskill0_notaul%d_sep%d_xgrowth0_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
         indic.spillovers,  plotts.regime_gov, indic.sep, indic.PV, plotts.sizeequ, plotts.GOV, etaa));
     opt_t_notaus=helper.opt_all_all(1:T,:)';
     RES_bench_noknspil=containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                                 {opt_t_notaus, opt_not_notaus});
     RES_bench_noknspil=add_vars(RES_bench_noknspil, list, params, indic, list.allvars, symms, MOM);
  
%- only optimal tauf
if plotts.xgr==0 
    helper=load(sprintf('COMPEquN_SIM_1110_taufopt1_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
        indic.noknow_spill,  indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_taul=helper.LF_COUNT';
        helper=load(sprintf('COMPEquN_SIM_1110_taufopt0_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
        indic.noknow_spill,  indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_tauf=helper.LF_COUNT';
        RES_count=containers.Map({'CountOnlyTauf', 'CountOnlyTaul'},{count_taul, count_tauf});
     RES_count=add_vars(RES_count, list, params, indic, list.allvars, symms, MOM);
  
end
%% Tables
if plotts.table==1
for xgr=0:1
for nsk= 0:1
    %- choose respective containers cell
    if xgr ==0 && nsk==0
    OTHERPOLL= OTHERPOL;
    elseif xgr ==1 && nsk==0
        OTHERPOLL= OTHERPOL_xgr;
    elseif xgr ==0 && nsk==1
        OTHERPOLL= OTHERPOL_nsk;
    elseif xgr ==1 && nsk==1
        OTHERPOLL= OTHERPOL_xgr_nsk;
    end
  %- discount vector
      betaa=params(list.params=='betaa');

      % version with T=11
      disc=repmat(betaa, 1,T+1);
      expp=0:T;
      vec_discount= disc.^expp;

     %- Table
     TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1));
     TableSWF_PV.Properties.VariableNames={'Allocation','Integrated', 'Scenario 1: as 0 but no taul', 'Scenario 2: Gov consu no taul', 'Scenario 3: Gov consu with taul', 'Scenario 4: Lump Sum with taul', 'Scenario 5: Lump Sum no taul'};

    %- all results
    for i =keys(RES)
         ii=string(i);
         count=0; % to keep track of which container is used
        for ccc=OTHERPOLL
            pp=ccc{1}; % to get from cell to its content
            count=count+1;
            allvars=pp(ii);
            TableSWF_PV{TableSWF_PV.Allocation==ii,count+1}=vec_discount*allvars(find(varlist=='SWF'),:)'+indic.PV*allvars(find(varlist=='PV'),1);
        end
    end
    %%
    save(sprintf('Table_SWF_%s_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d_sizeequ%d_GOV%d.mat',date, indic.sep, nsk, etaa, xgr, indic.PV, indic.extern, plotts.sizeequ, plotts.GOV), 'TableSWF_PV');
end
end
end

%% Pick main policy version for plots
if plotts.xgr ==0 && plotts.nsk==0
    OTHERPOLL= OTHERPOL;
elseif plotts.xgr ==1 && plotts.nsk==0
    OTHERPOLL= OTHERPOL_xgr;
elseif plotts.xgr ==0 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_nsk;
elseif plotts.xgr ==1 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_xgr_nsk;
end

%% table CEV
if plotts.cev==1
    %- calculate CEV for a pair of policy regimes each
    if plotts.regime_gov==0
        h1= OTHERPOLL{1}; % taul can be used
        h2= OTHERPOLL{2}; % taul cannot be used
    elseif plotts.regime_gov==3
        h1= OTHERPOLL{4}; % taul can be used
        h2= OTHERPOLL{3}; % taul cannot be used
    elseif plotts.regime_gov==4
        h1= OTHERPOLL{5}; % taul can be used
        h2= OTHERPOLL{6}; % taul cannot be used
    end
        
    [COMP, COMPTable] = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);   
    save(sprintf('Table_CEV_%s_regime%d_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d.mat',date, plotts.regime_gov,  indic.sep, plotts.nsk, etaa, plotts.xgr, indic.PV, indic.extern), 'COMPTable');
end
%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));



%% Comparison model versions in one graph
% in levels
if plotts.count_modlev==1
    
    fprintf('plott counterfactual model level efficient and optimal policy')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};
        allvarseff=RES("SP_T");
        allvars= RES("OPT_T_NoTaus");
        allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarseff(find(varlist==varr),1:T),time,allvars(find(varlist==varr),1:T) ,time,allvarscount(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'k'; orrange} )   
            if lgdind==1
               lgd=legend('social planner' , 'optimal policy', 'only $\tau_F^*$',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end-1)])

           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTAUF_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr ,plotts.nsk, plotts.xgr, indic.noknow_spill, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% 
if plotts.count_devs==1
    
    fprintf('plott counterfactual deviation from optimal pol => role of taul')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};
        allvars= RES("OPT_T_NoTaus");
        allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        perdif= 100*(allvars-allvarscount)./allvarscount;
        
    for l =["Res"] %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdif(find(varlist==varr),1:T),time,zeros(size(perdif(find(varlist==varr),1:T))), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
            
            
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

           if (varr =="sg"|| varr =="sff" ||varr =="sn" ) && plotts.nsk==1
               ylim ([-0.0201, 0.0201])
               ytickformat('%.3f')
           elseif varr =="GFF"
                ylim ([-0.25, 0.1001])
                ytickformat('%.2f')
           elseif varr =="sgsff"
                ylim ([-1, 0.4001])
                ytickformat('%.1f')
           else
                ytickformat('%.2f')
           end
            ax=gca;
            ax.FontSize=13;
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTAUFPerDif_Opt_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f.png',date, varr ,plotts.nsk, plotts.xgr, indic.noknow_spill, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
end

%% comparison only carbon tax (reg 5) to only tauf from combined
if plotts.count_devs_fromcto==1
    
    fprintf('plott counterfactual deviation from carbon tax only pol => role of adjustment in tauf')

    %- read in variable container of chosen regime
    if plotts.regime_gov==4
        RES=OTHERPOLL{5+1};
    end
        allvars= RES("OPT_T_NoTaus");
        allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        perdif= 100*(allvarscount-allvars)./allvars;
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdif(find(varlist==varr),1:T),time,zeros(size(perdif(find(varlist==varr),1:T))), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
            
            
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

            ax=gca;
            ax.FontSize=13;
            if varr== "sgsff"
                ytickformat('%.1f')
            elseif varr == "Hagg"
                ytickformat('%.2f')
                ylim([-0.12, 0.061])
            else
                 ytickformat('%.2f')
            end
            
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CountTAUF_CTOPer_Opt_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f.png',date, varr ,plotts.nsk, plotts.xgr, indic.noknow_spill, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
end
%% comparison only carbon tax (reg 5) to only tauf from combined
if plotts.count_devs_both==1
    
    fprintf('plott counterfactual deviation from carbon tax only pol => role of adjustment in tauf')

    %- read in variable container of chosen regime
    if plotts.regime_gov==4
        RESntaul=OTHERPOLL{5+1};
        RES=OTHERPOLL{4+1};
    end
        allvarsntaul= RESntaul("OPT_T_NoTaus");
        allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        perdiftauf= 100*(allvarscount-allvarsntaul)./allvarsntaul;
        
         
        allvars= RES("OPT_T_NoTaus");
        %allvarscount=RES_count("CountOnlyTauf"); % version with only tauf
        perdiftaul= 100*(allvars-allvarscount)./allvarsntaul; % relative to carbon-tax-only policy
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdiftauf(find(varlist==varr),1:T), time,perdiftaul(find(varlist==varr),1:T), time,zeros(size(perdiftauf(find(varlist==varr),1:T))), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'; '--'}, {'color'}, {'k';orrange; grrey} )   
            
            
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

            ax=gca;
            ax.FontSize=13;
            if varr== "sgsff"
                ytickformat('%.1f')
            elseif varr == "Hagg"
                ytickformat('%.2f')
                ylim([-0.12, 0.061])
            else
                 ytickformat('%.2f')
            end
            
            if lgdind==1
               lgd=legend('effect $\tau_F$' , 'effect $\tau_\iota$', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            xticklabels(Year10)
            if indic.new
            path=sprintf('figures/all_%s/CountTAUF_Both_Opt_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr ,plotts.nsk, plotts.xgr, indic.noknow_spill, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% Comparison model versions in one graph
% in levels
if plotts.compnsk_xgr1==1
    
    fprintf('plotting comp nsk xgr')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
    RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};  
    
    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
        allvarsnsk=RESnsk(ii);
        allvarsxgr=RESxgr(ii);
        allvarsnknk= RES_bench_noknspil(ii);
        
    fprintf('plotting %s',ii );
    for lgdind=1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsxgr(find(varlist==varr),1:T) ,...
                time,allvarsnsk(find(varlist==varr),1:T), time,allvarsnknk(find(varlist==varr),1:T),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'; '-.'}, {'color'}, {'k'; 'b'; orrange; grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills', 'no knowledge spillovers', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'eastoutside', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CompMod1_%s_%s_regime%d_spillover%d_knspil%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, indic.noknow_spill, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
    
%% Comparison model versions
% in levels
if plotts.compnsk_xgr==1
    
    fprintf('plotting comp nsk xgr')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    for kk="xgr_nsk" %["nsk" "xgr" ]
        if kk =="nsk" % no skill
            RESalt=OTHERPOL_nsk{plotts.regime_gov+1};
        elseif kk =="xgr"
            RESalt=OTHERPOL_xgr{plotts.regime_gov+1};
        elseif kk=="xgr_nsk"
            RESalt=OTHERPOL_xgr_nsk{plotts.regime_gov+1};
        end
  
    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
        allvarsalt=RESalt(ii);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsalt(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
                if kk == "nsk"
                    lgd=legend('benchmark' , 'homogeneous skills',  'Interpreter', 'latex');
                elseif kk=="xgr"
                    lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
                elseif kk=="xgr_nsk"
                    lgd=legend('benchmark' , ['exogenous growth,' newline  'homogeneous skills'],  'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/CompMod_%s_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
    end % loop over model versions
end
%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.compEff_mod_dev1==1
    
    fprintf('plotting comp models efficient 1')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};

  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvarseff=RES(ef);
        allvarsnskeff=RESnsk(ef);
        allvarsxgreff=RESxgr(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarseff(find(varlist==varr),1:T),time,allvarsxgreff(find(varlist==varr),1:T),...
                        time,allvarsnskeff(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
                lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/EFF_CompMod_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% Comparison model versions: in deviation from efficient
% in levels
if plotts.compnsk_xgr_dev==1
    
    fprintf('plotting comp nsk xgr percent from efficient')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
    for kk= "xgr_nsk" %["nsk" "xgr"]
        if kk =="nsk" % no skill
            RESalt=OTHERPOL_nsk{plotts.regime_gov+1};
        elseif kk =="xgr"
            RESalt=OTHERPOL_xgr{plotts.regime_gov+1};
        elseif kk=="xgr_nsk"
            RESalt=OTHERPOL_xgr_nsk{plotts.regime_gov+1};
        end
  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvars= RES(ii);
        allvarseff=RES(ef);
        allvarsalt=RESalt(ii);
        allvarsalteff=RESalt(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,100*(allvars(find(varlist==varr),1:T)-allvarseff(find(varlist==varr),1:T))./allvarseff(find(varlist==varr),1:T),time,100*(allvarsalt(find(varlist==varr),1:T)-allvarsalteff(find(varlist==varr),1:T))./allvarsalteff(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
                if kk == "nsk"
                    lgd=legend('benchmark' , 'homogeneous skills',  'Interpreter', 'latex');
                elseif kk=="xgr"
                    lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
                elseif kk=="xgr_nsk"
                    lgd=legend('benchmark' , ['exogenous growth,' newline  'homogeneous skills'],  'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Per_CompMod_%s_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
    end % loop over model versions
end
%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.compnsk_xgr_dev1==1
    
    fprintf('plotting comp nsk xgr percent from efficient 1')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};

  
    %- loop over economy versions
    for ind = 1:2
        opt=["OPT_T_NoTaus" "OPT_NOT_NoTaus"];% only plotting polcies separately
        eff=["SP_T" "SP_NOT"];% only plotting polcies separately

        ii=string(opt(ind));
        ef=string(eff(ind));
        
        allvars= RES(ii);
        allvarseff=RES(ef);
        allvarsnsk=RESnsk(ii);
        allvarsnskeff=RESnsk(ef);
        allvarsxgr=RESxgr(ii);
        allvarsxgreff=RESxgr(ef);
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,100*(allvars(find(varlist==varr),1:T)-allvarseff(find(varlist==varr),1:T))./allvarseff(find(varlist==varr),1:T),time,100*(allvarsxgr(find(varlist==varr),1:T)-allvarsxgreff(find(varlist==varr),1:T))./allvarsxgreff(find(varlist==varr),1:T),...
                        time,(allvarsnsk(find(varlist==varr),1:T)-allvarsnskeff(find(varlist==varr),1:T))./allvarsnskeff(find(varlist==varr),1:T)*100, 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
            if lgdind==1
                lgd=legend('benchmark' , 'exogenous growth', 'homogeneous skills',  'Interpreter', 'latex');
                
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Per1_CompMod_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% emission limit
if plotts.ems==1
    gcf=figure('Visible','off');
    main= plot([0,time(percon+1:end)],zeros(1,T+1), [0,time(percon+1:end)],[StatsEms.netems_sum1519,Ems(1:T)],[time(percon+1:end)],[Ems(1:T)],  'LineWidth', 1.1);  
    set(main, {'LineStyle'},{'--';'--'; '-'}, {'color'}, { grrey;'k'; 'k'} ) 
    xticks(txx)
    xlim([0, time(end-1)])
    ylim([-0.5, 25.5])
    ytickformat('%.0f')
    ax=gca;
    ax.FontSize=13;
    xticklabels(Year10)
    path=sprintf('figures/all_%s/Emnet.png',date);
    
    exportgraphics(gcf,path,'Resolution', 400)
    close gcf

end
%% ems goals
if plotts.ems_goals==1
    for lgdind=0:1
        for oo=0:1
    gcf=figure('Visible','off');
    if oo==0
        main= plot( time(percon+1:end),StatsEms.netems_sum1519*ones(1,T), ...
             time(percon+1:end),StatsEms.Emslimit_constantEmsRat_Budget(1:T), 'LineWidth', 1.1);  
    set(main, {'LineStyle'},{'--';'-'}, {'color'}, {grrey; 'k'} ) 
    else
         main= plot( time(percon+1:end),StatsEms.netems_sum1519*ones(1,T), ...
             time(percon+1:end),StatsEms.Emslimit_constantEmsRat_Budget(1:T),...
             time(percon+1+2:end),[StatsEms.targetBidenGTCO2_modelperiod*ones(1,4), zeros(1,6)],...
             time(percon+1:end),Ems, 'LineWidth', 1.1);  
        set(main, {'LineStyle'},{'--';'-'; '--'; '--'}, {'color'}, {grrey; 'k'; 'b'; orrange} ) 
    end
    xticks(txx)
    xlim([1, time(end-1)])
    ylim([-0.5, 25.5])
    ytickformat('%.1f')
    ax=gca;
    ax.FontSize=13;
    xticklabels(Year10)
    if lgdind==1
        lgd=legend('2015-2019 US net emissions', 'IPCC carbon budget, equal \% reduction', 'Biden target','IPCC carbon budget, equal-per-capita','Interpreter', 'latex');
        set(lgd, 'Interpreter', 'latex', 'Location', 'bestoutside', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
    end
    path=sprintf('figures/all_%s/Emnet_goals_o%d_lgd%d.png',date,oo, lgdind);
    
    exportgraphics(gcf,path,'Resolution', 400)
    close gcf
        end
    end
end

%% All figures single
if plotts.single_pol==1
    
    fprintf('plotting single graphs')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};

    %- loop over economy versions
    for i = ["OPT_T_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
    fprintf('plotting %s',ii );
    for l ="Add" %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
             if ll=="HH" && varr=="Emnet"
                 main=plot(time,allvars(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T), 'LineWidth', 1.1);  
                 set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                 lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

             else
                main=plot(time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

             end
           xticks(txx)
            xlim([1, time(end-1)])
            xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
        
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            else
                ytickformat('%.2f')
            end
            xticklabels(Year10)
            if indic.count_techgap==0
                path=sprintf('figures/all_%s/Single_periods%d_%s_%s_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date, T, ii,varr , plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, plotts.sizeequ, plotts.GOV,   etaa);
            else
                path=sprintf('figures/all_%s/Single_periods%d_%s_%s_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_extern%d_countec_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date, T, ii, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV,plotts.sizeequ, plotts.GOV, etaa);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% figures single overlayed
if plotts.singov==1
    fprintf('plotting single overlayed graphs')
    RES=OTHERPOLL{plotts.regime_gov+1};
    varl=varlist;
    for lgdind=0:1
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]
        ii=string(i);
        allvars =RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lissComp) % loop over variable groups
        ll=string(l);
        plotvars=lissComp(ll); % here plotvars is a group of variable names which are to be plotted in the same graph

        gcf=figure('Visible','off');

        if length(plotvars)==2
            main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T), 'LineWidth', 1.1);    % plot vectors!        
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} ) 
    %    elseif ll=="LabourInp"
    %           main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),time,allvars(find(varlist==plotvars(3)),:), 'LineWidth', 1.1);    % plot vectors!        
    %           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; 'k'} )   
       elseif ll=="Growth"
              main=plot(time(1:T),allvars(find(varl==plotvars(1)),1:T), time(1:T),allvars(find(varl==plotvars(2)),1:T),...
              time(1:T),allvars(find(varl==plotvars(3)),1:T),time(1:T),allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
              set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   

        elseif ll=="Science"
              main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T),...
                  time,allvars(find(varl==plotvars(3)),1:T), time,allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
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
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, plotts.regime_gov, ii,ll, indic.spillovers,  indic.noknow_spill,  plotts.nsk, indic.sep, plotts.xgr,indic.extern, indic.PV,  etaa, lgdind);
        else
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_countec_extern%d_PV%d_etaa%.2f_lgd%d.png',date,  plotts.regime_gov, ii,ll, indic.spillovers, indic.noknow_spill,  plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV, etaa, lgdind);
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
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    for nt =  bb
        pp = OTHERPOLL{nt};
        count = nt-1; % subtract 1 to get indicator of notaul
        
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]

     ii=string(i);
     varl=varlist; 
     %- benchmark policy
     helpp=OTHERPOLL{plotts.regime_gov+1};
     allvars= helpp(ii);
     %- comparison
     allvarsnt=pp(ii); 
     
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst)
        ll=string(l);
      
        plotvars=lisst(ll);
      
        for v=1:length(plotvars)
           gcf=figure('Visible','off');

               varr=string(plotvars(v));
               main=plot(time,allvars(find(varl==varr),1:T), time,allvarsnt(find(varlist==varr),1:T), 'LineWidth', 1.1);

               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'b'} )   
               xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
               ax=gca;
               ax.FontSize=13;
               ytickformat('%.2f')
               xticklabels(Year10)

            if lgdind==1
                if count==0 % integrated policy
                    lgd=legend('benchmark', 'consolidated budget, with income tax', 'Interpreter', 'latex');
                elseif count == 1 % integrated without income tax
                     lgd=legend('benchmark', 'consolidated budget, no income tax', 'Interpreter', 'latex');
                elseif count == 2 % notaul =2 gov=tauf*F*pf, without income tax
                     lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
                elseif count == 3 % notaul=3
                     lgd=legend('benchmark', 'no redistribution, with income tax', 'Interpreter', 'latex');
                elseif count == 4 % Tls and taul
                     lgd=legend('benchmark', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                elseif count == 5 % Tls no taul
                     lgd=legend('benchmark', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            path=sprintf('figures/all_%s/comp_benchregime%d_notaul%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, plotts.regime_gov, count, ii, varr, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end
    end
end

%% comparison with and without target    
%- string to loop over 
if plotts.comptarg==1
    fprintf('plotting comparison target graphs')
    RES=OTHERPOLL{plotts.regime_gov+1};
    ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,3] % to filter benchmark: policy with target
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));

         allvars= RES(t);
         allvarsnot=RES(nt); 


    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time,allvars(find(varlist==varr),1:T), time,allvarsnot(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('wih emission limit', 'no emission limit', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_TargetComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, varr, ii, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end

%% comparison social planner and optimal policy benchmark and comparison policy
% graph incorporates with and without laissez faire allocation 
if plotts.compeff==1
    fprintf('plotting comparison efficient-optimal and non benchmark graphs')   

    %- read in container of results
    RES=OTHERPOLL{plotts.regime_gov+1};
%     bb=1:length(OTHERPOLL);
%     bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    
    for withlff=0
        lff=RES('LF');
        eff= string({'SP_T'});
        opt=string({'OPT_T_NoTaus'});

    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        varl=varlist;
        allvarseff=RES(ie); 

     %- comparison         
     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{2+1}; % version without taul
     elseif plotts.regime_gov == 0
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{1+1}; % version without taul 
     elseif plotts.regime_gov == 4         
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{5+1}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        allvarsnt =RESnt('OPT_T_NoTaus');
  
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varl==varr),1:T),  time,allvars(find(varl==varr),1:T), time,allvarsnt(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot( time,allvars(find(varl==varr),1:T), time,allvarsnt(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; orrange} )   
           end
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1

                  lgd=legend('laissez-faire', 'with income tax', 'without income tax',  'social planner', 'Interpreter', 'latex');
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                  lgd=legend( 'with income tax', 'without income tax',  'social planner', 'Interpreter', 'latex');
                                 
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_benchregime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png',date, varr, io, plotts.regime_gov, indic.spillovers,indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
    end
end


%% only social planner
if plotts.compeff1==1
    fprintf('plotting efficient')   
        RES=OTHERPOLL{plotts.regime_gov+1}; % same for all regimes
        allvarseff=RES('SP_T'); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if varr ~="Emnet"
               main=plot( time,allvarseff(find(varlist==varr),1:T));            
               set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'--'}, {'color'}, {'k'} )  
           else
               main=plot( time,allvarseff(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T));            
               set(main,{'LineWidth'}, {1.2; 1.2},  {'LineStyle'},{'--'; '--'}, {'color'}, {'k'; orrange} )  
           end
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
              lgd=legend('social planner', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff_Target_onlyeff_reg%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d.png', date, varr,plotts.regime_gov, indic.spillovers,indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
end

%% comparison social planner and non-benchmark policy
if plotts.compeff2==1
    %- only efficient and no income tax
    fprintf('plotting comparison efficient-non benchmark optimal graphs')   
    
    for withlff=0
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
   
    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        io=opt(i);
    for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        lff=RES('LF');
        allvarsnotaul =RES(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),1:T), time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
          
           if lgdind==1
               if withlff==1
                   if count ==0
                        lgd=legend('laissez-faire', 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire', 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               else
                  if count ==0
                        lgd=legend( 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend( 'efficient', 'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend( 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'efficient',  'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'efficient',  'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_noopt_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, count, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end
%% LF versus optimal policy in levels
if plotts.comp_LFOPT==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels with LF and no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(revall(find(varlist==varr),1:T)), time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-.';'-'; '--'}, {'color'}, {grrey; 'k'; 'b'} )   
           xticks(txx)
            xlim([1, time(end-1)])
          
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('laissez-faire', 'with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_LF_OPT_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end

%% optimal with and without taul in levels 
if plotts.comp_OPT==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt and no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{2+1}; % version without taul
     elseif plotts.regime_gov == 0
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{1+1}; % version without taul 
     elseif plotts.regime_gov == 4         
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{5+1}; % version without taul
     end
     for i ={'OPT_T_NoTaus'} %'OPT_T_NoTaus'
         ii=string(i);
     if plotts.extern==0
        allvars= RES(ii);
        allvarsnt =RESnt(ii);
     elseif plotts.extern==1
        allvars= RES('Ext');
        allvarsnt =RESnt('Ext');
     end
    %% 
    for l = "Add" %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

            ax=gca;
            ax.FontSize=13;
            if varr=="dTaulAv"
                ytickformat('%.1f')
            elseif varr=="Tauf"
                ytickformat('%.0f')
            else
                ytickformat('%.2f')
            end
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

           if indic.extern==0
                path=sprintf('figures/all_%s/%s_%s_COMPtaul_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
           else
                path=sprintf('figures/all_%s/%s_OPT_COMPtaul_regime%d_extern%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov,plotts.extern, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
           end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
     end
     end
end
%% optimal with and without taul in levels 
if plotts.comp_OPTPer==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage opt and no taul no eff') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; % version without taul
     elseif plotts.regime_gov == 0
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{1+1}; % version without taul 
     elseif plotts.regime_gov == 4         
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{5+1}; % version without taul
     end
     
     for i=  {'OPT_T_NoTaus'} % 'OPT_T_NoTaus'
         ii=string(i);
     if plotts.extern==0
        allvars= RES(ii);
        allvarsnt =RESnt(ii);
     elseif plotts.extern==1
        allvars= RES('Ext');
        allvarsnt =RESnt('Ext');
     end
        Perdif=100*(allvars-allvarsnt)./allvarsnt;
    %% 
    for l = "Add" %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(Perdif(find(varlist==varr),1:T)),time,zeros(size(time)));            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey},{'LineWidth'}, {1.1; 1} )   
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
           
            ax=gca;
            ax.FontSize=13;
            if varr=="SWF" || varr== "sn"
               ytickformat('%.0f')
            elseif varr=="sff" || varr=="sg" ||  varr == "sgsff" || varr =="GFF" || varr =="sffsg"
               ytickformat('%.1f')
            else
               ytickformat('%.2f')
            end
            xticklabels(Year10)
%            if lgdind==1
%               lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
%               set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%            end
        if plotts.extern==0
            path=sprintf('figures/all_%s/%s_%s_COMPtaulPer_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png',date, varr, ii, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa);
        elseif plotts.extern ==1
             path=sprintf('figures/all_%s/%s_COMPtaulPer_regime%d_extern_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa);
        end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
     end
     end
end



%% counterfactual optimal policy in NKS model in benchmark
if plotts.comp_Bench_CountNK==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting opt pol without kn spills in benchmark model') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        allvars= RES('OPT_T_NoTaus');
        allvarsCOUNT =RES_count_NKinBen('all');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsCOUNT(find(varlist==varr),1:T)),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'--'; '-'}, {'color'}, {grrey; 'k'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('optimal policy', 'counterfactual policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_KNCOUNT_FullMod_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end
%% comparison kn in levels
if plotts.comp_OPT_NK==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt benchmark and no knowledge spillovers') 

        RES=OTHERPOLL{plotts.regime_gov+1};
    for i={'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
            ii=string(i);
        allvars= RES(ii);
        allvarsNK =RES_noknspil(ii);
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if ismember(varr, ["dTaulHh" "dTaulHl" "dTaulAv" "taul"]) 
                main=plot(time,(allvars(find(varlist==varr),1:T)), time,(allvarsNK(find(varlist==varr),1:T)), time, zeros(size(time)), 'LineWidth', 1.1);            
                set(main, {'LineStyle'},{'-'; '--'; '--'}, {'color'}, {'k'; 'k';grrey} )   
           else
               main=plot(time,(allvars(find(varlist==varr),1:T)), time,(allvarsNK(find(varlist==varr),1:T)),'LineWidth', 1.1);            
               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
           end
           
           xticks(txx)
           xlim([1, time(end-1)])
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with knowledge spillovers', 'no knowledge spillovers',  'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_%s_KN_FullMod_sizeequ%d_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_GOV%d_etaa%.2f_lgd%d.png',date, varr, ii, plotts.sizeequ, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,plotts.GOV,  etaa, lgdind);
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
    RES=OTHERPOLL{plotts.regime_gov+1};

    for withlff=1
        lff=RES('LF');
    
        eff= string({'SP_T'});
        opt=string({'OPT_T_NoTaus'});

    for i =[1]

        ie=eff(i);
        io=opt(i);
        
        allvars= RES(io);
        varl=varlist;
        allvarseff=RES(ie); 

    for l =["Add"]%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varl==varr),1:T),  time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));  
            if indic.slides ==0
               set(main, {'LineWidth'}, {0.6; 1.2; 1.2}, {'LineStyle'},{'--';'--'; '-'}, {'color'}, {grrey; 'k'; 'k'} )   
            else
               set(main, {'LineWidth'}, {0.6; 1.2; 1.2}, {'LineStyle'},{':';'--'; '-'}, {'color'}, {'k'; 'k'; orrange} )   
            end
      else
            main=plot( time,allvars(find(varl==varr),1:T),time,allvarseff(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
      end
      
      xticks(txx)
      xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
      xlim([1, time(end-1)])
      
            ax=gca;
            ax.FontSize=13;
            if (varr == "sff" && plotts.xgr==0 && plotts.nsk ==0) || varr == "sg" 
                ytickformat('%.2f')
                ylim([-0.1, 0.2])
           
            elseif varr == "sgsff" 
                ytickformat('%.0f')
                ylim([-1e5, 5e5])
             elseif varr == "sffsg" && indic.noknow_spill==0
                 ytickformat('%.0f')
                 ylim([-1, 8.01])
                 
             elseif varr == "GFF" && indic.noknow_spill==0 
                 ytickformat('%.0f')
                 ylim([-5, 30.01])
            
            elseif varr =="hh" || varr =="hl" || varr =="Hagg"
                ytickformat('%.3f')
%             elseif varr == "C"  && plotts.xgr==0 && plotts.nsk ==0 && indic.noknow_spill==0
%                 ytickformat('%.1f')
%                  ylim([0.309999, 0.500001])
%                 
            elseif varr == "C"  
                ytickformat('%.1f')
                 ylim([0.59999, 1.4001])
            elseif varr == "sn"  
                ytickformat('%.1f')                 
            else
               ytickformat('%.2f')
            end

            xticklabels(Year10)

           if lgdind==1
               if withlff==1
                   if indic.slides ==0
                        lgd=legend('laissez-faire',  'social planner', 'optimal policy', 'Interpreter', 'latex');
                   else
                        lgd=legend('laissez-faire', 'first-best',  'optimal policy', 'Interpreter', 'latex');
                   end
               else
%                    if varr =="tauf"
%                       lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
%                    else
                     lgd=legend( 'optimal policy', 'social planner', 'Interpreter', 'latex');
%                    end
               end
               if varr~="sff"
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
               else
                   if indic.noknow_spill==0
                    set(lgd, 'Interpreter', 'latex', 'Location', 'west', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                   else
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                   end
               end
           end
           if indic.slides==0
              path=sprintf('figures/all_%s/%s_CompEff%s_regime%d_opteff_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
           elseif indic.slides ==1
              path=sprintf('figures/all_%s/%s_slides_CompEff%s_regime%d_opteff_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
           end
        exportgraphics(gcf,path,'Resolution', 400)
        
        close gcf
        end
        end
      end
    end
    end
end

%% comparison social planner and optimal policy benchmark
if plotts.compeff4==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient and lf graphs')   
    RES=OTHERPOLL{plotts.regime_gov+1};

        lff=RES('LF');
        eff= string({'SP_T'});

    for i =[1]

        ie=eff(i);
        
        varl=varlist;
        allvarseff=RES(ie); 

    for l =["Add"]%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));  
            set(main, {'LineWidth'}, {0.6; 1.2}, {'LineStyle'},{':'; '--'}, {'color'}, {'k'; 'k'} )   
         
      
      xticks(txx)
      xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
      xlim([1, time(end-1)])
      
            ax=gca;
            ax.FontSize=13;
            if (varr == "sff" && plotts.xgr==0 && plotts.nsk ==0) || varr == "sg" 
                ytickformat('%.2f')
                ylim([-0.1, 0.2])
           
            elseif varr == "sgsff" 
                ytickformat('%.0f')
                ylim([-1e5, 5e5])
             elseif varr == "sffsg" && indic.noknow_spill==0
                 ytickformat('%.0f')
                 ylim([-1, 8.01])
                 
             elseif varr == "GFF" && indic.noknow_spill==0 
                 ytickformat('%.0f')
                 ylim([-5, 30.01])
            
            elseif varr =="hh" || varr =="hl" || varr =="Hagg"
                ytickformat('%.3f')
%             elseif varr == "C"  && plotts.xgr==0 && plotts.nsk ==0 && indic.noknow_spill==0
%                 ytickformat('%.1f')
%                  ylim([0.309999, 0.500001])
%                 
            elseif varr == "C"  
                ytickformat('%.1f')
                 ylim([0.59999, 1.4001])
            elseif varr == "sn"  
                ytickformat('%.1f')                 
            else
               ytickformat('%.2f')
            end

            xticklabels(Year10)

           if lgdind==1
              lgd=legend('laissez-faire', 'first-best', 'Interpreter', 'latex');
           
               
               if varr~="sff"
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
               else
                   if indic.noknow_spill==0
                    set(lgd, 'Interpreter', 'latex', 'Location', 'west', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                   else
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                   end
               end
           end
        path=sprintf('figures/all_%s/%s_slides_CompEffLF_regime%d_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d.png', date, varr,  plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        
        close gcf
        end
      end
    end
    end
end

%% comparison to BAU
if plotts.bau==1
    fprintf('plotting bau graphs') 
for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
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
            main=plot(time,allvars(find(varlist==varr),1:T), time,bau(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'status quo', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'status quo', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_BAUComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end      
end
end

%% comparison to LF
if plotts.lf==1
    fprintf('plotting LF graphs') 
for nt =plotts.regime_gov+1 %  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        bau=RES('LF');

    for i ={'SP_T', 'OPT_T_NoTaus'} % SP_T, 'OPT_NOT_NoTaus'
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
            main=plot(time,bau(find(varlist==varr),1:T), time,allvars(find(varlist==varr),1:T),  'LineWidth', 1.1);            
           set(main, {'LineStyle'},{':'; '-'}, {'color'}, {'k';'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])
          
            ax=gca;
            ax.FontSize=13;
           
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('laissez-faire', 'social planner',  'Interpreter', 'latex');
               else
                  lgd=legend('laissez-faire', 'optimal policy',  'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_LFComp%s_regime%d_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, count,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end      
end
end
%%
if plotts.per_BAUt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end

%%
if plotts.per_LFt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end
%%
 if plotts.per_effopt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to eff vs opt first period') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('OPT_NOT_NoTaus');
        revalleff =RES('SP_NOT');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            revopt = revall(find(varlist==varr), 1);
             reveff = revalleff(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-reveff)./reveff,time,(allvars(find(varlist==varr),1:T)-revopt)./revopt, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageEffOptFirstPeriod_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

%%
 if plotts.per_optd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to opt dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        revall =RES('OPT_NOT_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend( 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageOptDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_baud==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('BAU');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end   
 
 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');

    %% 
    for l = "Add"% keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if varr~="GFF"
                main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,100*(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           else
               main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)/10,time,1/10*(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           end
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
          
            ax=gca;
            ax.FontSize=13;
            if varr=="hl" || varr== "hh" || varr =="Hagg"
                ytickformat('%.1f')
            elseif varr=="GFF"
                ytickformat('%.0f')
            elseif (varr=="C" && indic.noknow_spill==0)  
                ytickformat('%.0f')
                ylim([-0.3, 0.31]*100)

            elseif (varr=="C" && indic.noknow_spill==1)  
                ytickformat('%.0f')
                ylim([-0.3, 0.31]*100)
            elseif  varr =="S"
                ytickformat('%.0f')

            else
               ytickformat('%.0f')
            end

            xticklabels(Year10)
           if lgdind==1
              lgd=legend('optimal policy', 'social planner',  'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDyn_Target_regime%d_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
     end
    
end   
 
%%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd_nt==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to LF dynamic plus no taul') 

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{2+1}; 
     elseif plotts.regime_gov==0
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{1+1};  
     elseif plotts.regime_gov==4
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{5+1}; 
     end
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
           main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,(allvarsnt(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100, time,100*(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; grrey ; orrange} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'social planner', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDynNT_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
end
%% relative change lF dynamic no efficient
if plotts.per_LFd_ne_nt==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to LF dynamic plus no taul no eff') 
     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{3}; % version without taul
     elseif plotts.regime_gov == 0
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{1+1}; % version without taul 
     elseif plotts.regime_gov == 4         
        RES=OTHERPOLL{plotts.regime_gov+1};
        RESnt =OTHERPOLL{5+1}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        revall =RES('LF');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,(allvarsnt(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDynNT_noeff_Target_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
end
end     