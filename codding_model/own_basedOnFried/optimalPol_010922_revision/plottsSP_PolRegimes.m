function []=plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon, MOM)

% this script plots results

date="13Sept22";
if ~isfile(sprintf('figures/all_%s', date ))
    mkdir(sprintf('figures/all_%s', date));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.6 0.6 0.6];

%- variables
syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl wlf wlg wln ws wsg wsn wsf lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An A S real
syms analyTaul Hagg PV CEVv CEVvPV CEVvDy Tauf dTaulHh dTaulHl dTaulAv AgAf sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg  S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul Hagg PV Tauf dTaulHh dTaulHl dTaulAv AgAf sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAagg gAg gAf gAn Utilcon Utillab Utilsci];
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
TAUFS={};
TAUFS_xgr={};
TAUFS_nsk={};
TAUFS_xgr_nsk={};
% baseline results 
for lablab =0:1 % with equal and non-equal labor supply % MAKE sure this parameter is not used anymore further down
    indic.labshareequ=lablab; % relevant for additional variables
for xgr=0
for nsk =0
    indic.xgrowth=xgr;
    indic.noskill=nsk;
%- other results
    for i=[5] % loop over policy versions
        helper=load(sprintf('COMP_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill,  indic.noskill, indic.xgrowth,lablab, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
        all_TaulC=helper.COMP';
        helper=load(sprintf('COMP_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth,lablab, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
        all_Taul0=helper.COMP';
        helper=load(sprintf('TAUF_1409_taulZero0_knspil%d_limit%d_EmsBase%d_xgr%d_nsk%d_labequ%d_countec%d_GovRev%d_sep%d',...
            indic.noknow_spill, indic.limit_LF,indic.emsbase, indic.xgrowth, indic.noskill,lablab, indic.count_techgap, indic.GOV, indic.sep), 'TAUF');
        tauff_TaulC=helper.TAUF; % tauff is independent of notaul!
        helper=load(sprintf('TAUF_1409_taulZero1_knspil%d_limit%d_EmsBase%d_xgr%d_nsk%d_labequ%d_countec%d_GovRev%d_sep%d',...
           indic.noknow_spill, indic.limit_LF,indic.emsbase, indic.xgrowth, indic.noskill,lablab, indic.count_techgap, indic.GOV, indic.sep), 'TAUF');
        tauff_Taul0=helper.TAUF; % tauff is independent of notaul!
        helper = load(sprintf('BAU_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, lablab, indic.sep, i, indic.count_techgap, indic.GOV,  params(list.params=='etaa')));
        BAU = helper.COMP';
        helper = load(sprintf('BAU_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_countec%d_GovRev%d_etaa%.2f.mat', ...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, lablab, indic.sep, i, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
        LF = helper.COMP';
        RES = containers.Map({'TaulCalib', 'Taul0', 'BAU', 'LF'},{all_TaulC, all_Taul0, BAU, LF});
        %- add additional variables
        if lablab ==0
            if xgr==0 && nsk==0
                OTHERPOL{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
                TAUFS=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==0 && nsk==1
                OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
                TAUFS_nsk=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==1 && nsk==0
                OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
                TAUFS_xgr=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==1 && nsk==1
                OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
                TAUFS_xgr_nsk=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            end
        else
            if nsk==0 && xgr ==0
                OTHERPOL_EQULab{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
            elseif nsk==1 && xgr ==0
                OTHERPOL_EQULab_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
             elseif xgr==1 && nsk==0
                 OTHERPO_EQULabL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
            elseif xgr==1 && nsk==1
                OTHERPOL_EQULab_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms, MOM);
            end

        end
    end

end
end
end


%- for comparison model without knowledge spillovers
%     helper=load(sprintf('COMP_1409_taulZero0_spillovers%d_knspil1_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers, plotts.nsk, plotts.xgr,0, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%         all_TaulC=helper.COMP';
%         helper=load(sprintf('COMP_1409_taulZero1_spillovers%d_knspil1_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers,  plotts.nsk, plotts.xgr,0, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%         all_Taul0=helper.COMP';
%      RES_NK = containers.Map({'TaulCalib', 'Taul0'},{all_TaulC, all_Taul0});
%      RES_NK=add_vars(RES_NK, list, params, indic, list.allvars, symms, MOM);

 %- earmarking
%    helper=load(sprintf('COMP_1409_taulZero0_spillovers%d_knspil0_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers, plotts.nsk, plotts.xgr,0, indic.sep, 7,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%         all_TaulC=helper.COMP';
%          helper=load(sprintf('COMP_1409_taulZero1_spillovers%d_knspil0_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%              indic.spillovers,  plotts.nsk, plotts.xgr,0, indic.sep, 7,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%          all_Taul0=helper.COMP';
%      RES_GS = containers.Map({'TaulCalib', 'Taul0'},{all_TaulC, all_Taul0});
%      RES_GS=add_vars(RES_GS, list, params, indic, list.allvars, symms, MOM);
 %- partial equilibrium sep=2
%         helper=load(sprintf('COMP_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep2_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers, indic.noknow_spill, plotts.nsk, plotts.xgr,0, indic.sep, 0,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%         all_TaulC=helper.COMP';
%         helper=load(sprintf('COMP_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep2_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers, indic.noknow_spill, plotts.nsk, plotts.xgr,0,  0,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
%         all_Taul0=helper.COMP';
%         helper = load(sprintf('BAU_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep2_notaul%d_countec%d_GovRev%d_etaa%.2f.mat',...
%             indic.spillovers, indic.noknow_spill, plotts.nsk, plotts.xgr, 0, 0, indic.count_techgap, indic.GOV,  params(list.params=='etaa')));
%        
%         BAU = helper.COMP';
%         helper = load(sprintf('BAU_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep2_notaul%d_countec%d_GovRev%d_etaa%.2f.mat', ...
%             indic.spillovers, indic.noknow_spill, plotts.nsk, plotts.xgr,0, 0, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
%         LF = helper.COMP';
%         RES_Par = containers.Map({'TaulCalib', 'Taul0', 'BAU', 'LF'},{all_TaulC, all_Taul0, BAU, LF});
%         RES_Par=add_vars(RES_Par, list, params, indic, list.allvars, symms, MOM);
% 
%% Pick main policy version for plots
if plotts.xgr ==0 && plotts.nsk==0
    OTHERPOLL= OTHERPOL;
    OTHERPOLLEL = OTHERPOL_EQULab;
    TAUFF =TAUFS;
elseif plotts.xgr ==1 && plotts.nsk==0
    OTHERPOLL= OTHERPOL_xgr;
    OTHERPOLLEL =  OTHERPO_EQULabL_xgr;
    TAUFF =TAUFS_xgr;
elseif plotts.xgr ==0 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_nsk;
    OTHERPOLLEL = OTHERPOL_EQULab_nsk;
    TAUFF =TAUFS_nsk;
elseif plotts.xgr ==1 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_xgr_nsk;
    OTHERPOLLEL = OTHERPOL_EQULab_xgr_nsk;    
    TAUFF =TAUFS_xgr_nsk;
end

 %% table CEV
% if plotts.cev==1
%     %- calculate CEV for a pair of policy regimes each
%     if plotts.regime_gov==0
%         h1= OTHERPOLL{1}; % taul can be used
%         h2= OTHERPOLL{2}; % taul cannot be used
%     elseif plotts.regime_gov==3
%         h1= OTHERPOLL{4}; % taul can be used
%         h2= OTHERPOLL{3}; % taul cannot be used
%     elseif plotts.regime_gov==4
%         h1= OTHERPOLL{5}; % taul can be used
%         h2= OTHERPOLL{6}; % taul cannot be used
%     end
%         
%     [COMP, COMPTable] = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);   
%     save(sprintf('Table_CEV_%s_regime%d_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d.mat',date, plotts.regime_gov,  indic.sep, plotts.nsk, etaa, plotts.xgr, indic.PV, indic.extern), 'COMPTable');
% end
%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%% plotting necessary tauf 
if plotts.tauf_comp==1
    
    fprintf('plott tauf needed for emission limit')
    for k=keys(TAUFF)
        kk=string(k);
        TaufP=TAUFF(kk);  
    for lgdind=0:1
            gcf=figure('Visible','off');
            main=plot(time,TaufP(1:T, 4+1), time,TaufP(1:T, 0+1), time,TaufP(1:T, 7+1), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'-.'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
            if lgdind==1
               lgd=legend('lump-sum transfers', 'consolidated budget' , 'green subsidies',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/TAUFCO2_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date,kk, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap,indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
    end
end
%% TAUF comparison with and without taul
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.tauf_compTaul==1
    
    fprintf('plott tauf needed for emission limit Comp levels with ')
    TaufTC=TAUFF('TaulCalib');  
    TaufT0=TAUFF('Taul0');  
    TaufDif=(TaufTC-TaufT0)./TaufT0.*100;
    for lgdind=0:1
            gcf=figure('Visible','off');
            main=plot(time,TaufDif(1:T, 4+1), time,TaufDif(1:T, 0+1), time,TaufDif(1:T, 7+1), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'-.'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'; grrey} )   
            if lgdind==1
               lgd=legend('lump-sum transfers', 'consolidated budget' , 'green subsidies',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/TAUFCO2_PerDifTAUL_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
end

%% plotting necessary tauf 
if plotts.tauf_comp_Byregime==1
    
    fprintf('plott tauf needed for emission limit with and without taul')
    TaufTC=TAUFF('TaulCalib');  
    TaufT0=TAUFF('Taul0');  

    for lgdind=0:1
            gcf=figure('Visible','off');
            main=plot(time,TaufTC(1:T, plotts.regime+1),time,TaufT0(1:T, plotts.regime+1), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
            if lgdind==1
               lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$' , 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/TAUFCO2_LevDifTAUL_regime%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, plotts.regime, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
end
%% TAUF comparison with and without taul
% compare effect of only taul in benchmark and in exogenous growth model

if plotts.tauf_compTaul_BYregime==1
    
    fprintf('plott tauf needed for emission limit with and without taul ')
    TaufTC=TAUFF('TaulCalib');  
    TaufT0=TAUFF('Taul0');  
    TaufDif=(TaufTC-TaufT0)./TaufT0.*100;
    for lgdind=0
            gcf=figure('Visible','off');
            main=plot(time,TaufDif(1:T, plotts.regime+1), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
            if lgdind==1
               lgd=legend( 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/TAUFCO2_PerDifTAUL_regime%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, plotts.regime, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Effect of tauf allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% effect of tauf percent
if plotts.perDif_notauf==1
       fprintf('plott effect tauf percent')
 for lablab =0
     if lablab ==0
        Cons=OTHERPOLL{plotts.regime+1};
     else
         Cons=OTHERPOLLEL{plotts.regime+1};
     end
    for k={'TaulCalib', 'Taul0'}% keys are with and without taul
        kk=string(k);
        allvarsCons=Cons(kk);
        if kk=='TaulCalib'
            bench=Cons('BAU');
        else
            bench=Cons('LF');
        end

        perdif=(allvarsCons-bench)./bench*100;
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdif(find(varlist==varr),1:T),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'}, {'color'}, {'k'}  )  

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
            if lablab ==0
                path=sprintf('figures/all_%s/PerdifNoTauf_regime%d_%s_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f.png',date, plotts.regime, kk, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa);
            else
                path=sprintf('figures/all_%s/PerdifNoTauf_Equlab_regime%d_%s_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f.png',date, plotts.regime, kk, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
end
%% effect of tauf by taul in percent
if plotts.perDif_notauf_compTaul==1
    
    fprintf('plott effect tauf by taul')
for lablab =1
    if lablab ==1
        Cons=OTHERPOLL{plotts.regime+1};
    else
        Cons=OTHERPOLLEL{plotts.regime+1};
    end
        allvarsConsTAUL=Cons('TaulCalib');
        allvarsConsTAUL0=Cons('Taul0');

        benchTAUL=Cons('BAU');
        benchTAUL0=Cons('LF');

        perdifCalib=(allvarsConsTAUL-benchTAUL)./benchTAUL*100;
        perdif0=(allvarsConsTAUL0-benchTAUL0)./benchTAUL0*100;
    for lgdind=0

    for l ="Res"%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdifCalib(find(varlist==varr),1:T),time,perdif0(find(varlist==varr),1:T),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k';'k'}  )  

           xticks(txx)
           xlim([1, time(end-1)])
         
            if lgdind==1
                   lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
           if (varr == "F" && indic.noknow_spill==3) || varr == "sff"...
                ||(varr=="EY" && indic.noknow_spill==1 && lablab==0) || (varr=="EY" && indic.noknow_spill==0 && lablab==1)...
                || (varr == "snS" && indic.noknow_spill==1 && lablab==0) ||  (varr == "sn" && indic.noknow_spill==1  && lablab==1)...
                || (varr == "sn" && indic.noknow_spill==1 && lablab==0)  
                ytickformat('%.1f')
           elseif varr == "G" || varr == "sg" || varr == "GFF" ...
                || (varr == "F" && indic.noknow_spill==1 && lablab==0) || (varr == "F" && indic.noknow_spill==0 && lablab==1)...
                || (varr == "F" && indic.noknow_spill==1 && lablab==1) || (varr == "snS" && indic.noknow_spill==0 && lablab==1) 
                
                ytickformat('%.0f')  
           elseif (varr == "snS" && indic.noknow_spill==1 && lablab==1)
                ytickformat('%.1f')  
           elseif (varr == "sn" && indic.noknow_spill==3  && lablab==1)
                ytickformat('%.2f')
           elseif (varr == "snS" && indic.noknow_spill==0 && lablab==0) || varr == "hhhl" 
                ytickformat('%.3f') 
            else
                ytickformat('%.2f')
            end
            xticklabels(Year10)
            if lablab ==0
                path=sprintf('figures/all_%s/PerdifNoTauf_regime%d_CompTaul_%s_Sun%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, plotts.regime, varr, indic.Sun, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            else
                path=sprintf('figures/all_%s/PerdifNoTauf_Equlab_regime%d_CompTaul_%s_Sun%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, plotts.regime, varr,indic.Sun, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        
    end
    end
end
end

%% effect of tauf in levels

if plotts.tauf_notauf==1
for lablab=0:1  
    fprintf('plott effect tauf in levels')

    if lablab==0
        Cons=OTHERPOLL{plotts.regime+1};
    else
        Cons=OTHERPOLLEL{plotts.regime+1};
    end
for k={'TaulCalib', 'Taul0'}% keys are with and without taul
        kk=string(k);
        allvarsCons=Cons(kk);
        if kk=='TaulCalib' % choose correct benchmark
            bench=Cons('BAU');
        else
            bench=Cons('LF');
        end
 
for lgdind=0:1

for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time, allvarsCons(find(varlist==varr),1:T),time,bench(find(varlist==varr),1:T),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k';grrey}  )  

           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            if lgdind==1
                   lgd=legend('with $\tau_f=SCC$', '$\tau_f=0$',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if lablab==0
                path=sprintf('figures/all_%s/LevTaufNoTauf_%s_regime%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, kk, plotts.regime, varr, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            else
                path=sprintf('figures/all_%s/LevTaufNoTauf_%s_Equlab_regime%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, kk, plotts.regime, varr, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
end
end
end
end
end
%% Allocation tauf in levels with and without taul
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.compTauf_Lev==1
    
    fprintf('plott comp tauf by preexisting taul')
   for reg=plotts.regime%[0,2,4,7]
   for lablab=0
           if lablab ==0
               allvars=OTHERPOLL{reg+1};
           else
               allvars=OTHERPOLLEL{reg+1};
           end
        allvarsTaulCalib=allvars('TaulCalib');
        allvarsTaul0=allvars('Taul0');
        
    for lgdind=1
    for l ="Add"% keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if varr=="Emnet"
                main=plot(time,allvarsTaulCalib(find(varlist==varr),1:T), time, allvarsTaul0(find(varlist==varr),1:T),time, Ems);   
                set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; 'k'; 'k'}, {'LineWidth'}, {1.1; 1.1; 1} )   
              
            else
                main=plot(time,allvarsTaulCalib(find(varlist==varr),1:T), time, allvarsTaul0(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'k'} )   
               
            end
           xticks(txx)
           xlim([1, time(end-1)])
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey, 'HandleVisibility','off')

         
                 if lgdind==1
                    if varr=="Emnet"
                        lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$', 'net emission limit','',  'Interpreter', 'latex');
                        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                    else
                       lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$' ,'',  'Interpreter', 'latex');
                       set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                    end
                 end

            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf" || varr=="gAg" || varr=="Emnet" || varr=="GFF" 
                ytickformat('%.0f')
            elseif varr == "gAn"
                ytickformat('%.1f')
            else
                ytickformat('%.2f')
            end
            xticklabels(Year10)
            if lablab ==0
                path=sprintf('figures/all_%s/CompTauf_bytaul_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            else
                path=sprintf('figures/all_%s/CompTauf_bytaul_Equlab_Reg%d_%s_Sun%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.Sun, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            end                
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
   end
   end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Difference in percent between regimes with and without taul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Allocation tauf in levels with and without taul
% compare effect of tauf with and without taul (relevant in model where
% fossil tax is endoegenous, i.e. limitLF=1 ) 
if plotts.compTauf_PER==1
    
    fprintf('plott comp tauf by preexisting taul in percent')
   for reg=plotts.regime%[0,2,4,7]
   for lablab=1
           if lablab ==0
               allvars=OTHERPOLL{reg+1};
           else
               allvars=OTHERPOLLEL{reg+1};
           end
        allvarsTaulCalib=allvars('TaulCalib');
        allvarsTaul0=allvars('Taul0');
        
        %percentage difference between with and without taul allocation
        Perdif = 100*(allvarsTaulCalib-allvarsTaul0)./allvarsTaul0;
        
    for lgdind=0
    for l ="Add"% keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                main=plot(time,Perdif(find(varlist==varr),1:T),  'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
                if lgdind==1
                   lgd=legend('percentage difference with taul to no taul', 'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                end
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey, 'HandleVisibility','off' )
           xlim([1, time(end-1)])

           
            ax=gca;
            ax.FontSize=13;
            if varr=="gAg" || varr =="sgsff"|| varr=="Tauf" 
                ytickformat('%.0f')
            elseif varr == "gAn" 
                ytickformat('%.1f')
            else
                ytickformat('%.2f')
            end
            xticklabels(Year10)
            if lablab ==0
                path=sprintf('figures/all_%s/CompTaufPER_bytaul_Reg%d_%s_Sun%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr,indic.Sun, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            else
                path=sprintf('figures/all_%s/CompTaufPER_bytaul_Equlab_Reg%d_%s_Sun%d_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr,indic.Sun, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
            end                
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
   end
   end
end
%% contrast to no spillover
% fossil tax is endoegenous, i.e. limitLF=1 ) 
if plotts.compTauf_Lev_NK==1
    
    fprintf('plott comp tauf by preexisting taul in levels comparison NK')
   for reg=plotts.regime%[0,2,4,7]
   
        allvars=OTHERPOLL{reg+1};
        allvarsTaulCalib=allvars('TaulCalib');
        allCAlibNK=RES_NK('TaulCalib');
  
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                main=plot(time,allvarsTaulCalib(find(varlist==varr),1:T), time,allCAlibNK(find(varlist==varr),1:T),  'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
                if lgdind==1
                   lgd=legend('with knowledge spillovers', 'no knowledge spillovers', 'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                end
           xticks(txx)
           xlim([1, time(end-1)])

           
            ax=gca;
            ax.FontSize=13;
            if varr =="Tauf" || varr=="gAg"
                 ytickformat('%.0f')
            elseif varr =="gAf" || varr=="GFF"
                 ytickformat('%.1f')                 
            else
                 ytickformat('%.2f')
            end
            xticklabels(Year10)
        path=sprintf('figures/all_%s/CompTauf_bytaul_KN_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
   end
end
%% contrast to no spillover
% fossil tax is endoegenous, i.e. limitLF=1 ) 
if plotts.compTauf_PER_NK==1
    
    fprintf('plott comp tauf by preexisting taul in percent comparison NK')
   for reg=plotts.regime%[0,2,4,7]
   
        allvars=OTHERPOLL{reg+1};
        allvarsTaulCalib=allvars('TaulCalib');
        allvarsTaul0=allvars('Taul0');
        
        allCAlibNK=RES_NK('TaulCalib');
        all0NK=RES_NK('Taul0');
        %percentage difference between with and without taul allocation
        Perdif = 100*(allvarsTaulCalib-allvarsTaul0)./allvarsTaul0;
        PerdifNK = 100*(allCAlibNK-all0NK)./all0NK;
        
    for lgdind=0:1
    for l ="Add" %keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                main=plot(time,Perdif(find(varlist==varr),1:T), time,PerdifNK(find(varlist==varr),1:T),  'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
                if lgdind==1
                   lgd=legend('with knowledge spillovers', 'no knowledge spillovers', 'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
                end
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            if varr =="Tauf" || (varr=="gAg" && indic.limit_LF==1)
                 ytickformat('%.0f')

            elseif (varr =="gAf" && indic.limit_LF==1)|| (varr=="GFF"&& indic.limit_LF==1)
                 ytickformat('%.1f')
            elseif    (varr=="gAg" && indic.limit_LF==0) || (varr=="GFF" && indic.limit_LF==0)
                 ytickformat('%.3f')
            else
                 ytickformat('%.2f')
            end
            xticklabels(Year10)
        path=sprintf('figures/all_%s/CompTaufPER_bytaul_KN_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
   end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Effect of tauf: compariosn BAU and LF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LF and BAU by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.LF_BAU==1
    
    fprintf('plott comp taul LF and BAU')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLL{reg+1};
        
        allvarsBAU=allvars('BAU');
        allvarsLF=allvars('LF');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if varr=="Emnet"
                main=plot(time,allvarsBAU(find(varlist==varr),1:T), time, allvarsLF(find(varlist==varr),1:T),time, Ems);   
                set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; grrey; 'k'}, {'LineWidth'}, {1.1; 1.1; 0.8} )   
                if lgdind==1
                   lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$', 'net emission limit',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
            else
                main=plot(time,allvarsBAU(find(varlist==varr),1:T), time, allvarsLF(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; grrey} )   
                if lgdind==1
                   lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
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
            path=sprintf('figures/all_%s/CompTaul_LFBAU_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep, indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% LF and BAU by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.LF_BAU_PER==1
    
    fprintf('plott comp taul LF and BAU')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLL{reg+1};
        
        allvarsBAU=allvars('BAU');
        allvarsLF=allvars('LF');
        Perdif = (allvarsBAU-allvarsLF)./allvarsLF*100;
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                
            main=plot(time,Perdif(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
%             if lgdind==1
%                lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$',  'Interpreter', 'latex');
%                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%             end
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
            path=sprintf('figures/all_%s/CompTaul_LFBAUPer_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_countec%d_GovRev%d_etaa%.2f.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep, indic.count_techgap, indic.GOV,  etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
 
    end
end

%% LF and BAU by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.LF_BAU_equlab==1
    
    fprintf('plott comp taul LF and BAU')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLLEL{reg+1};
        
        allvarsBAU=allvars('BAU');
        allvarsLF=allvars('LF');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if varr=="Emnet"
                main=plot(time,allvarsBAU(find(varlist==varr),1:T), time, allvarsLF(find(varlist==varr),1:T),time, Ems);   
                set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; grrey; 'k'}, {'LineWidth'}, {1.1; 1.1; 0.8} )   
                if lgdind==1
                   lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$', 'net emission limit',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
            else
                main=plot(time,allvarsBAU(find(varlist==varr),1:T), time, allvarsLF(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; grrey} )   
                if lgdind==1
                   lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
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
            path=sprintf('figures/all_%s/CompTaul_Equlab_LFBAU_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep, indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% LF and BAU by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.LF_BAU_PER_equlab==1
    
    fprintf('plott comp taul LF and BAU')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLLEL{reg+1};
        
        allvarsBAU=allvars('BAU');
        allvarsLF=allvars('LF');
        Perdif = (allvarsBAU-allvarsLF)./allvarsLF*100;
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                
            main=plot(time,Perdif(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
%             if lgdind==1
%                lgd=legend('BAU, $\tau_{\iota}=0.181$', 'Laissez-faire, $\tau_{\iota}=0$',  'Interpreter', 'latex');
%                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%             end
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
            path=sprintf('figures/all_%s/CompTaul_Equlab_LFBAUPer_Reg%d_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_countec%d_GovRev%d_etaa%.2f.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill,    indic.sep, indic.count_techgap, indic.GOV,  etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
 
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison redistribution regimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% effect of redistribution regime
if plotts.compRed==1
    
    fprintf('plott comp policy regimes')
  
        LS=OTHERPOLL{5+1};
        Cons=OTHERPOLL{0+1};
        GS=RES_GS;
       
    for k={'TaulCalib', 'Taul0'}%keys(LS) % keys are with and without taul
        kk=string(k);
        allvarsCons=Cons(kk);
        allvarsGS=GS(kk);
        allvarsLS=LS(kk);

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            if varr=="Emnet"
                main=plot(time,(allvarsLS(find(varlist==varr),1:T)), time,allvarsCons(find(varlist==varr),1:T),...
                    time,allvarsGS(find(varlist==varr),1:T), time, Ems);   
                set(main, {'LineStyle'},{'-';'-.'; '--'; ':'}, {'color'}, {'k'; orrange; grrey; 'k'} , {'LineWidth'}, {1.1;1.1;1.1;0.8})  

                if lgdind==1
                   lgd=legend('lump-sum transfers', 'consolidated budget' , 'green subsidies', 'net emissions limit',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
            else
                main=plot(time,(allvarsLS(find(varlist==varr),1:T)), time,allvarsCons(find(varlist==varr),1:T),...
                time,allvarsGS(find(varlist==varr),1:T),'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'-.'; '--'}, {'color'}, {'k'; orrange; grrey}  )  

                if lgdind==1
                   lgd=legend('lump-sum transfers', 'consolidated budget' , 'green subsidies',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
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
            path=sprintf('figures/all_%s/CompRed_%s_%s_spillover%d_knspil%d_nsk%d_xgr%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers,indic.noknow_spill, plotts.nsk, plotts.xgr, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

end     