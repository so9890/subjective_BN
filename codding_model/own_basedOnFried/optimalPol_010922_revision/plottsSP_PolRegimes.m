function []=plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon)

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
syms analyTaul PV CEVv CEVvPV CEVvDy AgAf sgsff snS GFF EY CY hhhl whwl LgLf pgpftf pepn gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg  S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul PV AgAf sgsff snS  GFF EY CY hhhl whwl LgLf pgpftf pepn gAagg gAg gAf gAn Utilcon Utillab Utilsci];
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
for lablab =0:1 % with equal and non-equal labor supply
    indic.labshareequ=lablab; % relevant for additional variables
for xgr=0:1
for nsk =0:1
    indic.xgrowth=xgr;
    indic.noskill=nsk;
%- other results
    for i=[0,5] % loop over policy versions
        helper=load(sprintf('COMP_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth,lablab, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
        all_TaulC=helper.COMP';
        helper=load(sprintf('COMP_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth,lablab, indic.sep, i,indic.limit_LF, indic.emsbase, indic.count_techgap, indic.GOV,  etaa));
        all_Taul0=helper.COMP';
        helper=load(sprintf('TAUF_1409_taulZero0_knspil%d_limit%d_EmsBase%d_xgr%d_nsk%d_labequ%d_countec%d_GovRev%d',...
            indic.noknow_spill, indic.limit_LF,indic.emsbase, indic.xgrowth, indic.noskill,lablab, indic.count_techgap, indic.GOV), 'TAUF');
        tauff_TaulC=helper.TAUF; % tauff is independent of notaul!
        helper=load(sprintf('TAUF_1409_taulZero1_knspil%d_limit%d_EmsBase%d_xgr%d_nsk%d_labequ%d_countec%d_GovRev%d',...
           indic.noknow_spill, indic.limit_LF,indic.emsbase, indic.xgrowth, indic.noskill,lablab, indic.count_techgap, indic.GOV), 'TAUF');
        tauff_Taul0=helper.TAUF; % tauff is independent of notaul!
        helper = load(sprintf('BAU_1409_taulZero0_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_countec%d_GovRev%d_etaa%.2f.mat',...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, lablab, indic.sep, indic.notaul, indic.count_techgap, indic.GOV,  params(list.params=='etaa')));
        BAU = helper.COMP';
        helper = load(sprintf('BAU_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_countec%d_GovRev%d_etaa%.2f.mat', ...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, lablab, indic.sep, indic.notaul, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
        LF = helper.COMP';
        RES = containers.Map({'TaulCalib', 'Taul0', 'BAU', 'LF'},{all_TaulC, all_Taul0, BAU, LF});
        %- add additional variables
        if lablab ==0
            if xgr==0 && nsk==0
                OTHERPOL{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
                TAUFS=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==0 && nsk==1
                OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
                TAUFS_nsk=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==1 && nsk==0
                OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
                TAUFS_xgr=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            elseif xgr==1 && nsk==1
                OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
                TAUFS_xgr_nsk=containers.Map({'TaulCalib', 'Taul0'},{tauff_TaulC, tauff_Taul0});
            end
        else
            if nsk==0 && xgr ==0
                OTHERPOL_EQULab{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
            elseif nsk==1 && xgr ==0
                OTHERPOL_EQULab_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
             elseif xgr==1 && nsk==0
                 OTHERPO_EQULabL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
            elseif xgr==1 && nsk==1
                OTHERPOL_EQULab_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
            end

        end
    end

end
end
end

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
%% effect of tauf single regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.perDif_notauf==1
    
    fprintf('plott effect tauf by policy regime')
  

        Cons=OTHERPOLL{plotts.regime+1};
       
        
    for k={'TaulCalib', 'Taul0'}%keys(LS) % keys are with and without taul
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
            path=sprintf('figures/all_%s/PerdifNoTauf_regime%d_%s_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f.png',date, plotts.regime, kk, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% effect of tauf single regime
if plotts.perDif_notauf_compTaul==1
    
    fprintf('plott effect tauf by policy regime')
  

        Cons=OTHERPOLL{plotts.regime+1};
       
        
        allvarsConsTAUL=Cons('TaulCalib');
        allvarsConsTAUL0=Cons('Taul0');

        benchTAUL=Cons('BAU');
        benchTAUL0=Cons('LF');

        perdifCalib=(allvarsConsTAUL-benchTAUL)./benchTAUL*100;
        perdif0=(allvarsConsTAUL0-benchTAUL0)./benchTAUL0*100;
    for lgdind=0:1

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,perdifCalib(find(varlist==varr),1:T),time,perdif0(find(varlist==varr),1:T),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k';grrey}  )  

           xticks(txx)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            if lgdind==1
                   lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            path=sprintf('figures/all_%s/PerdifNoTauf_regime%d_CompTaul_%s_spillover%d_nsk%d_xgr%d_knspil%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, plotts.regime, varr, indic.spillovers, plotts.nsk, plotts.xgr,indic.noknow_spill, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        
    end
    end
    end

%% effect of redistribution regime
if plotts.compRed==1
    
    fprintf('plott comp policy regimes')
  
        LS=OTHERPOLL{4+1};
        Cons=OTHERPOLL{0+1};
        GS=OTHERPOLL{7+1};
       
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
            path=sprintf('figures/all_%s/CompRed_%s_%s_spillover%d_nsk%d_xgr%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% Effect Allocation taul by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.compTaul_Red==1
    
    fprintf('plott comp taul by policy regime')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLL{reg+1};
        
        allvarsTaulCalib=allvars('TaulCalib');
        allvarsTaul0=allvars('Taul0');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if varr=="Emnet"
                main=plot(time,allvarsTaulCalib(find(varlist==varr),1:T), time, allvarsTaul0(find(varlist==varr),1:T),time, Ems);   
                set(main, {'LineStyle'},{'-';'--'; ':'}, {'color'}, {'k'; grrey; 'k'}, {'LineWidth'}, {1.1; 1.1; 0.8} )   
                if lgdind==1
                   lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$', 'net emission limit',  'Interpreter', 'latex');
                    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
                end
            else
                main=plot(time,allvarsTaulCalib(find(varlist==varr),1:T), time, allvarsTaul0(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; grrey} )   
                if lgdind==1
                   lgd=legend('$\tau_{\iota}=0.181$', '$\tau_{\iota}=0$' ,  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CompTaul_Reg%d_%s_spillover%d_nsk%d_xgr%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, reg, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% LF and BAU by policy regime
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.LF_BAU==1
    
    fprintf('plott comp taul LF and BAU')
   for reg=plotts.regime%[0,2,4,7]
        allvars=OTHERPOLL{reg+1};
        
        allvarsBAU=allvars('BAU');
        allvarsLF=allvars('LF');
        
    for lgdind=0:1
    for l ="Add"%keys(lisst) % loop over variable groups
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
        
    for l ="Add"%keys(lisst) % loop over variable groups
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
    for l ="Add"%keys(lisst) % loop over variable groups
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
        
    for l ="Add"%keys(lisst) % loop over variable groups
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Below not yet updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% effect of taul in model percentage deviation
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.compRed_TaulPer==1
    
    fprintf('plott comp taul all policy regimes percent')
  
        LS=OTHERPOLL{4+1};
        Cons=OTHERPOLL{0+1};
        GS=OTHERPOLL{7+1};
        NR=OTHERPOLL{2+1};
       
        DifLS  =(LS('TaulCalib')-LS('Taul0'))./LS('Taul0').*100;
        DifCons=(Cons('TaulCalib')-Cons('Taul0'))./Cons('Taul0').*100;
        DifGS  =(GS('TaulCalib')-GS('Taul0'))./GS('Taul0').*100;
        DifNR  =(NR('TaulCalib')-NR('Taul0'))./NR('Taul0').*100;

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,DifLS(find(varlist==varr),1:T), time,DifCons(find(varlist==varr),1:T),...
                time,DifGS(find(varlist==varr),1:T), time, DifNR(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'-.'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'; grrey} )   
            if lgdind==1
               lgd=legend('lump-sum transfers', 'consolidated budget' , 'green subsidies', 'no redistribution',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CompRedTaulPer_%s_spillover%d_nsk%d_xgr%d_sep%d_LFlimit%d_emsbase%d_countec%d_GovRev%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, plotts.nsk, plotts.xgr, indic.sep,indic.limit_LF, indic.emsbase,  indic.count_techgap, indic.GOV,  etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end

end     