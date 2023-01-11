function []=plottsSP_tidiedUp(list, T, etaa, weightext,indic, params, Ems, plotts, percon, MOM)

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
syms analyTaul Hagg PV CEVv CEVvPV CEVvDy Tauf dTaulHh dTaulHl dTaulS dTaulAv dTaulAvS AgAf sffsg sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul Hagg PV Tauf dTaulHh dTaulHl dTaulS dTaulAv dTaulAvS AgAf sffsg sgsff snS  sffS sgS GFF EY CY hhhl whwl LgLf pgpftf pepn gAagg gAg gAf gAn Utilcon Utillab Utilsci];
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
if indic.count_techgap~=1
for xgr=0:1
for nsk =0:1
    indic.xgrowth=xgr;
    indic.noskill=nsk;
    %- sp solution independent of policy

    helper=load(sprintf('SP_target_1008_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
        indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep,  indic.xgrowth, indic.PV, indic.sizeequ, etaa));
    sp_t=helper.sp_all';
    helper=load(sprintf('SP_notarget_1008_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', ...
        indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth,  indic.PV, indic.sizeequ, etaa));
    sp_not=helper.sp_all';

%- other results
    for i=[0,1,4,5] % loop over policy versions
      
       helper = load(sprintf('BAU_1409_taulZero1_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ0_sep%d_notaul0_countec%d_GovRev%d_etaa%.2f.mat', ...
            indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, indic.sep, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
        LF = helper.COMP';
        
        helper=load(sprintf('OPT_notarget_0509_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',... 
            indic.spillovers,indic.noknow_spill, indic.noskill,i,  indic.sep, indic.xgrowth,indic.PV,  plotts.sizeequ, plotts.GOV, etaa));
        opt_not_notaus=helper.opt_all';
        
        helper=load(sprintf('OPT_target_0509_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', ... 
            indic.spillovers,indic.noknow_spill, indic.noskill, i, indic.sep, indic.xgrowth,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
        opt_t_notaus=helper.opt_all';

        RES = containers.Map({'LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                                { LF, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
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
end
%% counterfac technology

    helper=load('SP_target_1008_countec_spillover0_knspil0_noskill0_sep0_xgrowth0_zero1_PV0_sizeequ7.900000e-01_etaa.mat');
    sp_t=helper.sp_all';
    helper=load('SP_notarget_1008_countec_spillover0_knspil0_noskill0_sep0_extern0_xgrowth0_PV1_sizeequ0_etaa0.79.mat');
    sp_not=helper.sp_all';
RES_count=containers.Map({'SP_T', 'SP_NOT' },...
                                {  sp_t, sp_not});
RES_count=add_vars(RES_count, list, params, indic, list.allvars, symms, MOM);

%% counterfac emission limit
for ii=[4,5]
     helper=load(sprintf('OPT_target_0509_emnet1_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,plotts.nknk, plotts.nsk, ii, indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
     opt_t_notaus=helper.opt_all';
     RES_Ems=containers.Map({'OPT_T_NoTaus', },{opt_t_notaus});
     RES_Ems=add_vars(RES_Ems, list, params, indic, list.allvars, symms, MOM);
     OTHERPOL_Ems{ii+1}=RES_Ems;
end

%% counetrfac NEW calibration: phi=3/4, upper bar Scientists
helper=load(sprintf('SP_target_2112_emnet%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
       indic.targetWhat, indic.spillovers,plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,indic.sizeequ, params(list.params=='etaa')));
sp_t = helper.sp_all';

helper=load(sprintf('SP_notarget_2112_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
    indic.spillovers, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,indic.sizeequ, etaa));
sp_not = helper.sp_all';
%plus%d_0501_
helper=load(sprintf('OPT_target_plus30_0501_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul5_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_t_nt= helper.opt_all';
helper=load(sprintf('OPT_target_plus30_0501_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul4_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_t_wt= helper.opt_all';

helper=load(sprintf('OPT_notarget_plus30_0501_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul5_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_not_nt= helper.opt_all';
helper=load(sprintf('OPT_notarget_plus30_0501_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul4_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_not_wt= helper.opt_all';

% with research subsidies
helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul8_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_t_nt_rs= helper.opt_all';
helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul7_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_t_wt_rs= helper.opt_all';
helper=load(sprintf('OPT_notarget_2112_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul8_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_not_nt_rs= helper.opt_all';
helper=load(sprintf('OPT_notarget_2112_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul7_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
opt_not_wt_rs= helper.opt_all';


% laissez faire and BAU
helper = load(sprintf('BAU_2112_taulZero1_spillovers%d_knspil%d_Sun%d_noskill%d_xgrowth%d_labequ0_sep%d_notaul5_countec%d_GovRev%d_etaa%.2f.mat', ...
    indic.spillovers, plotts.nknk, indic.Sun, plotts.nsk,  plotts.xgr, indic.sep, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
LF = helper.COMP';
helper = load(sprintf('BAU_2112_taulZero0_spillovers%d_knspil%d_Sun%d_noskill%d_xgrowth%d_labequ0_sep%d_notaul5_countec%d_GovRev%d_etaa%.2f.mat', ...
    indic.spillovers, plotts.nknk, indic.Sun, plotts.nsk, plotts.xgr, indic.sep, indic.count_techgap, indic.GOV, params(list.params=='etaa')));
BAU = helper.COMP';

RES_NCalib=containers.Map({'LF', 'BAU', 'SP_T', 'SP_NoT', 'OPT_T_NOTaul', 'OPT_T_WithTaul',  'OPT_NoT_NOTaul', 'OPT_NoT_WithTaul'},...
                            {LF, BAU, sp_t, sp_not, opt_t_nt, opt_t_wt, opt_not_nt, opt_not_wt});
RES_NCalib=add_vars(RES_NCalib, list, params, indic, list.allvars, symms, MOM);
RES_NCalib_RS=containers.Map({'LF', 'BAU', 'SP_T', 'SP_NoT', 'OPT_T_NOTaul', 'OPT_T_WithTaul',  'OPT_NoT_NOTaul', 'OPT_NoT_WithTaul'},...
                            {LF, BAU, sp_t, sp_not, opt_t_nt_rs, opt_t_wt_rs, opt_not_nt_rs, opt_not_wt_rs});
RES_NCalib_RS=add_vars(RES_NCalib_RS, list, params, indic, list.allvars, symms, MOM);
%- only optNewCalib_pol_indic.noskillLn_emnet1_Sun2_spillover0_knspil0_xgr0_nsk0_sep0_extern0_PV1_etaa0.79_lgd0.pngimal tauf
if plotts.xgr==0 
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt1_newCalib_target_emnet1_Sun2_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
         plotts.nknk, indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_taul_T=helper.LF_COUNT';
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt0_newCalib_target_emnet1_Sun2_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
          plotts.nknk, indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_tauf_T=helper.LF_COUNT';
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt1_newCalib_notarget_emnet1_Sun2_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
         plotts.nknk, indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_taul_NoT=helper.LF_COUNT';
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt0_newCalib_notarget_emnet1_Sun2_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
          plotts.nknk, indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa));
    count_tauf_NoT=helper.LF_COUNT';
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt5_newCalib_target_emnet%d_Sun%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
          indic.targetWhat, indic.Sun, plotts.nknk,  indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV,etaa));
    count_tauf_TaulFixed=helper.LF_COUNT';
    helper=load(sprintf('COMPEquN_SIM_0501_taufopt6_newCalib_target_emnet%d_Sun%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',...
          indic.targetWhat, indic.Sun, plotts.nknk,  indic.spillovers, plotts.regime_gov, plotts.nsk, indic.sep, plotts.xgr, indic.PV,etaa));
    count_tauf_TaulFixed_TaufJoint=helper.LF_COUNT';
   
    RES_count_NCalib=containers.Map({'CountOnlyTauf_T', 'CountOnlyTaul_T', 'CountOnlyTauf_NoT', 'CountOnlyTaul_NoT' , 'Count_TaulFixed', 'Count_TaulFixed_TaufJoint'},...
            {count_taul_T, count_tauf_T, count_taul_NoT, count_tauf_NoT, count_tauf_TaulFixed, count_tauf_TaulFixed_TaufJoint});
     RES_count_NCalib=add_vars(RES_count_NCalib, list, params, indic, list.allvars, symms, MOM);
  
end
%% counterfactual model

% no knowledge spillovers 
 helper=load(sprintf('OPT_notarget_0509_spillover%d_knspil1_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.spillovers, plotts.nsk,plotts.regime_gov,  indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
 opt_not_notaus=helper.opt_all';
 helper=load(sprintf('OPT_target_0509_spillover%d_knspil1_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.spillovers, plotts.nsk, plotts.regime_gov, indic.sep, plotts.xgr,indic.PV, plotts.sizeequ, plotts.GOV, etaa));
 opt_t_notaus=helper.opt_all';
 RES_noknspil=containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                             {opt_t_notaus, opt_not_notaus});
 RES_noknspil=add_vars(RES_noknspil, list, params, indic, list.allvars, symms, MOM);

  % higher elasticity of substition eppse==10
  helper=load(sprintf('OPT_target_0509_elas10_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.noknow_spill , plotts.nsk,plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_ela=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_ela=add_vars(RES_ela, list, params, indic, list.allvars, symms, MOM);
  % more elastic labor supply
  helper=load(sprintf('OPT_target_0509_sigmaW1_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.noknow_spill ,plotts.nsk, plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_sigma=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_sigma=add_vars(RES_sigma, list, params, indic, list.allvars, symms, MOM);  
  % income effect dominates! households work more when income reduces
  helper=load(sprintf('OPT_target_0509_theta1_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.noknow_spill ,plotts.nsk,plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_bop=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_bop=add_vars(RES_bop, list, params, indic, list.allvars, symms, MOM);    
  % counterfactual technology gap 
  helper=load(sprintf('OPT_target_countec_0509_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.noknow_spill , plotts.nsk,plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_tec=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_tec=add_vars(RES_tec, list, params, indic, list.allvars, symms, MOM);    
 
  % alternative knowledge spillovers
  %knspil==2 smaller knkn (path dependency more important); expect smaller
  %use of labor tax
  helper=load(sprintf('OPT_target_0509_spillover0_knspil2_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', plotts.nsk,plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_smkn=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_smkn=add_vars(RES_smkn, list, params, indic, list.allvars, symms, MOM);
  % knspil ==3 higher knkn => expect to see higher labor tax/exploitation
  % of fossil research
  helper=load(sprintf('OPT_target_0509_spillover0_knspil3_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', plotts.nsk,plotts.regime_gov, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, etaa));
  opt_t_notaus=helper.opt_all';  
  RES_bgkn=containers.Map({'OPT_T_NoTaus'},{opt_t_notaus});
  RES_bgkn=add_vars(RES_bgkn, list, params, indic, list.allvars, symms, MOM);
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
        h1=  OTHERPOLL{5}; % taul can be used
        h2= OTHERPOLL{6}; % taul cannot be used
    
    end
        
%     [COMP, COMPTable] = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);   
    [COMP, COMPTable] = comp_CEV(RES_NCalib('OPT_T_WithTaul'),RES_NCalib('OPT_T_NOTaul') , varlist, varlist, symms, list, params, T, indic);   

    save(sprintf('Table_CEV_Newcalib_%s_regime%d_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d.mat',date, plotts.regime_gov,  indic.sep, plotts.nsk, etaa, plotts.xgr, indic.PV, indic.extern), 'COMPTable');
end
%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%'LF', 'BAU', 'SP_T', 'SP_NoT' 'OPT_T_NOTaul', 'OPT_T_WithTaul',
%'OPT_NoT_NOTaul', 'OPT_NoT_WithTaul', 'Count_TaulFixed' count_tauf_TaulFixed_TaufJoint


%% comparison taul fixed and carbon tax optimal vs carbon tax jont opt and taul fixed

if plotts.taulFixedtaufJoint_newcalib_polPer==1
    fprintf('plott new calib taul fixed pol per')

    nt= RES_count_NCalib(sprintf("Count_TaulFixed"));
    wt=RES_count_NCalib(sprintf("Count_TaulFixed_TaufJoint"));
    per= 100*(wt-nt)./nt;
        
    joint=RES_NCalib(sprintf("OPT_T_WithTaul"));
    perj= 100*(joint-wt)./nt;
    
   for lgdind=0:1
    for l =["Add"]%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(per(find(varlist==varr),1:T)),time,(perj(find(varlist==varr),1:T)), time,zeros(size(per(find(varlist==varr),1:T))),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; '--'}, {'color'}, {'k';orrange; grrey} )   

           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           if lgdind==1
               lgd=legend('effect $\Delta\ \tau_f$','effect $\Delta\ \tau_\iota$', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
             if varr=="Tauf" || varr=="GFF" || varr =="Hagg"
                 ytickformat('%.2f')
             elseif varr =="dTaulAv" || varr =="sffsg"
                 ytickformat('%.1f')
             else
                 ytickformat('%.1f')
             end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_polTaulFixedTaufJointPer_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr,indic.Sun, indic.targetWhat, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% comparison joint optimal versus taul given only tauf optimal
if plotts.taulFixed_newcalib==1
    fprintf('plott new calib taul fixed')

    for s=["T", "NoT"]
        ss=string(s);
    nt= RES_count_NCalib(sprintf("Count_TaulFixed"));
    wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
    eff=RES_NCalib(sprintf("SP_%s", ss));
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {'b';'k';orrange} )   

            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)           
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           if lgdind==1
               lgd=legend('$\tau_\iota$ fixed', 'benchmark', 'first best',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_eff2polTaulFixed_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, varr,indic.Sun, indic.targetWhat, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison taul fixed policy
if plotts.taulFixed_newcalib_pol==1
    fprintf('plott new calib taul fixed pol levels')

    for s=["T", "NoT"]
        ss=string(s);
    nt= RES_count_NCalib(sprintf("Count_TaulFixed"));
    wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'}, {'color'}, {'b';'k'} )   
 
            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
          
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
            if lgdind==1
               lgd=legend('$\tau_\iota$ fixed', 'benchmark',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_polTaulFixed_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, varr,indic.Sun, indic.targetWhat, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.taulFixed_newcalib_polPer==1
    fprintf('plott new calib taul fixed pol per')

    for s=["T", "NoT"]
        ss=string(s);
    nt= RES_count_NCalib(sprintf("Count_TaulFixed"));
    wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
    per= 100*(wt-nt)./nt;
        
  %  for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(per(find(varlist==varr),1:T)), time,zeros(size(per(find(varlist==varr),1:T))),'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )   
%             if lgdind==1
%                lgd=legend('$\tau_\iota$ fixed', 'benchmark',  'Interpreter', 'latex');
%                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%             end
            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
             if varr=="Tauf" || varr=="GFF"
                 ytickformat('%.2f')
             elseif varr =="dTaulAv" || varr =="sffsg"
                 ytickformat('%.1f')
             else
                 ytickformat('%.1f')
             end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_polTaulFixedPer_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f.png',date,ss, varr,indic.Sun, indic.targetWhat, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
     %   end
        end
    end
    end
end
%% All figures single
if plotts.single_pol_NC==1
    
    fprintf('plotting single graphs')

    %- loop over economy versions
    for nnt=[ "RS"]
        if nnt =="base"
            RES=RES_NCalib;
        else
            RES= RES_NCalib_RS;
        end
    for s=["T", "NoT"]
        ss=string(s);
        wt=RES(sprintf("OPT_%s_WithTaul", ss));
        
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ll=="HH" && varr=="Emnet"
                main=plot(time,wt(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T), 'LineWidth', 1.1);  
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

            else
                main=plot(time,wt(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

            end
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)

           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="dTaulAv" || varr=="dTaulAvS"
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Single_NC_%s_%s_emnet%d_Sun%d_regime%s_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',...
               date, ss, varr, indic.targetWhat, indic.Sun, nnt, indic.spillovers,plotts.nknk, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, plotts.sizeequ, plotts.GOV, etaa);
         
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_LF_newcalib==1
    fprintf('plott new calib levels with lf')

    LFall= RES_NCalib("LF");
    for s=["T", "NoT"]
        ss=string(s);
        nt= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
        wt=RES_NCalib(sprintf("OPT_%s_WithTaul", ss));
        eff=RES_NCalib(sprintf("SP_%s",ss) );
            
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--'; '--';'-'; '--'}, {'color'}, {'k';'b';'k';orrange} )   
            if lgdind==1
               lgd=legend('laissez-faire','carbon tax only', 'benchmark', 'first best',  'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_eff2pol_LF_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',...
                date, ss, varr,indic.Sun, indic.targetWhat, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_effBau_newcalib==1
    fprintf('plott new calib levels  eff with bau')

    LFall= RES_NCalib("BAU");
    for s=["T", "NoT"]
        ss=string(s);
        eff=RES_NCalib(sprintf("SP_%s",ss) );
            
    for lgdind=1
    for l =["Add"]%keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{':'; '--'}, {'color'}, {'k';'k'} )   
    
            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
        
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            elseif varr =="GFF"
                ytickformat('%.0f')
                 ylim([-5, 40.1])
            elseif varr =="sffsg"
                ytickformat('%.1f')
                 ylim([-0.2, 1.21])
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            if lgdind==1
               lgd=legend('business as usual', 'first best',  'Interpreter', 'latex');
               if varr=="GFF"
                   set(lgd, 'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');               
               else
                   set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
               end
            end
            path=sprintf('figures/all_%s/NewCalib_effBAU_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',...
                date, ss, varr,indic.Sun, indic.targetWhat, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_effBauOPt_newcalib==1
    fprintf('plott new calib levels  eff with bau')

    LFall= RES_NCalib("BAU");
    for s=["T", "NoT"]
        ss=string(s);
        eff=RES_NCalib(sprintf("SP_%s",ss) );
        wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{':'; '--'; '-'}, {'color'}, {'k';'k'; orrange} )   
      
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
        
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            elseif varr =="GFF"
                ytickformat('%.0f')
                 ylim([-5, 40.1])
            elseif varr =="sffsg"
                ytickformat('%.1f')
                 ylim([-0.2, 1.21])
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            if lgdind==1
               lgd=legend('business as usual', 'first best', 'optimal policy',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            
            path=sprintf('figures/all_%s/NewCalib_effBauOpt_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',...
                date, ss, varr,indic.Sun, indic.targetWhat, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_effLF_newcalib==1
    fprintf('plott new calib levels  eff with lf')

    LFall= RES_NCalib("LF");
    for s=["T", "NoT"]
        ss=string(s);
        eff=RES_NCalib(sprintf("SP_%s",ss) );
            
    for lgdind=1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{':'; '--'}, {'color'}, {'k';'k'} )   
    
            
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
        
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            elseif varr =="GFF"
                ytickformat('%.0f')
                 ylim([-5, 40.1])
            elseif varr =="sffsg"
                ytickformat('%.1f')
                 ylim([-0.2, 1.21])
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
                    if lgdind==1
               lgd=legend('laissez-faire', 'first-best',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            path=sprintf('figures/all_%s/NewCalib_effLF_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',...
                date, ss, varr,indic.Sun, indic.targetWhat, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_effLFOPt_newcalib==1
    fprintf('plott new calib levels  eff with lf')

    LFall= RES_NCalib("LF");
    for s=["T", "NoT"]
        ss=string(s);
        eff=RES_NCalib(sprintf("SP_%s",ss) );
        wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(LFall(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{':'; '--'; '-'}, {'color'}, {'k';'k'; orrange} )   
      
           xticks(txx)
           xline(7, 'LineStyle', ':', 'LineWidth', 0.8, 'color', grrey)
        
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            elseif varr =="GFF"
                ytickformat('%.0f')
                 ylim([-5, 40.1])
            elseif varr =="sffsg"
                ytickformat('%.1f')
                 ylim([-0.2, 1.21])
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            if lgdind==1
               lgd=legend('laissez-faire', 'first-best', 'optimal policy',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            
            path=sprintf('figures/all_%s/NewCalib_effLFOpt_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',...
                date, ss, varr,indic.Sun, indic.targetWhat, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% comparison parameter alternatives
if plotts.phi_newcalib==1
    fprintf('plott new calib')

    for s=["T", "NoT"]
        ss=string(s);
    nt= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
    wt=RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
    eff=RES_NCalib(sprintf("SP_%s", ss));
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)), time,(eff(find(varlist==varr),1:T)), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {'b';'k';orrange} )   
            if lgdind==1
               lgd=legend('carbon tax only', 'benchmark', 'first best',  'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_eff2pol_%s_%s_Sun%d_emnet%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, varr,indic.Sun, indic.targetWhat, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% no efficient
if plotts.phi_newcalib_noeff==1
    fprintf('plott new calib no eff')

    for s =["NoT", "T"]
        ss=string(s);
    nt = RES_NCalib(sprintf("OPT_%s_NOTaul",ss));
    wt = RES_NCalib(sprintf("OPT_%s_WithTaul", ss));
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(nt(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'}, {'color'}, {'b';'k'} )   
            if lgdind==1
               lgd=legend('carbon tax only', 'benchmark', 'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_pol_%s_%s_emnet%d_Sun%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, varr,indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% no efficient
if plotts.phi_newcalib_TvsNoT==1
    fprintf('plott new calib no T vs NoT')

for nnt=["RS"]
    if nnt=="base"
        RES=RES_NCalib;
    else
        RES=RES_NCalib_RS;
    end
    wt = RES(sprintf("OPT_T_WithTaul"));
    not = RES(sprintf("OPT_NoT_WithTaul"));

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(not(find(varlist==varr),1:T)), time,(wt(find(varlist==varr),1:T)),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'}, {'color'}, {'b';'k'} )   
            if lgdind==1
               lgd=legend('no target', 'with target', 'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/NewCalib_pol_TvsNoT_%s_%s_emnet%d_Sun%d_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,varr,nnt, indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end

%% comparison only carbon tax (reg 5) to only tauf from combined
if plotts.count_devs_both_NC==1
    
    fprintf('plott counterfactual deviation from carbon tax only pol => role of adjustment in tauf')

    %- read in variable container of chosen regime
         for s=["NoT"]
             ss=string(s);
        allvarsntaul= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
        allvarscount=RES_count_NCalib(sprintf("CountOnlyTauf_%s",ss)); % version with only tauf
        perdiftauf= 100*(allvarscount-allvarsntaul)./allvarsntaul;
        
         
        allvars= RES_NCalib(sprintf("OPT_%s_WithTaul",ss));
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
%             if varr== "sgsff"
%                 ytickformat('%.1f')
%             elseif varr == "Hagg"
%                 ytickformat('%.2f')
%                 ylim([-0.12, 0.061])
%             else
%                  ytickformat('%.2f')
%             end
            
            if lgdind==1
               lgd=legend('effect $\tau_F$' , 'effect $\tau_\iota$', 'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 21,'Orientation', 'vertical');
            end
            xticklabels(Year10)

            path=sprintf('figures/all_%s/CountTAUF_Both_Opt_NewCalib_%s_emnet%d_Sun%d_target_%s_nsk%d_xgr%d_knspil%d_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date,ss, indic.targetWhat, indic.Sun, varr ,plotts.nsk, plotts.xgr, plotts.nknk, plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
    end
end

%% optimal with and without taul in levels 
if plotts.comp_OPTPer_NCalib==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage opt and no taul no eff') 
    
   for s=["T", "NoT"]
       ss=string(s);
        allvarsnt= RES_NCalib(sprintf("OPT_%s_NOTaul", ss));
        allvars=RES_NCalib(sprintf("OPT_%s_WithTaul", ss));
        Perdif=100*(allvars-allvarsnt)./allvarsnt;
    %% 
    for l = keys(lisst) % loop over variable groups
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
%             if varr=="SWF" || varr== "sn"
%                ytickformat('%.0f')
%             elseif varr=="sff" || varr=="sg" ||  varr == "sgsff" || varr =="GFF" || varr =="sffsg"
%                ytickformat('%.1f')
%             else
%                ytickformat('%.2f')
%             end
            xticklabels(Year10)
%            if lgdind==1
%               lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
%               set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
%            end
        path=sprintf('figures/all_%s/%s_COMPtaulPerNewCalib_%s_regime%d_emnet%d_Sun%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png',date, varr,ss, plotts.regime_gov,indic.targetWhat, indic.Sun, indic.spillovers, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa);
        
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
    end
end

%% comparison parameter alternatives
if plotts.phi_sens==1
    
    fprintf('plott sensitivity phi')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        ben= RES("OPT_T_NoTaus");
        nokn=RES_noknspil("OPT_T_NoTaus");
        smkn=RES_smkn("OPT_T_NoTaus");
        bgkn=RES_bgkn("OPT_T_NoTaus");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(bgkn(find(varlist==varr),1:T)), time,(ben(find(varlist==varr),1:T)), time,(smkn(find(varlist==varr),1:T)),time,nokn(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {'b';'k';orrange; grrey} )   
            if lgdind==1
               lgd=legend('$\phi=0.6$', '$\phi=0.5$', '$\phi=0.25$','$\phi=0$',  'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Phi_Sens_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,indic.noknow_spill, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% comparison parameter alternatives
if plotts.sens_other==1
    
    fprintf('plott sensitivity other')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        ben= RES("OPT_T_NoTaus");
        sig=RES_sigma("OPT_T_NoTaus");
        bop=RES_bop("OPT_T_NoTaus");
        ela=RES_ela("OPT_T_NoTaus");
        tec = RES_tec("OPT_T_NoTaus");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(ben(find(varlist==varr),1:T)), time,(sig(find(varlist==varr),1:T)), ...
                time,(bop(find(varlist==varr),1:T)),time,ela(find(varlist==varr),1:T),  'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'; '--'; '--'}, {'color'}, {'k';'b';orrange; grrey} )   
            if lgdind==1
               lgd=legend('baseline', '$\sigma=2/3$', '$\theta=2$','$\varepsilon_e=10$',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            
           xticks(txx)
           xlim([1, time(end-1)])

           
            ax=gca;
            ax.FontSize=13;
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.0f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/SensOther_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,indic.noknow_spill, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_nsk_LF==1
    
    fprintf('plott only taul nsk model devLF')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1}; 

        allvarslf= RES("LF");
        allvarstlf= RESnsk("LF");
        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('nsk');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarslf(find(varlist==varr),1:T),...
                time,allvarst(find(varlist==varr),1:T), time, allvarstlf(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';':'; '--'; ':'}, {'color'}, {'k'; orrange; 'b';grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'benchmark, LF', 'homogeneous skill', 'homogeneous skill, LF',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountNskTaulLF_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_xgr_LF==1
    
    fprintf('plott only taul xgr model devLF')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1}; 

        allvarslf= RES("LF");
        allvarstlf= RESxgr("LF");
        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('xgr');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarslf(find(varlist==varr),1:T),...
                time,allvarst(find(varlist==varr),1:T), time, allvarstlf(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';':'; '--'; '--'}, {'color'}, {'k';orrange; 'b';grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'benchmark, LF', 'exogenous growth', 'exogenous growth, LF',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountXgrTaulLF_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% Comparison Emission limits
if plotts.Ems_limit==1
    
    fprintf('plott levels with other emission limit')

        RES=OTHERPOL{plotts.regime_gov+1}; 
        RES_ems=OTHERPOL_Ems{plotts.regime_gov+1}; 

        allvars= RES("OPT_T_NoTaus");
        allvarsems= RES_ems("OPT_T_NoTaus");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(allvars(find(varlist==varr),1:T)), time,allvarsems(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k';orrange} )   
            if lgdind==1
               lgd=legend('benchmark' , 'equal \% reduction' ,  'Interpreter', 'latex');
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
            if varr=="Tauf"
                ytickformat('%.0f')
            elseif varr =="dTaulAv"
                ytickformat('%.1f')
            else
                ytickformat('%.1f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/Ems_Sens_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% policy decomposition with less strict emission limit; comparison to benchmark decomposition
if plotts.Ems_limit_Decomp==1
    
    fprintf('plott decomposition ems limit')

        RES=OTHERPOL{plotts.regime_gov+1};
        RES_ems=OTHERPOL_Ems{plotts.regime_gov+1}; 

        if plotts.regime_gov==4
            RESnt=OTHERPOL{plotts.regime_gov+1+1}; 
            RES_emsnt=OTHERPOL_Ems{plotts.regime_gov+1+1};
        end
        allvars= RES("OPT_T_NoTaus");
        allvarsnt= RESnt("OPT_T_NoTaus");
        allvarsems= RES_ems("OPT_T_NoTaus");
        allvarsemsnt= RES_emsnt("OPT_T_NoTaus");

        Perdif=100*(allvars-allvarsnt)./allvarsnt;
        Perdifems=100*(allvarsems-allvarsemsnt)./allvarsemsnt;
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,(Perdif(find(varlist==varr),1:T)), time,Perdifems(find(varlist==varr),1:T),time,zeros(1,T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-'; '--'; '--'}, {'color'}, {'k';orrange; grrey} )   
            if lgdind==1
               lgd=legend('benchmark' , 'equal \% reduction' ,  'Interpreter', 'latex');
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
            if varr=="sffsg"
                ytickformat('%.0f')
            elseif varr=="Tauf"
                ytickformat('%.1f')                
            else
                ytickformat('%.0f')
            end
            xticklabels(Year10)
            path=sprintf('figures/all_%s/EmsDecomp_CTO_Sens_%s_spillover%d_knspil%d_xgr%d_nsk%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers,plotts.nknk, plotts.xgr, plotts.nsk, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% effect of taul in model without xgr and with xgr
% compare effect of only taul in benchmark and in exogenous growth model
if plotts.count_taul_xgr_lev==1
    
    fprintf('plott only taul xgr model level')

        allvars= RES_count_SAMETAUL_onlytaul("test");
        allvarst=RES_count_SAMETAUL_onlytaul('xgr');
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('benchmark' , 'exogenous growth',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountXgrTaul_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev==1
    
    fprintf('plott only tauf model level')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars= RES("LF");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only environmental tax',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountTauf_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev_Ben==1
    
    fprintf('plott only tauf model level and benchmark pol')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvarsLF= RES("LF");
        allvars=RES("OPT_T_NoTaus");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarsLF(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T),time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-.';'--'; '-'}, {'color'}, {grrey; 'b'; 'k'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only environmental tax', 'both taxes optimal', 'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountTauf_Ben_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_tauflev_Ben_noLF==1
    
    fprintf('plott only tauf model level and benchmark pol no lf')

    %- loop over economy versions
    for k=keys(RES_count_onlytauf)
        kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars=RES("OPT_T_NoTaus");
        allvarst= RES_count_onlytauf(kk);
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvarst(find(varlist==varr),1:T),time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'--'; '-'}, {'color'}, {'b'; 'k'} )   
            if lgdind==1
               lgd=legend( 'only environmental tax', 'both taxes optimal', 'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountTauf_Ben_noLF_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end
%% effect tauf relative to laissez faire
% in levels
if plotts.count_taullev==1
    
    fprintf('plott only taul model level')
%- loop over economy versions
    for k=keys(RES_count_onlytaul)
      kk=string(k);
            %- read in variable container of chosen regime
            if kk=="test"
                RES=OTHERPOL{plotts.regime_gov+1}; 
            elseif kk=="nsk"
                RES=OTHERPOL_nsk{plotts.regime_gov+1}; 
            elseif kk=="xgr"
                RES=OTHERPOL_xgr{plotts.regime_gov+1}; 
            elseif kk=="xgr_nsk"
                RES=OTHERPOL_xgr_nsk{plotts.regime_gov+1}; 
            end
        allvars= RES("LF");
        allvarst=RES_count_onlytaul(kk);

    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
     
            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarst(find(varlist==varr),1:T), 'LineWidth', 1.1);   
            set(main, {'LineStyle'},{'-';'--'}, {'color'}, {'k'; 'b'} )   
            if lgdind==1
               lgd=legend('laissez-faire' , 'only income tax',  'Interpreter', 'latex');
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
            path=sprintf('figures/all_%s/CountTaul_mod%s_target_%s_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, kk, varr, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end



%% Comparison model versions in one graph
% in levels
if plotts.count_modlev==1
    
    fprintf('plott counterfactual model level')

    %- read in variable container of chosen regime
    RES=OTHERPOL{plotts.regime_gov+1};
%     RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
%     RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};  
    %- loop over economy versions
        allvars= RES("OPT_T_NoTaus");
        allvarsnsk=RES_count("nsk");
        allvarsxgr=RES_count("xgr");
        
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsxgr(find(varlist==varr),1:T) ,time,allvarsnsk(find(varlist==varr),1:T), 'LineWidth', 1.1);   
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
            path=sprintf('figures/all_%s/CountMod1_target_%s_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr , plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end

%% Comparison model versions: in deviation from efficient 1 graph
% in levels
if plotts.count_modlev_eff
    
    fprintf('plotting comp nsk xgr percent from efficient 1 counterfact')

    %- read in variable container of chosen regime
        RES=OTHERPOL{plotts.regime_gov+1};
        RESnsk=OTHERPOL_nsk{plotts.regime_gov+1};
        RESxgr=OTHERPOL_xgr{plotts.regime_gov+1};
        
        allvars= RES("OPT_T_NoTaus");
        allvarseff=RES("SP_T");
        allvarsnsk=RES_count("nsk");
        allvarsnskeff=RESnsk("SP_T");
        allvarsxgr=RES_count("xgr"); % counterfactual
        allvarsxgreff=RESxgr("SP_T");
        
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
            path=sprintf('figures/all_%s/Per1_CountMod_%s_regime%d_spillover%d_sep%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, varr , plotts.regime_gov, indic.spillovers, indic.sep,indic.extern, indic.PV, etaa, lgdind);
    
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
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
        
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));

            main=plot(time,allvars(find(varlist==varr),1:T),time,allvarsxgr(find(varlist==varr),1:T) ,time,allvarsnsk(find(varlist==varr),1:T), 'LineWidth', 1.1);   
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
            path=sprintf('figures/all_%s/CompMod1_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, etaa, lgdind);
    
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

%% All figures single
if plotts.single_pol==1
    
    fprintf('plotting single graphs')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};
%     RES=ccc;
  

    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lisst) % loop over variable groups
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
           if ismember(varr, list.growthrates)
                xlim([1, time(end-1)])
           else             
                xlim([1, time(end)])
           end
           
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if indic.count_techgap==0
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date,  ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, indic.PV, plotts.sizeequ, plotts.GOV,   etaa);
            else
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_countec_PV%d_sizeequ%d_GOV%d_etaa%.2f.png',date,  ii, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV,plotts.sizeequ, plotts.GOV, etaa);
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
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,indic.extern, indic.PV,  etaa, lgdind);
        else
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_countec_extern%d_PV%d_etaa%.2f_lgd%d.png',date,  plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV, etaa, lgdind);
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
                    lgd=legend('benchmark', 'integrated policy, with income tax', 'Interpreter', 'latex');
                elseif count == 1 % integrated without income tax
                     lgd=legend('benchmark', 'integrated policy, no income tax', 'Interpreter', 'latex');
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
    ssr= string({'SP_T', 'SP_NOT' });%,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'

    for i =[1] % to filter benchmark: policy with target
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));
        if indic.count_techgap==0
                RES=OTHERPOLL{plotts.regime_gov+1};

                 allvars= RES(t);
                 allvarsnot=RES(nt); 
        else
                 allvars =RES_count(t);
                 allvarsnot=RES_count(nt); 
        end

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time,allvars(find(varlist==varr),1:T), time,allvarsnot(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('wih emission limit', 'no emission limit', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
           end
           if indic.count_techgap==0
                path=sprintf('figures/all_%s/%s_TargetComp%s_regime%d_knspil%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, varr, ii, plotts.regime_gov,indic.noknow_spill, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
           else
                path=sprintf('figures/all_%s/%s_TargetComp%s_countec_knspil0_spillover0_noskill0_sep0_xgrowth0_PV1_etaa%.2f_lgd%d.png', date, varr, ii, etaa, lgdind);
           end
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
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    
    for withlff=0
        lff=RES('LF');
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        varl=varlist;

         %- benchmark policy
         allvars= RES(io);
         allvarseff=RES(ie); 

     %- comparison         
     for nt = 3% bb % loop over policy scenarios but benchmark
        RES_help=OTHERPOLL{nt};
        count=nt-1;
        allvarsnotaul =RES_help(io);

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varl==varr),1:T),  time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot( time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
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
                        lgd=legend('laissez-faire',  'benchmark policy',  'integrated policy, with income tax', 'efficient','Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire','benchmark policy', 'integrated policy, no income tax',   'efficient', 'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'with income tax', 'without income tax', 'efficient',  'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'benchmark policy', 'no redistribution, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'benchmark policy', 'lump-sum transfers, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire',  'benchmark policy', 'lump-sum transfers, no income tax', 'efficient',  'Interpreter', 'latex');                        
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if count ==0
                        lgd=legend(  'benchmark policy', 'integrated policy, with income tax',  'efficient', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'benchmark policy', 'integrated policy, no income tax', 'efficient',   'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend(  'with income tax', 'without income tax',  'efficient', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend(  'benchmark policy', 'no redistribution, with income tax', 'efficient',  'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'benchmark policy', 'lump-sum transfers, with income tax',  'efficient', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'benchmark policy', 'lump-sum transfers, no income tax', 'efficient',  'Interpreter', 'latex');                        
                   end
                   
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_benchregime%d_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png',date, varr, io, plotts.regime_gov, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind, withlff);
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
    fprintf('plotting efficient')   

    %- read in container of results: any fine for social planner
    
    eff= ["SP_T" "SP_NOT"];   
    for i =[1,2]

        ie=eff(i);
        if indic.count_techgap==0
            
            RES=OTHERPOLL{plotts.regime_gov+1};
            allvarseff=RES(ie); 
        else
            allvarseff =RES_count(ie);
        end

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           
           main=plot( time,allvarseff(find(varlist==varr),1:T));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end-1)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
              lgd=legend('efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
           if indic.count_techgap==0
                path=sprintf('figures/all_%s/%s_CompEff%s_onlyeff_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d.png', date, varr, ie, indic.spillovers,indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind);
           else
                path=sprintf('figures/all_%s/%s_CompEff%s_count_onlyeff_spillover0_knspil0_noskill0_sep0_xgrowth0_PV1_etaa%.2f_lgd%d.png', date, varr, ie, etaa, lgdind);

           end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
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

%% level differences: NEW CALIBRATION with and without taul
if plotts.comp_OPT_NCAlib==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt and no taul no eff') 

        allvars= RES_NCalib('OPT_T_WithTaul');
        allvarsnt =RES_NCalib('OPT_T_NOTaul');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
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
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_OPT_COMPtaul_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nknk, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
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
        RESnt =OTHERPOLL{3}; % version without taul
     end
        allvars= RES('OPT_T_NoTaus');
        allvarsnt =RESnt('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvars(find(varlist==varr),1:T)),time,(allvarsnt(find(varlist==varr),1:T)),'LineWidth', 1.1);            
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
              lgd=legend('with income tax', 'without income tax', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_OPT_COMPtaul_regime%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, indic.noknow_spill, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
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

     if plotts.regime_gov==3
        RES=OTHERPOLL{plotts.regime_gov+1};
     end
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
%% optimal with and without taul in levels 
if plotts.comp_OPT_NK==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting levels opt benchmark and no knowledge spillovers') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        allvars= RES('OPT_T_NoTaus');
        allvarsNK =RES_noknspil('OPT_T_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot( time,(allvarsNK(find(varlist==varr),1:T)),time,(allvars(find(varlist==varr),1:T)),'LineWidth', 1.1);            
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
              lgd=legend('no knowledge spillovers', 'benchmark', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_KN_FullMod_sizeequ%d_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_GOV%d_etaa%.2f_lgd%d.png',date, varr, plotts.sizeequ, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,plotts.GOV,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
   end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
    RES=OTHERPOLL{plotts.regime_gov+1};

    for withlff=0
         lff=RES('LF');
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        
        allvars= RES(io);
        varl=varlist;
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
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
        path=sprintf('figures/all_%s/%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_PV%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
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
for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        bau=RES('LF');

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
                  lgd=legend('Social planner', 'laissez-faire', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'laissez-faire', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_LFComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, ii, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
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

        path=sprintf('figures/all_%s/%s_PercentageLFDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV,  etaa, lgdind);
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
        RESnt =OTHERPOLL{3}; 
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
           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'b'; orrange} )   
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
              lgd=legend('with income tax', 'without income tax', 'efficient', 'Interpreter', 'latex');
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
        RESnt =OTHERPOLL{3}; 
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
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,time,(allvarsnt(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T)*100,'LineWidth', 1.1);            
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