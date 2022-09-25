%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The role of fiscal %%%
%%%%% policies in the env %%
%%%%% policy %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Sonja Dobkowitz
% building on Lint Barrage's code ReStud2019
% Version: August 2022
clear
%cd C:\Users\lintb\Desktop\COMET\ReStud
cd '/home/sonja/Documents/projects/subjective_BN/codding_model/own_basedOnFried/optimalPol_010922_revision'
%%
%M-File Outline%
%%%%%%%%%%%%%%%%
%General Notes
%Section 1: Select Fiscal Scenario
%....          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 1: Select Scenario        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 12;  % Direct optimization period time horizon: 2020-2080
         % one period = 5 years

lengthh = 5; % number of zears per period         
indic.util =0; % ==0 log utility, otherwise as  in Boppart

indic.Bop=0; % indicator ==1 then uses version as discussed in Boppart: 
                 % income effect stronger than substitution effect and
                 % thetaa > 1
indic.sep =0; % ==0 one joint market (in calibration very low fossil and green scientists to satisfy wage clearing 
              % ==1 3 separate markets 
              % ==2 if partial equbm; relative to joint market
              % ==3 energy market joint and non-energy market separate
indic.target =0; % ==1 if uses emission target
indic.noknow_spill =0; % ==0 then there are knowledge spillovers (benchmark model)
indic.sizeequ=0; %==1 then research sectors have same size => is there still a higher progressive tax when there are spillovers?
indic.spillovers =0; % ==1 then there are positive spillover effects of scientists within sectors! 
indic.noskill = 0; % == 1 if no skill calibration of model
indic.notaul=0; % Indicator of policy
                % ==0 benchmark redistribution via income tax scheme
                % ==1 no taul but env tax revenues redistributed via income tax scheme; 
                % ==2 gov consumes env. taxes; no taul  (No Income Tax scheme)
                % ==3 gov consumes env taxes; taul can be used 
                % ==4 taul is an option and env tax revenues are redistributed lump sum
                % ==5 lump sum trans, no taul
                % ==6 consolidated budget and taul adjusts (lambdaa fixed)
                % ==7 earmarking
indic.limit_LF =0; % ==1 then tauf is determined by meeting limit in each period
                   %  set by a planner who knows how economy works but each
                   %  period; not dynamic! (in optimal policy taking dynamics into account)
indic.taus =0; %==0 then no subsedy on green 
indic.xgrowth=0;
indic.extern=0; % extern==0 when uses no externality in utility
% but ensure no externality when target is used 
if indic.target==1
    indinc.extern=0;
end
indic.count_techgap=0; % if ==1 then uses technology gap as in Fried
indic.subs = 0; %==1 eppsy>1 (energy and neutral good are substitutes)
indic.PV = 1; % ==1 if continuation value is added to planners problem
indic.PVwork =0; %==0 then disutility of work is not in 
indic.emsbase=0; % ==0 then uses emission limits as calculated
indic.zero = 0;  % ==1 then version without growth
indic.GOV=0; % ==0 then no gov revenues
indic.taul0=0; %==0 then calibrated value for taul; ==1 then 0
indic.labshareequ =0; %==1 then equal capital shares in fossil green and non-energy sectors
indic

percon = 0;  % periods nonconstrained before 50\% constrained

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 2: Parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 1) sets direct parameters, 
%          2) calibrates model to indirect params.

if isfile(sprintf('params_0209_sep%d.mat', indic.sep))
    fprintf('loading parameter values sep%d', indic.sep)
    load(sprintf('params_0209_sep%d', indic.sep),...
        'params', 'Sparams', 'polCALIB', 'init201014', 'init201519', 'list', 'symms', 'Ems', 'Sall', 'x0LF', 'MOM', 'indexx')
else
    fprintf('calibrating model')
    indic.notaul=0; indic.limit_LF=0;indic.labshareequ =0; indic.sizeequ=0; indic.taul0=0; indic.GOV=0; indic.noknow_spill=0; indic.noskill=0; indic.xgrowth=0; 
    [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx, StatsEms]...
        =get_params_Base( T, indic, lengthh);
    save(sprintf('params_0209_sep%d', indic.sep))
end
if indic.spillovers==1
    params(list.params=='etaa')=1.2;
    Sparams=1.2;
end

%% 
%%%%% LIST for separate markets

     symms.sepchoice=symms.choice(list.choice~='S'&list.choice~='ws'& list.choice~='gammas');
     syms wsg wsn wsf gammasg gammasf gammasn real
    symms.sepchoice=[symms.sepchoice wsg wsn wsf gammasg gammasf gammasn];
    list.sepchoice=string(symms.sepchoice);
    
    %- allvars 
     symms.sepallvars=symms.allvars(list.allvars~='ws');
    symms.sepallvars=[symms.sepallvars wsg wsn wsf gammasg gammasf gammasn]; 
    list.sepallvars=string(symms.sepallvars);
    
    %- without growth
     symms.choice_xgrowth=symms.choice(list.choice~='S'&list.choice~='gammas'&list.choice~='sff'&list.choice~='sn'&list.choice~='sg'...
         &list.choice~='Af'&list.choice~='An'&list.choice~='Ag'&list.choice~='ws'&list.choice~='wsg'&list.choice~='wsn'&list.choice~='wsf'&list.choice~='gammasf'&list.choice~='gammasn'&list.choice~='gammasg');
     list.choice_xgrowth=string(symms.choice_xgrowth);
     symms.allvars_xgrowth=symms.allvars(list.allvars~='S'& list.allvars~='sff'&list.allvars~='sn'&list.allvars~='sg'...
         &list.allvars~='wsf'&list.allvars~='wsn'&list.allvars~='wsg'&list.allvars~='ws'&list.allvars~='gammasf'&list.allvars~='gammasn'&list.allvars~='gammasg');
     list.allvars_xgrowth=string(symms.allvars_xgrowth);

    % to save additional government variables
    syms GovCon Tls real;
    symms.addgov = [GovCon Tls];
    list.addgov = string(symms.addgov);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 3: BAU Simulation        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this section I simulate the economy starting from 2015-2019
% order of variables in LF_SIM as in list.allvars
indic.sizeequ=0;
indic.sep=2;
% indic.labshareequ=1;
% indic.noknow_spill=0; % ==1 then version without knowledge spillovers : should find different effect of tauf on growth and emissions?
POL=polCALIB; % tauf chosen in code; or updated below of limit-LF=0
% indic.tauff="BAU"; %=> also: to get BAU policy, run with tata=0, and indic.tauff=='BAU', GOV=='1', indic.limit_LF=='0'
                   % => to get LF run with BAU and tata=1
indic
for lablab =0:1
    indic.labshareequ=lablab;
for nknk =0
    indic.noknow_spill=nknk; 
for bb =["BAU", "SCC"]
    indic.tauff=bb;
for tata=0:1 % ==0 then uses calibrated taul
    indic.taul0=tata; %==1 then sets taul =0
for GG=0
    indic.GOV=GG; %==1 then lambdaa chosen to match government expenditures; ==0 then GOV=0

for cc=0
    indic.count_techgap=cc;
for ff=0:1
    indic.limit_LF=ff;
% choose environmental tax fixed
if indic.limit_LF==0
    if indic.tauff=='SCC'
        scc=185; % rff estimate per ton of carbon
        scc_giga22=185*1e9; % 185*1e9: price per gigaton in 2022 prices
        %=> in 2019 prices
        scc_giga19= scc_giga22/1.12; % 12% inflation from 2019 to 2022
        % => get in units of GDP over 2015-2019 horizon
        tauftilde= scc_giga19/(1e6*MOM.GDP1519MILLION); % in units of total GDP from 15 to 19 => price as percent of gdp
        POL(list.pol=='tauf')=tauftilde*Sparams.omegaa; % => tauf per unit of fossil in model
    elseif indic.tauff=='BAU'
        POL(list.pol=='tauf')=0;
    end
end

for ee=0
    indic.emsbase=ee;
    if indic.emsbase==1
        Ems_alt=x0LF(list.choice=='F')*0.7*ones(size(Ems))*Sparams.omegaa-Sparams.deltaa;
    else
        Ems_alt=Ems;
    end
% full model
for nsk=0
    indic.noskill=nsk;
    for xgr=0
        indic.xgrowth=xgr;
        % to save tauf
        TAUF=zeros(T,7); % 7= number of scenarios

    for nnt=[0,5]
        indic.notaul=nnt;

        if xgr==0
            if indic.count_techgap==0
                [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, POL, params, symms, x0LF, init201014, indexx, indic, Sall, Ems_alt);
                helper.LF_SIM=LF_SIM;
                [COMP]=solve_LF_VECT(T, list,  params,symms, init201519, helper, indic, Ems_alt, MOM);
            else
                iin=load('init_techgap.mat');
                [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, POL, params,   symms, x0LF, iin.initcount, indexx, indic, Sall, Ems_alt);
                 helper.LF_SIM=LF_SIM;
                [COMP]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic, Ems_alt, MOM);
            end
        else   
            if indic.count_techgap==0
               [COMP, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, POL, params,  symms, x0LF, init201014, indexx, indic, Sall, MOM, Ems_alt);
            else
                iin=load('init_techgap.mat');
               [COMP, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, POL, params,  symms, x0LF,iin.initcount, indexx, indic, Sall, MOM, Ems_alt);
            end
        end
        % umrechnung tauf in per ton of carbon in 2014-19 us dollars
        tauf=COMP(:,list.allvars=='tauf');
        tauf_CO2=tauf./Sparams.omegaa;
        tauf_perton2019 = tauf_CO2*(MOM.GDP1519MILLION*1e6)./(1e9); % denominator to go from gigaton to ton in 2019 prices
        TAUF(:,nnt+1)=tauf_perton2019*1.12; % to have it in 2022 prices
%         tauf_alt=tauf_CO2/x0LF(list.choice=='pf')*10.1554 *106.5850/100/0.07454;

    %-- save stuff
        if ~(indic.tauff=="BAU" && indic.limit_LF==0)
            save(sprintf('COMP_1409_taulZero%d_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_emlimit%d_Emsalt%d_countec%d_GovRev%d_etaa%.2f.mat',...
                indic.taul0, indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, indic.labshareequ, indic.sep, indic.notaul,indic.limit_LF,indic.emsbase, indic.count_techgap, indic.GOV,  params(list.params=='etaa')), 'COMP', 'tauf_perton2019')
            save(sprintf('TAUF_1409_taulZero%d_knspil%d_limit%d_EmsBase%d_xgr%d_nsk%d_labequ%d_countec%d_GovRev%d_sep%d',...
                indic.taul0, indic.noknow_spill, indic.limit_LF,indic.emsbase, indic.xgrowth, indic.noskill,indic.labshareequ, indic.count_techgap, indic.GOV, indic.sep), 'TAUF')
        else
            save(sprintf('BAU_1409_taulZero%d_spillovers%d_knspil%d_size_noskill%d_xgrowth%d_labequ%d_sep%d_notaul%d_countec%d_GovRev%d_etaa%.2f.mat',...
                indic.taul0, indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, indic.labshareequ,  indic.sep, indic.notaul, indic.count_techgap, indic.GOV,  params(list.params=='etaa')), 'COMP', 'tauf_perton2019')
        end                

       clearvars COMP pol TAUF tauf_perton2019 FVAL
    end
    end
end
end
end
end
end
end
end
end
end


%% plotts policy regimes comparison
etaa=params(list.params=='etaa');
weightext=0.01;
indic.GOV=0;
indic.sep=0;
plotts.regime=0;
indic.sizeequ=0;
indic

% choose sort of plots to be plotted
plotts.tauf_comp                = 0;
plotts.tauf_comp_Byregime       = 0;
plotts.tauf_compTaul            = 0;
plotts.tauf_compTaul_BYregime   = 0;
%- allocations with and without tauf 
plotts.perDif_notauf            = 0; %
plotts.perDif_notauf_compTaul   = 1;
plotts.tauf_notauf              = 0; % plots allocation with and without tauf in levels with and without taul and with and without equal labor share
plotts.compTauf_Lev             = 0; % compares allocation with tauf in model with and without taul in levels
plotts.compTauf_PER             = 0;
%- plots: effect of taul
plotts.LF_BAU                   = 0;
plotts.LF_BAU_PER               = 0;
plotts.LF_BAU_equlab            = 0;
plotts.LF_BAU_PER_equlab        = 0;

%- comparison policy regime
plotts.compRed                  = 0;
plotts.compTaul_Red             = 0;
plotts.compRed_TaulPer          = 0;

plotts.compRed_noGS             = 0;

for ee=0 % ==0 then uses benchmark emission limit
    indic.emsbase=ee;
        
for ll=0 % no emission limit : 
    indic.limit_LF=ll;
for nknk=0:1 % nowledge spillovers
    for xgr =0:1
        for nsk=0:1
    plotts.xgr = xgr; % main version to be used for plots
    plotts.nsk = nsk;
    indic.noknow_spill=nknk;

    plotts
    %%
    plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon);
    %plottsSP_tidiedUp(list, T-1, etaa, weightext,indic, params, Ems, plotts, percon); 
        end
    end
end
end
end

%% to be continued updating
%- version with counterfactual technology gap

% An0=init201014(list.init=='An0');
% Ag0=0.9*An0;
% Af0=Ag0/0.4; 
% initcount= eval(symms.init); % vector or counterfactual technology 
% iin=load('init_techgap.mat');
% [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, polCALIB, params, Sparams,  symms, x0LF, iin.initcount, indexx, indic, Sall);
%  helper.LF_SIM=LF_SIM;
% %    helper=load(sprintf('LF_BAU_spillovers%d.mat', indic.spillovers));
% %- initial gap 2015/19
% % An0=LF_SIM(list.sepallvars=='An', 1);
% % Ag0=LF_SIM(list.sepallvars=='Ag', 1);
% % Af0=LF_SIM(list.sepallvars=='Af', 1);
% % init1519count=eval(symms.init);
% [LF_BAU]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic);
% save(sprintf('BAU_countec_spillovers%d_knspil%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_SIM', 'Sparams')

%% Laissez faire
% 13/09: if interesed in BAU: need to set taul = 0.181 (above code does not
% give bau solution as tauf neq 0. 
% taus=0;
% tauf=0;
% taul=0;
% lambdaa=1; % placeholder, determined in comp eqbm
% pol=eval(symms.pol);
% 
% %%
% %  if indic.noskill==0
% % note: taul0 indicator irrelevant if taul in pol vector =0
% indic.limit_LF=0; 
% indic.noknow_spill=0;
% indic.sizeequ=0;
% 
% for gg=0
%     indic.GOV=gg;
% for nnt=[0,1,2,3,4,5,7]
%       indic.notaul=nnt;
%   for i=0:1
%       indic.noskill=i;
%       for xgr=0:1
%           indic.xgrowth=xgr;
% 
%           if indic.xgrowth==0
%             [LF_SIM, polLF, FVAL] =solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall, Ems);
% 
%             helper.LF_SIM=LF_SIM;
%             indic.xgrowth=0;
%             [LF_SIM]=solve_LF_VECT(T, list, params,symms, init201519, helper, indic,Ems);
%           else
%             [LF_SIM , pol, FVAL]= solve_LF_nows_xgrowth(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall, MOM, Ems);
%           end
%         save(sprintf('LF_SIM_spillover%d_knspil%d_noskill%d_xgr%d_sep%d_notaul%d_GOV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers, indic.noknow_spill, indic.noskill, indic.xgrowth, indic.sep,indic.notaul, indic.GOV, indic.sizeequ, params(list.params=='etaa')),'LF_SIM', 'Sparams');
%         clearvars LF_SIM helper
%       
%      end
%  end
% end
% end
%- version LF with counterfac tec gap
     
%     indic.xgrowth=0; % does not exist with exogenous growth
%     [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, iin.initcount, indexx, indic, Sall);
%     helper.LF_SIM=LF_SIM;
%     [LF_COUNTTec]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic, Ems_alt);
%     save(sprintf('LF_countec_spillovers%d_knspil%d_noskill%d_sep%d_notaul%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.notaul, indic.GOV, params(list.params=='etaa')), 'LF_SIM', 'Sparams')
    
%% end
%%% Check swf value in LF
disc=repmat(Sparams.betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;
hhel= load(sprintf('BAU_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')));
sswfbau=vec_discount*hhel.LF_BAU( :, list.sepallvars=='SWF');
hhblf = load(sprintf('LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.notaul,  params(list.params=='etaa')));
sswf=vec_discount*hhblf.LF_SIM( :, list.sepallvars=='SWF');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 4: Sociel Planner allocation                             %%%
% Timing: starting from 2020-2025  as initial period                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indic.sizeequ=0;
indic.noknow_spill=1;
indic.sep=0; % Has to be equal to loaded parameters!!!
for xgr=0:1
    indic.xgrowth=xgr;
    for ns=0:1
        indic.noskill=ns;
%             if ~isfile(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, params(list.params=='etaa')))
%                 indic.target=1;
%                 fprintf('solving Social planner solution with target, noskill%d', indic.noskill);
       for tar=0
            indic.target=tar;
            indic               
            SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems, MOM, percon);
       end            
    end
end

%%
if indic.count_techgap==1
    SP_solve(list, symms, params, Sparams, x0LF, initcount, init1519count, indexx, indic, T, Ems, MOM, percon);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 5: Solve for Optimal Allocation       %%%
% Timing: starting from 2020-2025 the gov. chooses      %%
% the optimal allocation                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indic.taus  = 0; % with ==0 no taus possible!
indic.sep =0;
indic.extern=0;
indic.GOV=0; % ==0 then no gov revenues
indic.sizeequ=0; 
indic.noknow_spill=1;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T

for tr =1
    indic.target=tr;
for xgr=1
    indic.xgrowth=xgr;
for nsk=1
    indic.noskill=nsk;
 for nnt=[5]
     indic.notaul=nnt;
     indic
 if indic.count_techgap==0
     OPT_solve_sep(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon);
 else
     OPT_solve_sep(list, symms, params, x0LF, init1519count, indexx, indic, T, Ems, MOM, percon);
 end
 end
end
end
end

%%
%-- extend optimality for count

indic.taus  = 0; % with ==0 no taus possible!
indic.sep =0;
indic.extern=0;
indic.GOV=0; % ==0 then no gov revenues
indic.sizeequ=0; 
indic.noknow_spill=0;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T
indic
count=30;% addiitonal periods
for tr =1
    indic.target=tr;
for xgr=[1,0]
    indic.xgrowth=xgr;
for nsk=[0, 1]
    indic.noskill=nsk;
 for nnt=1
     indic.notaul=nnt;
     indic
[symms, list, opt_all]= OPT_solve_sep_ExtT(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon, count);
 end
end
end
end
%%
%-- extend optimality for count
% for 0,3,4 => code with taul generated directly from non-taul alternative
indic.taus  = 0; % with ==0 no taus possible!
indic.sep =0;
indic.extern=0;
indic.GOV=0; % ==0 then no gov revenues
indic.sizeequ=0; 
indic.noknow_spill=0;
indic.limit_LF=0; % no need to test this
indic.testT =0; % do not test value of T
indic
count=30;% addiitonal periods
for tr =1
    indic.target=tr;
for xgr=1
    indic.xgrowth=xgr;
for nsk=0:1
    indic.noskill=nsk;
 for nnt=5
     indic.notaul=nnt;
     indic
[symms, list, opt_all]= OPT_solve_sep_ExtT_direct(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon, count);
 end
end
end
end
%% TO BE UPDATED FOR NEW STRUCTURE OF TAX! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: Competitive equi 
%%%      counterfactual policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tf=4 %==4 (set optimal pol without no spil in benchmark model); ==2 then uses policy as in helper.opt_all; only benchmark taul, tauf=0 in other models
indic.tauf=tf; % ==1 uses version with optimal tauf but taul=0; ==0 uses version with tauf optimal but taul =0
indic.notaul=0;
indic.PV=1;
indic.noknow_spill=0; % counterfactuals so far only without knowledge spillovers
T=12;

for xgr=1
    indic.xgrowth=xgr;
    for nsk=0:1
        indic.noskill=nsk;
% load benchmark policy
if tf>=2 % read in benchmark model wrt skill xgr
    if tf~=4
        helper=load(sprintf('OPT_target_1008_spillover%d_knspil%d_taus0_noskill0_notaul%d_sep%d_xgrowth0_PV%d_etaa%.2f.mat',indic.spillovers, indic.noknow_spill, indic.notaul, indic.sep, indic.PV, Sparams.etaa));
    elseif tf==4
        helper=load(sprintf('OPT_target_1008_spillover%d_knspil1_taus0_noskill0_notaul%d_sep%d_xgrowth0_PV%d_etaa%.2f.mat',indic.spillovers, indic.notaul, indic.sep, indic.PV, Sparams.etaa));
    end
elseif tf<2 % load in same model wsrt skill xgr
     helper=load(sprintf('OPT_target_1008_spillover%d_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noknow_spill, indic.noskill, indic.notaul, indic.sep, indic.xgrowth, indic.PV, Sparams.etaa));
end
         if  indic.tauf==1 || indic.noskill==1 % get better starting values!
             if indic.tauf<2
                 T=11;
             end
             LF_SIM=helper.opt_all;       
                for ll=list.choice
                    x0LF(list.choice==ll)=LF_SIM(1, list.allvars==ll);
                end
                 if indic.tauf>=2
                     poll= [LF_SIM(:,list.sepallvars=='taul'), LF_SIM(:,list.sepallvars=='taus'),LF_SIM(:,list.sepallvars=='tauf'), LF_SIM(:,list.sepallvars=='lambdaa')];
                 elseif indic.tauf==0
                     poll= [LF_SIM(:,list.sepallvars=='taul'), LF_SIM(:,list.sepallvars=='taus'),zeros(size(LF_SIM(:,list.sepallvars=='tauf'))), LF_SIM(:,list.sepallvars=='lambdaa')];
                 elseif indic.tauf==1
                     poll= [zeros(size(LF_SIM(:,list.sepallvars=='taul'))), LF_SIM(:,list.sepallvars=='taus'),LF_SIM(:,list.sepallvars=='tauf'), LF_SIM(:,list.sepallvars=='lambdaa')];
                     
                 end
             if indic.xgrowth==0
                    
                 [LF_SIM, pol, FVAL] = solve_LF_nows(T, list,poll, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
             else
                 [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, poll, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
             end
             helper.LF_SIM=LF_SIM'; % but not correct policy!

             if indic.tauf<2
                helper.LF_SIM=[helper.LF_SIM;helper.LF_SIM(end,:)];
             end
                 
             for pp=list.pol(list.pol~='lambdaa')
                 helper.LF_SIM(:, list.allvars==pp)= helper.opt_all(:, list.allvars==pp);
             end
         else
           helper.LF_SIM=helper.opt_all;

         end
         T=12; 
        [LF_COUNT]=compequ(T, list, params, init201519, symms, helper.LF_SIM,indic);
        % helper.opt_all: as initial values and to deduce policy 
        save(sprintf('COMPEquN_SIM_taufopt%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.tauf,  indic.spillovers, indic.notaul, indic.noskill, indic.sep, indic.xgrowth, indic.PV, params(list.params=='etaa')),'LF_COUNT', 'Sparams');
    end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: PLOTS       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaa=params(list.params=='etaa');
weightext=0.01;
indic

% choose sort of plots to be plotted
plotts.regime_gov=  0; % = equals policy version to be plotted

plotts.table=       0;
plotts.cev  =       0; 
plotts.analyta =    0;
plotts.limit=       0; %==1 if plots emission target
plotts.robust=      0;
plotts.countcomp=   0;
plotts.countcomp2=  0;
plotts.countcomp3=  0;
plotts.extern=      0;
plotts.compEff_mod_dev1         = 0;
plotts.count_taul_nsk_LF        = 0;
plotts.count_taul_xgr_LF        = 0;
plotts.count_taul_xgr_lev       = 0;
plotts.count_tauflev            = 0; % counterfactual with only tauf in laissez faire
plotts.count_taullev            = 0; % counterfactual with only taul in laissez faire
plotts.count_tauflev_Ben        = 0; % laissez faire, only optimal tauf and benchmark policy
plotts.count_tauflev_Ben_noLF   = 0;
plotts.compnsk_xgr              = 0;
plotts.compnsk_xgr1             = 0;

plotts.compnsk_xgr_dev          = 0;
plotts.compnsk_xgr_dev1         = 0;
plotts.count_modlev             = 0; 

plotts.count_modlev_eff         = 0;
plotts.single_pol               = 1;
plotts.singov                   = 0;

plotts.notaul                   = 0; % policy comparisons; this one needs to be switched on to get complete table
plotts.bau                      = 0; % do plot bau comparison
plotts.lf                       = 0; % comparison to laissez faire allocation 

plotts.comptarg                 = 0; % comparison with and without target
plotts.compeff                  = 0; % efficient versus optimal benchmark and non-benchmark
plotts.compeff3                 = 0; % sp versus optimal benchmark
plotts.comp_LFOPT               = 0; % laissez faire and optimal with and without taul
plotts.compeff1=    0; %1; only social planner
plotts.compeff2=    0; %1; efficient and non benchmark
plotts.comp_OPT=    0; % laissez faire and optimal with and without taul
plotts.comp_OPT_NK= 1; % laissez faire and optimal with and without taul
plotts.comp_Bench_CountNK =0; % policy from model without knowledge spillovers in benchmark model
plotts.per_BAUt0 =  0;
plotts.per_effopt0= 0;
plotts.per_effoptd= 0;
plotts.per_baud =   0;
plotts.per_LFd  =   0; % dynamic lf as benchmark
plotts.per_LFd_nt=  0; % dynamic lf as benchmark plus no income tax
plotts.per_LFd_ne_nt=0; % dynamic lf as benchmark plus no income tax

plotts.per_LFt0  =  0; % 2020  lf as benchmark
plotts.per_optd =   0;

plotts.tauf_comp=0;
plotts.compREd=0;
    
for xgr =0
    for nsk=0
        for se=0
plotts.xgr = xgr; % main version to be used for plots
plotts.nsk = nsk;
plotts.sizeequ =se; % important for comparison of 
plotts.GOV =0;
indic.noknow_spill=0; % in the benchmark allocation there are kn spillovers

plotts
%%
%     plottsSP_PolRegimes(list, T, etaa, weightext,indic, params, Ems, plotts, percon);
plottsSP_tidiedUp(list, T-1, etaa, weightext,indic, params, Ems, plotts, percon); 
        end
    end
end
%%
for gg=0:1
    indic.xgrowth=gg;
for ns=0:1
    indic.noskill=ns;
    plottsSP(list, T, etaa, weightext,indic, params, Ems, plotts, percon); 
    % T-1 as time period as last period is dropped from regency 
end
end

%% tables
% constructed in plotsSP file

tt=load('Table_SWF_July22_sep1_noskill0_etaa0.79_xgrowth0_extern0.mat');
addpath('tools')
table2latex(tt.TableSWF_PV)
%% Analytical model files
% solves analytical model and compares policy regimes
solve_easy;

%% Simulate economy after last optimisation period using competitive economy files
% under the assumption of constant policy
% GOAL: check if emission target gets violated

indic.notaul =3;

indic.noskill=0;
indic.taus   =0;
indic.sep    =1;
indic.spillovers = 0;

for xgr =0
    indic.xgrowth=xgr;
    res_opt=load(sprintf('OPT_target_spillover%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth, params(list.params=='etaa')), 'opt_all', 'Sparams');
    % simulate in laissez faire for peroiods after 2075 with assumption of
    % constant policy (after 2080 because laissez faire code starts from first period)

    %read in policy
    taul = res_opt.opt_all(T, list.sepallvars=='taul');
    tauf = res_opt.opt_all(T, list.sepallvars=='tauf');
    taus = res_opt.opt_all(T, list.sepallvars=='taus');
    lambdaa = res_opt.opt_all(T, list.sepallvars=='lambdaa');
    polLFcon = eval(symms.pol);
    %- initial 2080 technology levels
    % => first simulation period is 2085
    Af0 = res_opt.opt_all(T, list.sepallvars=='Af');
    Ag0 = res_opt.opt_all(T, list.sepallvars=='Ag');
    An0 = res_opt.opt_all(T, list.sepallvars=='An');
    initcon=eval(symms.init);

    %- initial values as in x0LF; indicator: list.choice
    x0LFcon=x0LF;
    for ii =list.choice
        if ismember(ii, ["ws" "gammas"])
            x0LFcon(list.choice==ii) = res_opt.opt_all(T, list.sepallvars==sprintf('%sg', ii ));
        else
            x0LFcon(list.choice==ii) = res_opt.opt_all(T, list.sepallvars==ii);
        end
    end
    %- choose periods for which to simulate
    Tcon = 12; 

    if indic.xgrowth==0
        [LF_SIM, polLF, FVAL] =solve_LF_nows_continued(Tcon, list, polLFcon, params, Sparams,  symms, x0LFcon, initcon, indexx, indic, Sall);
        helper.LF_SIM=LF_SIM;
        indic.xgrowth=0;
        [LF_SIM]=solve_LF_VECT(T, list, params,symms, initcon, helper, indic);
    else
        [LF_SIM, polLF, FVAL, indexx] = solve_LF_nows_xgrowth_continued(Tcon, list, polLFcon, params, Sparams,  symms, x0LFcon, initcon, indexx, indic, Sall);
         LF_SIM = LF_SIM';   
    end

    save(sprintf('LF_CON_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, Sparams.etaa), 'Tcon', 'LF_SIM', 'Sparams', 'polLFcon');
% evaluate how economy evolves
Fcon= LF_SIM(:,list.sepallvars=='F');
Gcon= LF_SIM(:,list.sepallvars=='G');
GFFcon = Gcon./Fcon;
Ftargetcon = (Ems(end)+Sparams.deltaa)/Sparams.omegaa;
    if sum(Fcon >Ftargetcon)
        fprintf('The emission target gets violated in at least one period after intervention; exog growth %d', indic.xgrowth);
    end    
end

concon= load(sprintf('LF_CON_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, Sparams.etaa),'LF_SIM');
%=> given the model's assumptions there needs to be continued intervention
%   to meet the emission limit
%%
for BN=1
    indic.BN=BN;
    for inn=0:1
        indic.ineq=inn;
        indic
        plottsSP(list, T, etaa, weightext,indic, params, Ems, plotts);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find good initial starting point for target opt %%
%%% as a function of taul
taul= 0.5;
pf= 1;
exx = polExp(pf, params, polCALIB, list, taul, T, Ems, indexx);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Symbolic approach to solve Ramsey problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indic.target=1;
indic.noskill=0;
% indic.sep=0;
indic.target=1; 
%- with etaa=1
% params(list.params=='etaa')=1;
% Sparams.etaa=params(list.params=='etaa');

%1) get objective function 
if indic.sep==1
    [OB_RAM, list, symms, Ftarget]= model_ram_sep( list, params, T, init201519, indic, Ems, symms);
else
    [OB_RAM, list, symms, Ftarget]= model_ram( list, params, T, init201519, indic, Ems, symms);
end
%- x is a symbolic vector of choice variables! 

%2) take derivatives and write resulting equations as function
if indic.target==1
    [indexx, model, list]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_target_1905_KTS', list, indic, indexx);
else
    if indic.sep==0
        [indexx, model]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_1905_KTS', list, indic, indexx);
    else
        [indexx, model]=symmodel_eq_sep(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_sep_1905', list, indic, indexx);
    end
end

%3) solve model using fsolve
if indic.sep==1
    RAM = solve_sym_sep(symms, list, Ftarget, indic);
else
    RAM = solve_sym(symms, list, Ftarget, indic);
end
