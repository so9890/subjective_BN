%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%COMET Model M-File%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Sonja Dobkowitz
% building on Lint Barrage's code ReStud2019
% Version: August 2022
clear
%cd C:\Users\lintb\Desktop\COMET\ReStud
cd '/home/sonja/Documents/DocumentsSonja/projects/subjective_BN/codding_model/own_basedOnFried/optimalPol_190722_tidiedUp'
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
indic.sep =1; %==1 is the benchmark; when finalising should be dropped
indic.target =0; % ==1 if uses emission target
indic.noknow_spill =0; % ==0 then there are knowledge spillovers (benchmark model)
indic.sameSize =1; % ==0 then non-energy sector bigger (benchmark model)

indic.spillovers =0; % ==1 then there are positive spillover effects of scientists within sectors! 
indic.noskill = 0; % == 1 if no skill calibration of model
indic.notaul=0; % Indicator of policy
                % ==0 benchmark redistribution via income tax scheme
                % ==1 no taul but env tax revenues redistributed via income tax scheme; 
                % ==2 gov consumes env. taxes; taul cannot be used (No Income Tax scheme)
                % ==3 gov consumes env taxes; taul can be used 
                % ==4 taul is an option and env tax revenues are redistributed lump sum
                % ==5 lump sum trans, no taul
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
indic

percon = 0;  % periods nonconstrained before 50\% constrained

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 2: Parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 1) sets direct parameters, 
%          2) calibrates model to indirect params.
if isfile(sprintf('params_0308.mat'))
    fprintf('loading parameter values')
    load(sprintf('params_0308.mat'),...
        'params', 'Sparams', 'polCALIB', 'init201014', 'init201519', 'list', 'symms', 'Ems', 'Sall', 'x0LF', 'MOM', 'indexx')
else
    fprintf('calibrating model')
    [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx]=get_params_Base( T, indic, lengthh);
    save(sprintf('params_0308'))
end
if indic.spillovers==1
    params(list.params=='etaa')=1.2;
    Sparams=1.2;
end

%% 
%%%%% LIST for separate markets
%     symms.sepchoice=symms.choice(list.choice~='S'&list.choice~='ws'& list.choice~='gammas');
%     syms wsg wsn wsf gammasg gammasf gammasn real
    symms.sepchoice=symms.choice; %[symms.sepchoice wsg wsn wsf gammasg gammasf gammasn];
    list.sepchoice=string(symms.sepchoice);
    
    %- allvars 
%     symms.sepallvars=symms.allvars(list.allvars~='ws');
    symms.sepallvars=symms.allvars; %[symms.sepallvars wsg wsn wsf gammasg gammasf gammasn]; 
    list.sepallvars=string(symms.sepallvars);
    
    %- without growth
    symms.choice_xgrowth=symms.choice(list.choice~='sff'&list.choice~='sn'&list.choice~='sg'&list.choice~='Af'&list.choice~='An'&list.choice~='Ag'&list.choice~='wsg'&list.choice~='wsn'&list.choice~='wsf'&list.choice~='gammasf'&list.choice~='gammasn'&list.choice~='gammasg');
    list.choice_xgrowth=string(symms.choice_xgrowth);
    symms.allvars_xgrowth=symms.allvars(list.allvars~='sff'&list.allvars~='sn'&list.allvars~='sg'&list.allvars~='wsf'&list.allvars~='wsn'&list.allvars~='wsg'&list.allvars~='gammasf'&list.allvars~='gammasn'&list.allvars~='gammasg');
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
 
% full model
for nnt=[0,1,4,5]
    indic.notaul=nnt;
for i=0:1
    indic.noskill=i;
%     if ~isfile(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')))
        [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, polCALIB, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
        helper.LF_SIM=LF_SIM;
        [LF_BAU]=solve_LF_VECT(T, list,  params,symms, init201519, helper, indic);
       
        save(sprintf('BAU_spillovers%d_knspil%d_size_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_BAU', 'Sparams')
        clearvars LF_SIM pol FVAL
%     end
     fprintf('LF_BAU no skill %d exists', indic.noskill);


%%
%- version without growth
[LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, polCALIB, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
save(sprintf('BAU_xgrowth_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_SIM', 'Sparams')

%- version with counterfactual technology gap

% An0=init201014(list.init=='An0');
% Ag0=0.9*An0;
% Af0=Ag0/0.4; 
% initcount= eval(symms.init); % vector or counterfactual technology 
iin=load('init_techgap.mat');
[LF_SIM, pol, FVAL] = solve_LF_nows(T, list, polCALIB, params, Sparams,  symms, x0LF, iin.initcount, indexx, indic, Sall);
 helper.LF_SIM=LF_SIM;
%    helper=load(sprintf('LF_BAU_spillovers%d.mat', indic.spillovers));
%- initial gap 2015/19
% An0=LF_SIM(list.sepallvars=='An', 1);
% Ag0=LF_SIM(list.sepallvars=='Ag', 1);
% Af0=LF_SIM(list.sepallvars=='Af', 1);
% init1519count=eval(symms.init);
[LF_BAU]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic);
save(sprintf('BAU_countec_spillovers%d_knspil%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_SIM', 'Sparams')
end
end
%% Laissez faire
taus=0;
tauf=0;
taul=0;
lambdaa=1; % placeholder, determined in comp eqbm
pol=eval(symms.pol);

%%
  %  if indic.noskill==0
for nnt=[0,1,4,5]
      indic.notaul=nnt;
  for i=0:1
      indic.noskill=i;

        [LF_SIM, polLF, FVAL] =solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);

        helper.LF_SIM=LF_SIM;
        indic.xgrowth=0;
        [LF_SIM]=solve_LF_VECT(T, list, params,symms, init201519, helper, indic);
        save(sprintf('LF_SIM_NOTARGET_spillover%d_knspil%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep,indic.notaul, params(list.params=='etaa')),'LF_SIM', 'Sparams');
        clearvars LF_SIM helper

        %- exogenous growth
%          indic.xgrowth=1;
%          [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
%          % helper.LF_SIM=LF_SIM;
%          save(sprintf('LF_xgrowth_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_SIM', 'Sparams')

    %- version LF with counterfac tec gap
    indic.xgrowth=0; % does not exist with exogenous growth
    [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, iin.initcount, indexx, indic, Sall);
    helper.LF_SIM=LF_SIM;
    [LF_COUNTTec]=solve_LF_VECT(T, list, params,symms, iin.init1519count, helper, indic);
    save(sprintf('LF_countec_spillovers%d_knspil%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.notaul, params(list.params=='etaa')), 'LF_SIM', 'Sparams')
  end
end
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

for xgr=0
    indic.xgrowth=xgr;
    for ns=0:1
        indic.noskill=ns;
%             if ~isfile(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, params(list.params=='etaa')))
%                 indic.target=1;
%                 fprintf('solving Social planner solution with target, noskill%d', indic.noskill);
       for tar=1
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
indic.sep =1;
indic.extern=0;

for tr =0:1
    indic.target=tr;
for xgr=1
    indic.xgrowth=xgr;
for nsk=0:1
    indic.noskill=nsk;
 for nnt=[3]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: PLOTS       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaa=params(list.params=='etaa');
weightext=0.01;
indic

% choose sort of plots to be plotted
plotts.regime_gov=  3; % = equals policy version to be plotted

plotts.table=       0;
plotts.cev  =       0; 
plotts.analyta =    0;
plotts.limit=       0; %==1 if plots emission target
plotts.robust=      0;
plotts.countcomp=   0;
plotts.countcomp2=  0;
plotts.countcomp3=  0;
plotts.extern=      0;
plotts.compEff_mod_dev1=0;
 plotts.count_taul_nsk_LF=0;
plotts.count_taul_xgr_LF =0;
plotts.count_taul_xgr_lev =0;
plotts.count_tauflev =0; % counterfactual with only tauf in laissez faire
plotts.count_taullev =0; % counterfactual with only taul in laissez faire
plotts.count_tauflev_Ben =0; % laissez faire, only optimal tauf and benchmark policy
plotts.count_tauflev_Ben_noLF=1;
plotts.compnsk_xgr = 0;
plotts.compnsk_xgr1= 0;

plotts.compnsk_xgr_dev= 0;
plotts.compnsk_xgr_dev1 =0;
plotts.count_modlev= 0; 

plotts.count_modlev_eff= 0;
plotts.single_pol=  0;
plotts.singov=      0;

plotts.notaul=      0; % policy comparisons; this one needs to be switched on to get complete table
plotts.bau=         0; % do plot bau comparison
plotts.lf=          0; % comparison to laissez faire allocation 

plotts.comptarg=    0; % comparison with and without target
plotts.compeff=     0; % efficient versus optimal benchmark and non-benchmark
plotts.compeff3=    0; % sp versus optimal benchmark
plotts.comp_LFOPT=  0; % laissez faire and optimal with and without taul
plotts.compeff1=    0; %1; only social planner
plotts.compeff2=    0; %1; efficient and non benchmark
plotts.comp_OPT=    0; % laissez faire and optimal with and without taul
plotts.comp_OPT_NK= 0; % laissez faire and optimal with and without taul
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

for xgr =0
    for nsk=0
plotts.xgr = xgr; % main version to be used for plots
plotts.nsk = nsk;
plotts
%%
plottsSP_tidiedUp(list, T-1, etaa, weightext,indic, params, Ems, plotts, percon); 
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
