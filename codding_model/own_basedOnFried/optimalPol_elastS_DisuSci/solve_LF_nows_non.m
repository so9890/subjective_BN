function [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_non(T, list, pol, params, Sparams,  symms, x0LF, init, indexx, indic, Sall)
% simulate economy under laissez faire
indic.xgrowth = 0;
indic.ineq = 0;
indic.util = 0;

if indic.count_techgap==1
    iin=load('init_techgap.mat');
    init = iin.initcount; % 201014 values
end

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%% first run version with only one skill type
syms F pf Lg Af Ag sff sg wsf wsg gammalh gammasg gammasf real
symms.choice_small = [F pf Lg Af Ag sff sg wsf wsg gammalh gammasg gammasf];
list.choice_small=string(symms.choice_small);

F = x0LF(list.choice=='F');
pf = x0LF(list.choice=='pf');
Af = x0LF(list.choice=='Af');
Ag = x0LF(list.choice=='Ag');
sff = x0LF(list.choice=='sff');
sg = x0LF(list.choice=='sg');
wsf = x0LF(list.choice=='ws');
wsg = x0LF(list.choice=='ws');
alphag=Sparams.alphag;
Lg = x0LF(list.choice=='G')./(Ag*(x0LF(list.choice=='pg').*alphag).^(alphag./(1-alphag)));

gammalh = x0LF(list.choice=='gammalh');
gammasg = x0LF(list.choice=='gammas');
gammasf = x0LF(list.choice=='gammas');
x0=eval(symms.choice_small);

indexxLF.lab = boolean(zeros(size(list.choice_small)));
indexxLF.exp = boolean(zeros(size(list.choice_small)));
indexxLF.sqr = boolean(zeros(size(list.choice_small)));
indexxLF.oneab = boolean(zeros(size(list.choice_small)));

indexxLF.lab(  list.choice_small=='sff'  | list.choice_small=='sg')=1;
indexxLF.exp(list.choice_small~='sff'&list.choice_small~='sg'...
    &list.choice_small~='gammalh'&list.choice_small~='gammasg'& list.choice_small~='gammasf')=1;
indexxLF.sqr(list.choice_small=='gammalh'| list.choice_small=='gammasg'|list.choice_small=='gammasf' )=1;

indexx('LF_noneutral_sep_noskill')=indexxLF;

% transform data
guess_trans=trans_guess(indexx('LF_noneutral_sep_noskill'), x0, params, list.params);

% test
f=laissez_faire_nows_sep_non_noskillSmall(guess_trans, params, list, pol, init, indic);

modFF = @(x)laissez_faire_nows_sep_non_noskillSmall(x, params, list, pol, init, indic);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e4, 'MaxIter', 5e5,  'Algorithm', 'levenberg-marquardt', 'Display', 'Iter');%, );%, );%, );
[sol2, fval, exitf] = fsolve(modFF, guess_trans, options);
[sol2, fval, exitf] = fsolve(modFF, sol2, options);

options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e6, 'MaxIter', 5e5, 'Display', 'Iter');%, );%, );%, );
[sol2, fval, exitf] = fsolve(modFF, guess_trans, options);
[sol2, fval, exitf] = fsolve(modFF, sol2, options);
 
f=laissez_faire_nows_sep_non_noskillSmall(sol2, params, list, pol, init, indic);
save('noskill_noneutral_growth_orgtechgap_1706', 'sol2')

% transform
LF=trans_allo_out(indexx('LF_noneutral_sep_noskill'), sol2, params, list.params, indic);
cell_par=arrayfun(@char, symms.choice_small, 'uniform', 0);
SLF=cell2struct(num2cell(LF), cell_par, 2);
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare lists if using separate markets for scientists
     %- new initial point: separate markets
     
    
    x0=[x0LF(list.choice~='hln'  & list.choice~='sn' & list.choice~='pn'& list.choice~='pee'&list.choice~='pg'& ...
        list.choice~='An' & list.choice~='hhn'&list.choice~='ws' & list.choice~='S' & list.choice~='gammas'), ...
        x0LF(list.choice=='ws'), x0LF(list.choice=='ws'), 0,0];
   
    list.choice=list.choice_non;
    symms.choice=symms.choice_non;
    %- new index for tranformation
        
    indexxLF.lab = boolean(zeros(size(list.choice)));
    indexxLF.exp = boolean(zeros(size(list.choice)));
    indexxLF.sqr = boolean(zeros(size(list.choice)));
    indexxLF.oneab = boolean(zeros(size(list.choice)));

    indexxLF.lab( list.choice=='hh'|list.choice=='hl' | list.choice=='sff'  | list.choice=='sg')=1;
    indexxLF.exp(list.choice~='hl'&list.choice~='hh'&list.choice~='sff'&list.choice~='sg'...
        &list.choice~='gammalh'&list.choice~='gammall'& list.choice~='gammasg'& list.choice~='gammasf')=1;
    indexxLF.sqr(list.choice=='gammalh'|list.choice=='gammall'| list.choice=='gammasg'|list.choice=='gammasf' )=1;
    
    indexx('LF_nonsep')=indexxLF;
    
% initialise stuff
%-- to save results
% symms.allvars=[symms.allvars, gammasn, gammasg, gammasf];
% list.allvars=string(symms.allvars);
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);

%-- initialise values
laggs   = init; % (init should refer to 2010-2014 period)
t       = 1; % number of periods: t=1: 2015-2019 => does include base year period (in matrix on first row) but dont save!

%-- change list if noskill ==1 version
if indic.noskill==1
    x0=x0(list.choice~='hl'&list.choice~='hh'& list.choice~='hlf'&list.choice~='hlg'& list.choice~='gammall'&...
                 list.choice~='hhg'& list.choice~='hhf'& list.choice~='wl'& list.choice~='wh');
    syms h w Lf Lg real
    symms.choice= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'), h, w, Lf, Lg];
    list.choice=string(symms.choice);
    x0=[x0,Sall.hh, Sall.wh, Sall.Lf,Sall.hh-Sall.Lf]; % order has to match how variables are added to symms.choice!
   
    %- update initial values
    x0(list.choice=='h')= ((1-pol(list.pol=='taul'))/Sparams.chii)^(1/(1+Sparams.sigmaa));   
    x0(list.choice=='lambdaa')=(x0(list.choice=='w').*x0(list.choice=='h')+pol(list.pol=='tauf').*x0(list.choice=='pf').*x0(list.choice=='F'))./(x0(list.choice=='w').*x0(list.choice=='h')).^(1-pol(list.pol=='taul')); 
    x0(list.choice=='C')=x0(list.choice=='lambdaa').*(x0(list.choice=='w').*x0(list.choice=='h')).^(1-pol(list.pol=='taul')); 
    x0(list.choice=='Lf')=x0(list.choice=='F')./(((1-pol(list.pol=='tauf')).*Sparams.alphaf.*x0(list.choice=='pf')).^(Sparams.alphaf./(1-Sparams.alphaf)).*x0(list.choice=='Af'));
    x0(list.choice=='Lg')=x0(list.choice=='h')-x0(list.choice=='Lf');
    
    indexxLFsep.lab = boolean(zeros(size(list.choice)));
        indexxLFsep.exp = boolean(zeros(size(list.choice)));
        indexxLFsep.sqr = boolean(zeros(size(list.choice)));
        indexxLFsep.oneab = boolean(zeros(size(list.choice)));

        indexxLFsep.lab( list.choice=='h'| list.choice=='sff'  | list.choice=='sg' )=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='sff'&list.choice~='sg'...
            &list.choice~='gammalh'& list.choice~='gammasg'& list.choice~='gammasf')=1;
        indexxLFsep.sqr(list.choice=='gammalh'| list.choice=='gammasg'| list.choice=='gammasf' )=1;

        indexx('LF_noskill_sep_non')=indexxLFsep;
end


while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 
    if indic.noskill==0
        guess_trans=trans_guess(indexx('LF_nonsep'), x0, params, list.params);
    else
        guess_trans=trans_guess(indexx(sprintf('LF_noskill_sep_non')), x0, params, list.params);
    end
    %% test
    f=laissez_faire_nows_sep_non(guess_trans, params, list, pol, laggs, indic);
    
    %% - solving model
     lb=[];
     ub=[];
     objf=@(x)objectiveCALIBSCI(x);
     
     constrf = @(x)laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic);
     
options = optimset('algorithm','active-set','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
%[c,fval]=constrf(sol3);
% 
%- other solvers
    modFF = @(x)laissez_faire_nows_sep_non(x, params, list, pol, laggs, indic);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e4, 'MaxIter', 5e5,  'Algorithm', 'levenberg-marquardt', 'Display', 'Iter');%, );%, );%, );
    [sol2, fval, exitf] = fsolve(modFF, guess_trans, options);

    % pass to standard algorithm
%      options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,);%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, x1, options);
% 
     options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
     [sol3, fval, exitf] = fsolve(modFF, sol2, options);

    %- transform results to bounded variables
    if indic.noskill==0
        if indic.sep==0
            LF=trans_allo_out(indexx('LF'), sol3, params, list.params, indic);
        else
            if indic.BN==0
                LF=trans_allo_out(indexx('LF_sep'), sol3, params, list.params, indic);
            else
                LF=trans_allo_out(indexx('LF_sep_BN'), sol3, params, list.params, indic);
            end
        end
     else
        LF=trans_allo_out(indexx(sprintf('LF_noskill_sep%d', indic.sep)), sol3, params, list.params, indic);
    end
    %% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    if t>1
        if indic.sep==0
            LF_SIM(:,t-1)=aux_solutionLF(Sparams, SLF, pol, laggs, list, symms, indexx, params, indic);
        else
            LF_SIM(:,t-1)=aux_solutionLF_sep(Sparams, SLF, pol, laggs, list, symms, indexx, params, indic);
        end
        FVAL(t-1)=max(abs(fval));
    end
    %% - update for next round
    x0 = LF; % initial guess
        Af0= SLF.Af; % today's technology becomes tomorrow's lagged technology
        Ag0= SLF.Ag; 
        An0= SLF.An; 
    laggs=eval(symms.init);
    t=t+1;
end
end