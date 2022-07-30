function [LF_SIM, pol, FVAL, indexx] = solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init, indexx, indic, Sall)
% simulate economy under laissez faire
indic.xgrowth=0;
indic.ineq=0;
% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare lists if using separate markets for scientists
if indic.sep==1
    %- new initial point
    x0=[x0LF(list.choice~='ws' & list.choice~='S' & list.choice~='gammas'), ...
        x0LF(list.choice=='ws'), x0LF(list.choice=='ws'), x0LF(list.choice=='ws'), 0,0,0];
    %- new set of choice variables
    list.choice=list.sepchoice;
    symms.choice=symms.sepchoice;
    list.allvars=list.sepallvars;
    symms.allvars=symms.sepallvars;
    
    %- new index for tranformation
        
    indexxLF.lab = boolean(zeros(size(list.choice)));
    indexxLF.exp = boolean(zeros(size(list.choice)));
    indexxLF.sqr = boolean(zeros(size(list.choice)));
    indexxLF.oneab = boolean(zeros(size(list.choice)));

    indexxLF.lab( list.choice=='hh'|list.choice=='hl' | list.choice=='sff'  | list.choice=='sg' | list.choice=='sn')=1;
    indexxLF.exp(list.choice~='hl'&list.choice~='hh'&list.choice~='sff'&list.choice~='sg'& list.choice~='sn'...
        &list.choice~='gammalh'&list.choice~='gammall'& list.choice~='gammasg'& list.choice~='gammasf'& list.choice~='gammasn')=1;
    indexxLF.sqr(list.choice=='gammalh'|list.choice=='gammall'| list.choice=='gammasg'| list.choice=='gammasn'| list.choice=='gammasf' )=1;
    
    indexx('LF_sep')=indexxLF;
    
    %-version with BN
   indexxLF.BN = boolean(zeros(size(list.choice)));
   indexxLF.exp = boolean(zeros(size(list.choice)));

    indexxLF.exp(list.choice~='C' & list.choice~='hl'&list.choice~='hh'&list.choice~='sff'&list.choice~='sg'& list.choice~='sn'...
        &list.choice~='gammalh'&list.choice~='gammall'& list.choice~='gammasg'& list.choice~='gammasf'& list.choice~='gammasn')=1;

    indexxLF.BN(list.choice=='C')=1;
    indexx('LF_sep_BN')=indexxLF;

else
    x0      = x0LF;
end

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
    x0=x0(list.choice~='hl'&list.choice~='hh'& list.choice~='hln'& list.choice~='hlf'&list.choice~='hlg'& list.choice~='gammall'&...
                list.choice~='hhn'& list.choice~='hhg'& list.choice~='hhf'& list.choice~='wl'& list.choice~='wh');
    syms h w Lf Lg Ln real
    symms.choice= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'), h, w, Lf, Lg, Ln];
    list.choice=string(symms.choice);
    x0=[x0,Sall.hh, Sall.wh, Sall.Lf,Sall.Lg,Sall.Ln]; % order has to match how variables are added to symms.choice!

    if indic.sep==0
        indexxLFsep.lab = boolean(zeros(size(list.choice)));
        indexxLFsep.exp = boolean(zeros(size(list.choice)));
        indexxLFsep.sqr = boolean(zeros(size(list.choice)));
        indexxLFsep.oneab = boolean(zeros(size(list.choice)));

        indexxLFsep.lab( list.choice=='h' | list.choice=='S')=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='S'& list.choice~='gammalh'& list.choice~='gammas' )=1;
        indexxLFsep.sqr(list.choice=='gammalh'| list.choice=='gammas' )=1;
        indexx('LF_noskill_sep0')=indexxLFsep;
        
    elseif indic.sep==1
        
        indexxLFsep.lab = boolean(zeros(size(list.choice)));
        indexxLFsep.exp = boolean(zeros(size(list.choice)));
        indexxLFsep.sqr = boolean(zeros(size(list.choice)));
        indexxLFsep.oneab = boolean(zeros(size(list.choice)));

        indexxLFsep.lab( list.choice=='h'| list.choice=='sff'  | list.choice=='sg' | list.choice=='sn')=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='sff'&list.choice~='sg'& list.choice~='sn'...
            &list.choice~='gammalh'& list.choice~='gammasg'& list.choice~='gammasf'& list.choice~='gammasn')=1;
        indexxLFsep.sqr(list.choice=='gammalh'| list.choice=='gammasg'| list.choice=='gammasn'| list.choice=='gammasf' )=1;

        indexx('LF_noskill_sep1')=indexxLFsep;

    end
end


while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 
    if indic.noskill==0
        if indic.sep==0
            guess_trans=trans_guess(indexx('LF'), x0, params, list.params);
        else
            guess_trans=trans_guess(indexx('LF_sep'), x0, params, list.params);
        end
    else
        guess_trans=trans_guess(indexx(sprintf('LF_noskill_sep%d', indic.sep)), x0, params, list.params);
    end
    % test
    if indic.sep==0
        error('without separate markets not yet updated')
        f=laissez_faire_nows(guess_trans, params, list, pol, laggs, indic);
    else
        f=laissez_faire_nows_sep(guess_trans, params, list, pol, laggs, indic);
    end
    %% - solving model
     lb=[];
     ub=[];
     objf=@(x)objectiveCALIBSCI(x);
     if indic.sep==0
        constrf = @(x)laissez_faire_nows_fmincon(x, params, list, pol, laggs, indic);
     else
         constrf = @(x)laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic);
     end
options = optimset('algorithm','active-set','TolCon', 1e-11,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
%[c,fval]=constrf(sol3);
% 
%- other solvers
    if indic.sep==0
        modFF = @(x)laissez_faire_nows(x, params, list, pol, laggs, indic);
    else
        modFF = @(x)laissez_faire_nows_sep(x, params, list, pol, laggs, indic);
    end
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol2, fval, exitf] = fsolve(modFF, sol3, options);

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
                LF=trans_allo_out(indexx('LF_sep'), sol3, params, list.params, indic);
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