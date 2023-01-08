function [LF_SIM, poll, FVAL, indexx] = solve_LF_nows(T, list, poll, params, symms, x0LF, init, indexx, indic, Sall, Ems)
% simulate economy under laissez faire
indic.xgrowth=0;
indic.ineq=0;
% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%-- to save results
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);

%-- initialise values
x0      = x0LF; % initial guess from calibration 
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

    indexxLFsep.lab = boolean(zeros(size(list.choice)));
    indexxLFsep.exp = boolean(zeros(size(list.choice)));
    indexxLFsep.sqr = boolean(zeros(size(list.choice)));
    indexxLFsep.oneab = boolean(zeros(size(list.choice)));

    if indic.sep==1
        indexxLFsep.lab( list.choice=='h'| list.choice=='sff'  | list.choice=='sg' | list.choice=='sn')=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='sff'&list.choice~='sg'& list.choice~='sn'...
            &list.choice~='gammalh'& list.choice~='gammasg'& list.choice~='gammasf'& list.choice~='gammasn' ...
            & list.choice~='lambdaa')=1;
        indexxLFsep.sqr(list.choice=='gammalh'| list.choice=='gammasg'| list.choice=='gammasn'| list.choice=='gammasf' )=1;
    else
        indexxLFsep.lab( list.choice=='h'| list.choice=='S')=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='S'...
            &list.choice~='gammalh'& list.choice~='gammas' & list.choice~='lambdaa')=1;
        indexxLFsep.sqr(list.choice=='gammalh'| list.choice=='gammas' )=1;
    end
        indexx(sprintf('LF_noskill_sep%d', indic.sep))=indexxLFsep;
        
end

if indic.sep==2 % partial equilibrium
    x0=x0(list.choice~='wsf'&list.choice~='wsg'& list.choice~='wsn'& list.choice~='ws' &...
        list.choice~='gammasf'&list.choice~='gammasg'& list.choice~='gammasn'& list.choice~='gammas'& list.choice~='S');
    
    symms.choice= symms.choice(list.choice~='wsf'&list.choice~='wsg'& list.choice~='wsn'...
                     &list.choice~='gammasf'&list.choice~='gammasg'& list.choice~='gammasn' ...
                     & list.choice~='gammas' & list.choice~='ws' & list.choice~='S' );
                 % sff sg sn are still endogenous and chosen by machine
                 % producers
    list.choice=string(symms.choice);

    indexxLFsep2.lab = boolean(zeros(size(list.choice)));
    indexxLFsep2.exp = boolean(zeros(size(list.choice)));
    indexxLFsep2.sqr = boolean(zeros(size(list.choice)));
    indexxLFsep2.oneab = boolean(zeros(size(list.choice)));

    indexxLFsep2.lab( list.choice=='hh'|list.choice=='hl'| list.choice=='sff'  | list.choice=='sg' | list.choice=='sn')=1;
    indexxLFsep2.exp(list.choice~='hh'&list.choice~='hl'&list.choice~='sff'&list.choice~='sg'& list.choice~='sn'...
        &list.choice~='gammalh'&list.choice~='gammall'& list.choice~='lambdaa')=1;
    indexxLFsep2.sqr(list.choice=='gammalh' | list.choice=='gammall' )=1;

    indexx('LF_sep2')=indexxLFsep2;
end
    
if indic.sep==3 % free movement energy scientists
    x0=x0(list.choice~='gammasf');
    syms se real
    symms.choice= [symms.choice(list.choice~='gammasf'), se];
    list.choiceold=list.choice;
    list.choice=string(symms.choice);
    x0=[x0,Sall.sff+Sall.sg]; % order has to match how variables are added to symms.choice!
    hhelper=indexx('LF');
    hhelper.lab = [hhelper.lab(list.choiceold~='gammasf'), boolean(ones(1,1))];
    hhelper.exp = [hhelper.exp(list.choiceold~='gammasf'), boolean(zeros(1,1))];
    hhelper.sqr = [hhelper.sqr(list.choiceold~='gammasf'), boolean(zeros(1,1))];
    hhelper.oneab = [hhelper.oneab(list.choiceold~='gammasf'), boolean(zeros(1,1))];
    indexx('LF_sep3')=hhelper;

end
    
%%
%-- check size of policy matrix
[row]=size(poll);

%-- meeting emission limit 
% => extend lists etc. as tauf is choice variable
if indic.limit_LF==1
    
    syms tauf real
    symms.choice = [symms.choice, 'tauf'];
    
    list.choice = string(symms.choice);
    x0=[x0, 2];

    if indic.sep<=2
        if indic.noskill==0
            hhelper=indexx('LF');
        else
            hhelper=indexx(sprintf('LF_noskill_sep%d', indic.sep));
        end
    else
        hhelper=indexx(sprintf('LF_sep%d', indic.sep));
    end
    hhelper.lab = [hhelper.lab, boolean(zeros(1,1))];
    hhelper.exp = [hhelper.exp, boolean(zeros(1,1))];
    hhelper.sqr = [hhelper.sqr, boolean(zeros(1,1))];
    hhelper.oneab = [hhelper.oneab, boolean(zeros(1,1))];
    
    if indic.sep<=2
    if indic.noskill==0
        
        indexx('LF')=hhelper;
    else
        indexx(sprintf('LF_noskill_sep%d', indic.sep))=hhelper;
    end
    else
        indexx(sprintf('LF_sep%d', indic.sep))=hhelper;
    end
end

%update for Sun

%     hhelper=indexx('LF');
%     hhelper.Sc=boolean(zeros(size(hhelper.lab)));
%     hhelper.Sc(list.choice=='S' )=1;
%     hhelper.lab(list.choice=='S' )=0;
%     indexx('LF')=hhelper;

%%
while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %--- emission limit per period only relevant for t>=2; for t=1 set
    % tauf=0
    if t==1 % baseline period: no emission limit, take same as in second period
        Emlim=0;
    else % 
        Emlim = Ems(t-1);
    end
    %--- read in policy (can be vector or static)
    if row(1)>1
        if t<=T
            pol=poll(t,:);
        else
            pol=poll(t-1,:);
        end
    else
        pol=poll;
    end
    %% - transforming variables to unbounded variables
    %-- index for transformation 
    if indic.sep<=1
        if indic.noskill==0
            guess_trans=trans_guess(indexx('LF'), x0, params, list.params);
        else
            guess_trans=trans_guess(indexx(sprintf('LF_noskill_sep%d', indic.sep)), x0, params, list.params);
        end
    else
        guess_trans=trans_guess(indexx(sprintf('LF_sep%d',indic.sep)), x0, params, list.params);        
    end
    
    % test
    if indic.sep<=1
        f=laissez_faire_nows_sep(guess_trans, params, list, pol, laggs, indic, Emlim, t);
    elseif indic.sep == 2 % partial equilibrium version, wsf wsg and wsn are fixed
        f=laissez_faire_nows_partialS(guess_trans, params, list, pol, laggs, indic, Emlim, t);
    elseif indic.sep==3
        f=laissez_faire_nows_sepSe(guess_trans, params, list, pol, laggs, indic, Emlim, t);
    end
    %% - solving model
     lb=[];
     ub=[];
     if ~(indic.limit_LF==0 && indic.notaul==4 && indic.noknow_spill==0)
         objf=@(x)objectiveCALIBSCI(x);
    %     if indic.sep==0
      %      constrf = @(x)laissez_faire_nows_fmincon(x, params, list, pol, laggs, indic);
     %    else
             constrf = @(x)laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic, Emlim, t);
     %    end
         if indic.labshareequ==0 && (~indic.taul0==1)
            options = optimset('algorithm','active-set','TolCon', 1e-7,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         else
            options = optimset('algorithm','sqp','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         end
        [sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
     end
    %[c,fval]=constrf(sol3);
% 
%- other solvers
    if indic.sep<=1
        modFF = @(x)laissez_faire_nows_sep(x, params, list, pol, laggs, indic, Emlim, t);
    elseif indic.sep==2
        modFF = @(x)laissez_faire_nows_partialS(x, params, list, pol, laggs, indic, Emlim, t);
    elseif indic.sep==3
        modFF = @(x)laissez_faire_nows_sepSe(x, params, list, pol, laggs, indic, Emlim, t);
    end
    options = optimoptions('fsolve', 'TolFun', 10e-9,'Display','iter', 'MaxFunEvals',8e6, 'MaxIter', 3e6,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
        if (indic.limit_LF==0 && indic.notaul==4 && indic.noknow_spill==0)
            [sol3, fval, exitf] = fsolve(modFF, guess_trans, options);
        end
    [sol2, fval, exitf] = fsolve(modFF, sol3, options);

    % pass to standard algorithm
%      options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,);%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, x1, options);
% 
% if ~(indic.noskill==1 && indic.tauf==1 && indic.xgrowth==0)
if indic.sep<=1
    if ~(indic.sizeequ==1 && indic.GOV==0 && indic.noknow_spill==0 ) &&  (indic.labshareequ==1&& indic.GOV==0 && indic.noknow_spill==0 )
         options = optimoptions('fsolve', 'TolFun', 10e-9, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
         [sol3, fval, exitf] = fsolve(modFF, sol2, options);
    else
        sol3=sol2;
    end
else
    sol3=sol2;
end
 
if exitf<=0
    error('code did not solve')
end
%- transform results to bounded variables
if indic.sep<=1
    if indic.noskill==0
            LF=trans_allo_out(indexx('LF'), sol3, params, list.params, indic);
%         else
%                 LF=trans_allo_out(indexx('LF_sep'), sol3, params, list.params, indic);
%         end
     else
        LF=trans_allo_out(indexx(sprintf('LF_noskill_sep%d', indic.sep)), sol3, params, list.params, indic);
     end
else
         LF=trans_allo_out(indexx(sprintf('LF_sep%d', indic.sep)), sol3, params, list.params, indic);
end
%% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    if t>1
%         if indic.sep==0
%             LF_SIM(:,t-1)=aux_solutionLF( SLF, pol, laggs, list, symms, indexx, params, indic);
%         else
            LF_SIM(:,t-1)=aux_solutionLF_sep( SLF, pol, laggs, list, symms, indexx, params, indic, Emlim, t);
%         end
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