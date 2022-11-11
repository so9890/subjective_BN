function [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth(T, list, poll, params,   symms, x0LF, init, indexx, indic, Sall, MOM, Ems)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- new initial point: drop science and technology as choices
    x0=x0LF(list.choice~='ws'&list.choice~='S' &list.choice~='gammas' &list.choice~='wsf' & list.choice~='wsn'& list.choice~='wsg' & list.choice~='gammasn'& list.choice~='gammasf'& list.choice~='gammasg'...
        &list.choice~='sff' & list.choice~='sg' & list.choice~='sn'&list.choice~='Af' & list.choice~='Ag' & list.choice~='An');
    %- new set of choice variables
    list.choice=list.choice_xgrowth;
    symms.choice=symms.choice_xgrowth;
%     list.allvars=list.allvars_xgrowth;
%     symms.allvars=symms.allvars_xgrowth;
%     
    %- new index for tranformation
        
    indexxLF.lab = boolean(zeros(size(list.choice)));
    indexxLF.exp = boolean(zeros(size(list.choice)));
    indexxLF.sqr = boolean(zeros(size(list.choice)));
    indexxLF.oneab = boolean(zeros(size(list.choice)));

    indexxLF.lab( list.choice=='hh'|list.choice=='hl')=1;
    indexxLF.exp(list.choice~='lambdaa'&list.choice~='hl'&list.choice~='hh'&list.choice~='gammalh'&list.choice~='gammall')=1;
    indexxLF.sqr(list.choice=='gammalh'|list.choice=='gammall')=1;
    
    indexx('LF_xgrowth')=indexxLF;

%% initialise stuff
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

        indexxLFsep.lab = boolean(zeros(size(list.choice)));
        indexxLFsep.exp = boolean(zeros(size(list.choice)));
        indexxLFsep.sqr = boolean(zeros(size(list.choice)));
        indexxLFsep.oneab = boolean(zeros(size(list.choice)));

        indexxLFsep.lab( list.choice=='h' )=1;
        indexxLFsep.exp(list.choice~='h'&list.choice~='gammalh'& list.choice~='lambdaa' )=1;
        indexxLFsep.sqr(list.choice=='gammalh')=1;
        indexx('LF_noskill_xgrowth')=indexxLFsep;
end

% check size of policy matrix
    [row]=size(poll);
%-- meeting emission limit 
% => extend lists etc. as tauf is choice variable
if indic.limit_LF==1
    
    syms tauf real
    symms.choice = [symms.choice, 'tauf'];
    
    list.choice = string(symms.choice);
    x0=[x0, 2];

    if indic.noskill==0
        hhelper=indexx('LF_xgrowth');
    else
        hhelper=indexx('LF_noskill_xgrowth');
    end
    hhelper.lab = [hhelper.lab, boolean(zeros(1,1))];
    hhelper.exp = [hhelper.exp, boolean(zeros(1,1))];
    hhelper.sqr = [hhelper.sqr, boolean(zeros(1,1))];
    hhelper.oneab = [hhelper.oneab, boolean(zeros(1,1))];
    if indic.noskill==0
        indexx('LF_xgrowh')=hhelper;
    else
        indexx('LF_noskill_xgrowth')=hhelper;
    end
end

%%
while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %--- emission limit per period only relevant for t>=2; for t=1 set
    % tauf=0
    if t==1 % baseline period: no emission limit, take same as in second period
        Emlim=0;
    else % 
        Emlim = Ems (t-1);
    end
    % read in policy 
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
    if indic.noskill==0
        guess_trans=trans_guess(indexx('LF_xgrowth'), x0, params, list.params);
    else
        guess_trans=trans_guess(indexx(sprintf('LF_noskill_xgrowth')), x0, params, list.params);
    end
    % test
    f=laissez_faire_xgrowth(guess_trans, params, list, pol, laggs, indic, MOM,t, Emlim);
    
    %% - solving model
 lb=[];
 ub=[];
 objf=@(x)objectiveCALIBSCI(x);
 constrf = @(x)laissez_faire_xgrowth_fmincon(x, params, list, pol, laggs, indic, MOM,t, Emlim);

options = optimset('algorithm','sqp','TolCon', 1e-11,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
[sol3,fval,exitflag,output,lambda] = fmincon(objf,sol3,[],[],[],[],lb,ub,constrf,options);

%[c,fval]=constrf(sol3);
% 
%- other solvers
    modFF = @(x)laissez_faire_xgrowth(x, params, list, pol, laggs, indic, MOM,t, Emlim);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol2, fval, exitf] = fsolve(modFF, sol3, options);
    [sol2, fval, exitf] = fsolve(modFF,sol2, options);

    % pass to standard algorithmop
%      options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,);%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, x1, options);
% 
     options = optimoptions('fsolve', 'TolFun', 10e-11, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
     [sol3, fval, exitf] = fsolve(modFF, sol2, options);


    %- transform results to bounded variables
    if indic.noskill==0
        LF=trans_allo_out(indexx('LF_xgrowth'), sol3, params, list.params, indic);
     else
        LF=trans_allo_out(indexx('LF_noskill_xgrowth'), sol3, params, list.params, indic);
    end
    %% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    if t>1
        [LF_SIM(:,t-1), An0, Ag0,Af0] =aux_solutionLF_xgrowth( SLF, pol, laggs, list, symms, indexx, params, indic, MOM,t,Emlim);
        FVAL(t-1)=max(abs(fval));
    else
        [~, An0, Ag0,Af0] =aux_solutionLF_xgrowth( SLF, pol, laggs, list, symms, indexx, params, indic, MOM,t,Emlim);
    end
    %% - update for next round
    x0 = LF; % initial guess
    laggs=eval(symms.init);
    t=t+1;
end
LF_SIM=LF_SIM';
end