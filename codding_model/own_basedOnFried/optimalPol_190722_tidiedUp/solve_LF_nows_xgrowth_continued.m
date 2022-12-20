function [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_xgrowth_continued(T, list, pol, params, Sparams,  symms, x0LF, init, indexx, indic, Sall)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %- new initial point: drop science and technology as choices
    x0=x0LF(list.choice~='ws' & list.choice~='S' & list.choice~='gammas'&list.choice~='sff' &...
            list.choice~='sg' & list.choice~='sn'&list.choice~='Af' & list.choice~='Ag' & list.choice~='An');
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
    indexxLF.exp(list.choice~='hl'&list.choice~='hh'&list.choice~='gammalh'&list.choice~='gammall')=1;
    indexxLF.sqr(list.choice=='gammalh'|list.choice=='gammall')=1;
    
    indexx('LF_xgrowth')=indexxLF;
    
    %-version with BN
   indexxLF.BN = boolean(zeros(size(list.choice)));
   indexxLF.exp = boolean(zeros(size(list.choice)));

    indexxLF.exp(list.choice~='C' & list.choice~='hl'&list.choice~='hh'&list.choice~='gammalh'&list.choice~='gammall')=1;
    indexxLF.BN(list.choice=='C')=1;
    indexx('LF_xgrowth_BN')=indexxLF;

%% initialise stuff
%-- to save results
% symms.allvars=[symms.allvars, gammasn, gammasg, gammasf];
% list.allvars=string(symms.allvars);
if indic.sep==1
    list.allvars=list.sepallvars;
    symms.allvars=symms.sepallvars;
end
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
        indexxLFsep.exp(list.choice~='h'&list.choice~='gammalh' )=1;
        indexxLFsep.sqr(list.choice=='gammalh')=1;
        indexx('LF_noskill_xgrowth')=indexxLFsep;
end

%%
while t<=T % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 
    if indic.noskill==0
       guess_trans=trans_guess(indexx('LF_xgrowth'), x0, params, list.params);
            
    else
        guess_trans=trans_guess(indexx(sprintf('LF_noskill_xgrowth')), x0, params, list.params);
    end
    % test
    f=laissez_faire_xgrowth(guess_trans, params, list, pol, laggs, indic);
    
    %% - solving model
 lb=[];
 ub=[];
 objf=@(x)objectiveCALIBSCI(x);
 constrf = @(x)laissez_faire_xgrowth_fmincon(x, params, list, pol, laggs, indic);

options = optimset('algorithm','active-set','TolCon', 1e-11,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);
%[c,fval]=constrf(sol3);
% 
%- other solvers
    modFF = @(x)laissez_faire_xgrowth(x, params, list, pol, laggs, indic);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol2, fval, exitf] = fsolve(modFF, sol3, options);

    % pass to standard algorithmop
%      options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,);%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, x1, options);
% 
     options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
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
        LF_SIM(:,t)=aux_solutionLF_xgrowth(Sparams, SLF, pol, laggs, list, symms, indexx, params, indic);
        FVAL(t)=max(abs(fval));
    %% - update for next round
    x0 = LF; % initial guess
        read_in_params;
        An0=(1+vn)*laggs(list.init=='An0');
        Ag0=(1+vg)*laggs(list.init=='Ag0');
        Af0=(1+vf)*laggs(list.init=='Af0');
    laggs=eval(symms.init);
    t=t+1;
end
end