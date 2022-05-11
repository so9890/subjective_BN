function [symms, list, sp_all]=SP_solve_est(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems)
addpath('../optimalPol')
% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn real
symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn ];
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%


% helper=load('SP_solution_wse.mat');
helper=load('SP_notarget_405.mat');
sp_all=helper.sp_all; % use solution as starting value

    x0 = zeros(nn*T,1);
    Ftarget = (Ems+deltaa)/omegaa;
    kappaa = Ftarget'./sp_all(1:T, list.allvars=='F'); % ratio of targeted F to non-emission
    if indic.target==1
        kappaa = kappaa*(1-1e-10);
    else
        kappaa =1; % no scaling! 
    end
    x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =kappaa.*sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =kappaa.*sp_all(1:T, list.allvars=='hhg'); % hhg
    x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =kappaa.*sp_all(1:T, list.allvars=='hhn'); % hhg
    x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =kappaa.*sp_all(1:T, list.allvars=='hlf'); % hlf
    x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =kappaa.*sp_all(1:T, list.allvars=='hlg'); % hlg 
    x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =kappaa.*sp_all(1:T, list.allvars=='hln'); % hlg 
    x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =kappaa.*sp_all(1:T, list.allvars=='xf'); % hlf
    x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =kappaa.*sp_all(1:T, list.allvars=='xg'); % hlg 
    x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =kappaa.*sp_all(1:T, list.allvars=='xn'); % hlg 
    x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
    x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
    x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
    x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =kappaa.*sp_all(1:T, list.allvars=='hl');  % hl
    x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =kappaa.*sp_all(1:T, list.allvars=='hh');  % hh
    x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =kappaa.*sp_all(1:T, list.allvars=='C');  % C
    x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =kappaa.*sp_all(1:T, list.allvars=='F');  % C
    x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))     =kappaa.*sp_all(1:T, list.allvars=='sg');  % C
    x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))     =kappaa.*sp_all(1:T, list.allvars=='sn');  % C
    x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff')))     =kappaa.*sp_all(1:T, list.allvars=='sff');  % C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution; do not take init which refers to 2010-2014!
% this version here under assumption of first best policy in initial policy
% that is: tauf=0, taul=0, taus=0, lambdaa to balage budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ag0 = LF_SIM(list.allvars=='Ag',1);
% Af0 = LF_SIM(list.allvars=='Af',1);
% An0 = LF_SIM(list.allvars=='An',1);
% 
initOPT= init201519; % as calibrated under BAU policy

%%% Transform variables to unbounded vars => requires less constraints! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 guess_trans=log(x0);
 guess_trans(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))))./...
     x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))));
 guess_trans(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))))./...
     x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))));
if indic.target==1
    guess_trans(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))=log((Ftarget'-x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F'))))./...
     x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F'))));
end
lb=[];
ub=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic, initOPT);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic);

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic);

%  options = optimoptions('Algorithm','sqp','TolStep',1e-10,'TolFun',1e-16,'MaxFunEvals',500000,'MaxIter',6200,'Display','Iter','MaxSQPIter',10000);

options = optimset('algorithm','sqp', 'TolCon',1e-6, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
[x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
vv =output.bestfeasible.x
%gg=load('SP_solution_wse_withT');
% if abs(x-guess_trans)<1e-7
%     fprintf('In version target=%d, the LF and FB are the same.', indic.target);
% end

%
 %if exitflag==2  %(otherwise does not solve)
 options = optimset('algorithm','active-set','TolCon',1e-6,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     [xas,fval,exitflag,output,lambda] = fmincon(objfSP,gg.x,[],[],[],[],lb,ub,constfSP,options);
       [xas,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
       save('active-set_solu_target_test')
% end
x=ss
out_trans=exp(x);
out_trans((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)=upbarH./(1+exp(x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
out_trans((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)=upbarH./(1+exp(x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
if indic.target==1
    out_trans((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)=Ftarget'./(1+exp(x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)));
end
% additional variables 
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wse, wsn, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF]= SP_aux_vars_2S(out_trans, list, params, T, init201519);
gammall = zeros(size(pn));
gammalh = zeros(size(pn));

sp_all=eval(symms.allvars);

if indic.target==1
    save('SP_target_bestfeasible', 'sp_all')
else
    save('SP_notarget', 'sp_all')
end
end