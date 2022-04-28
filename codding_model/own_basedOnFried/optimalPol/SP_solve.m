function [symms, list, sp_all]=SP_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F real
symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F];
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

% FIRST read in efficient allocation without emission target
% in SECOND rescale so that emission target is reached.
% for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
% which is optimal as competitive economy is efficient
taus=0;
tauf=0;
taul=0;
lambdaa=1; % balances budget with tauf= taul=0
pol=eval(symms.pol);
if ~isfile('FB_LF_SIM_NOTARGET.mat')
    [LF_SIM, polLF] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init201519, indexx);   
    save('FB_LF_SIM_NOTARGET','LF_SIM');
    if pol~=polLF
        error('LF not solved under fb policy');
    end
else
     help=load('FB_LF_SIM_NOTARGET.mat');
     LF_SIM=help.LF_SIM;
end

if indic.target==0
    x0 = zeros(nn*T,1);
    Ftarget = 0; % placeholder
    x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
    x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
    x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =LF_SIM(list.allvars=='hhn',1:T); % hhg
    x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
    x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
    x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =LF_SIM(list.allvars=='hln',1:T); % hlg 
    x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =LF_SIM(list.allvars=='xf',1:T); % hlf
    x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =LF_SIM(list.allvars=='xg',1:T); % hlg 
    x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =LF_SIM(list.allvars=='xn',1:T); % hlg 
    x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
    x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
    x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An
    x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
    x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
    x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',1:T);  % C
    x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =LF_SIM(list.allvars=='F',1:T);  % C

    % Rescale if model with emission target %%
elseif indic.target==1  
    Ftarget = (Ems+deltaa)/omegaa;
    kappaa = Ftarget./LF_SIM(list.allvars=='F',1:T); % ratio of targeted F to non-emission
    kappaa = kappaa*(1-1e-10);
    x0 = zeros(nn*T,1);
    x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =kappaa.*LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
    x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =kappaa.*LF_SIM(list.allvars=='hhg',1:T); % hhg
    x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =kappaa.*LF_SIM(list.allvars=='hhn',1:T); % hhg
    x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =kappaa.*LF_SIM(list.allvars=='hlf',1:T); % hlf
    x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =kappaa.*LF_SIM(list.allvars=='hlg',1:T); % hlg 
    x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =kappaa.*LF_SIM(list.allvars=='hln',1:T); % hlg 
    x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =kappaa.*LF_SIM(list.allvars=='xf',1:T); % hlf
    x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =kappaa.*LF_SIM(list.allvars=='xg',1:T); % hlg 
    x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =kappaa.*LF_SIM(list.allvars=='xn',1:T); % hlg 
    x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
    x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
    x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An
    x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =kappaa.*LF_SIM(list.allvars=='hl',1:T);  % hl
    x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =kappaa.*LF_SIM(list.allvars=='hh',1:T);  % hh
    x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =kappaa.*LF_SIM(list.allvars=='C',1:T);  % C
    x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =kappaa.*LF_SIM(list.allvars=='F',1:T);  % C
end
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

 f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic);

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic);

%  options = optimoptions('Algorithm','sqp','TolStep',1e-10,'TolFun',1e-16,'MaxFunEvals',500000,'MaxIter',6200,'Display','Iter','MaxSQPIter',10000);

options = optimset('algorithm','sqp','TolX', 1e-8, 'Tolfun',1e-16,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
[x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);

if x==guess_trans
    fprintf('In version target=%d, the LF and FB are the same.', indic.target);
end
%   if exitflag==2 && indic.target==1 %(otherwise does not solve)
%     options = optimset('algorithm','active-set','Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     [x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);
%   end

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
            pn, pg, pf, pee, wh, wl, ws, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf]= SP_aux_vars(out_trans, list, params, T, init201519);
gammall = zeros(size(pn));
gammalh = zeros(size(pn));

sp_all=eval(symms.allvars);

if indic.target==1
    save('SP_target', 'sp_all')
else
    save('SP_notarget', 'sp_all')
end
end