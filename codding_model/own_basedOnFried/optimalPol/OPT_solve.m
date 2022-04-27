syms hhf hhg hlf hlg C F G Af Ag An hl hh real
symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh];
list.opt  = string(symms.opt); 

nn= length(list.opt); % number of variables

%%% Initial Guess %%%
%%%%%%%%%%%%%%%%%%%%%

% if indic.target==0
%     % for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
%     taus=0;
%     tauf=0;
%     taul=0;
%     lambdaa=1; % balances budget with tauf= taul=0
%     pol=eval(symms.pol);
%     [out_trans, pol] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init, indexx);   
% end
% 
% x0 = zeros(nn*T,1);
% 
% x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =out_trans(list.sp=='hhf',2:T+1); % hhf; first period in LF is baseline
% x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =out_trans(list.sp=='hhg',2:T+1); % hhg
% x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =out_trans(list.sp=='hlf',2:T+1); % hlf
% x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =out_trans(list.sp=='hlg',2:T+1); % hlg 
% x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =out_trans(list.sp=='C',2:T+1);   % C
% if indic.target==1
%     upF=(1-1e-10)*(Ems+params(list.params=='deltaa'))./params(list.params=='omegaa');
%     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))) =upF;   % F
% else
%     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))) =out_trans(list.sp=='F',2:T+1);
% end
% x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =out_trans(list.sp=='G',2:T+1);   % G
% x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =out_trans(list.sp=='Af',2:T+1);  % Af
% x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =out_trans(list.sp=='Ag',2:T+1);  % Ag
% x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =out_trans(list.sp=='An',2:T+1);  % An
% x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =out_trans(list.sp=='hl',2:T+1);  % hl
% x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =out_trans(list.sp=='hh',2:T+1);  % hh

if indic.target == 0
    load('SP_notarget')
else
    load('SP_target')
end
    

x0 = zeros(nn*T,1);

x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =out_trans(list.sp=='hhf',2:T+1); % hhf; first period in LF is baseline
x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =out_trans(list.sp=='hhg',2:T+1); % hhg
x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =out_trans(list.sp=='hlf',2:T+1); % hlf
x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =out_trans(list.sp=='hlg',2:T+1); % hlg 
x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =out_trans(list.sp=='C',2:T+1);   % C
if indic.target==1
    upF=(1-1e-10)*(Ems+params(list.params=='deltaa'))./params(list.params=='omegaa');
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))) =upF;   % F
else
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))) =out_trans(list.sp=='F',2:T+1);
end
continue finding an expression for G
x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =out_trans(list.sp=='G',2:T+1);   % G
x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =out_trans(list.sp=='Af',2:T+1);  % Af
x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =out_trans(list.sp=='Ag',2:T+1);  % Ag
x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =out_trans(list.sp=='An',2:T+1);  % An
x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =out_trans(list.sp=='hl',2:T+1);  % hl
x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =out_trans(list.sp=='hh',2:T+1);  % hh

%%% Linear constraints optimisation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
guess_trans=log(x0);
guess_trans(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))))./...
 x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))));
guess_trans(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))))./...
 x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))));
if indic.target==1
    guess_trans(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))=log((upF-x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))));
end

lb=[];
ub=[];

% lb = zeros((nn*T),1);
% ub = Inf*ones((nn*T),1);
% 
% ub((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = params(list.params=='upbarH');
% ub((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = params(list.params=='upbarH');
% 
% if indic.target==1
%   ub((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T) = (Ems'+params(list.params=='deltaa'))./params(list.params=='omegaa');
% end

% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution, not == init (which refers to 2010-2014)!
Ag0 = out_trans(list.sp=='Ag',1);
Af0 = out_trans(list.sp=='Af',1);
An0 = out_trans(list.sp=='An',1);
initOPT= eval(symms.init);
%%% Test Constraints and Objective Function %%%
f =  objective(x0,T,params, list);
[c, ceq] = constraints(x0, T, params, initOPT, list, Ems);

%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

objf=@(x)objective(x,T,params, list);
constf=@(x)constraints(x, T, params, initOPT, list, Ems);
options = optimset('algorithm','sqp','Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
% options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constf,options);

%- test if output equals LF solution (in no target version)
if indic.target==0
if min(x0==output.bestfeasible.x)~=1
    error('optimal policy is not the theoretically optimal one!')
end
end

x0 = x;
options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(@(x)COMET_Objective(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip),x0,[],[],[],[],lb,ub,@(x)COMET_Constraints(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,b3,B0,energy_wedge_fix),options);

%Note: Check exit flags to ensure convergence. If program stops for other reasons
%(e.g., max iterations exceeded), it needs to be re-started, either from
%the current guess, or, in case of failure, from a revised guess.
 
