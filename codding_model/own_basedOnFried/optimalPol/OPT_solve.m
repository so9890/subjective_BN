function [symms, list, op_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;

% symbilic variables and lists
syms hhf hhg hlf hlg C F G Af Ag An hl hh real
symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh];
list.opt  = string(symms.opt); 

nn= length(list.opt); % number of variables

%%% Initial Guess %%%
% uses social planner of LF allocation 
%%%%%%%%%%%%%%%%%%%%%
if indic.target==1

    helper=load('SP_target');
    sp_all=helper.sp_all;
    
    if ~isvarname('sp_all')
        error('did not load sp solution')
    end 

    x0 = zeros(nn*T,1);

    x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =sp_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =sp_all(:,list.allvars=='hhg'); % hhg
    x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =sp_all(:,list.allvars=='hlf'); % hlf
    x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =sp_all(:,list.allvars=='hlg'); % hlg 
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =sp_all(:,list.allvars=='C');   % C
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =sp_all(:,list.allvars=='F');

    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =sp_all(:,list.allvars=='G');   % G
    x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =sp_all(:,list.allvars=='Af');  % Af
    x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =sp_all(:,list.allvars=='Ag');  % Ag
    x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =sp_all(:,list.allvars=='An');  % An
    x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =sp_all(:,list.allvars=='hl');  % hl
    x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =sp_all(:,list.allvars=='hh');  % hh
   
    
    % initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
    
 
elseif indic.target==0
        % for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
        % to get initial guess
        taus=0;
        tauf=0;
        taul=0;
        lambdaa=1; % balances budget with tauf= taul=0
        pol=eval(symms.pol);

        if ~isfile('FB_LF_SIM_NOTARGET.mat')
            [LF_SIM, polLF] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx);   
            save('FB_LF_SIM_NOTARGET','LF_SIM');
        else
            help=load('FB_LF_SIM_NOTARGET.mat');
            LF_SIM=help.LF_SIM;
        end

    x0 = zeros(nn*T,1);

    x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
    x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
    x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
    x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =LF_SIM(list.allvars=='C',1:T);   % C
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =LF_SIM(list.allvars=='F',1:T);
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =LF_SIM(list.allvars=='G',1:T);   % G
    x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
    x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
    x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An
    x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
    x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh

%     % initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
%     % Laissez faire solution, not == init (which refers to 2010-2014)!
%     Ag0 = LF_SIM(list.allvars=='Ag',1);
%     Af0 = LF_SIM(list.allvars=='Af',1);
%     An0 = LF_SIM(list.allvars=='An',1);
%     init201519= eval(symms.init);
end

%%% Transform to unbounded variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- most of variables bounded by zero below
guess_trans=log(x0);

%- exceptions with upper bound; hl, hh, F in case of target
guess_trans(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))))./...
    x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))));
guess_trans(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))))./...
    x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))));
if indic.target==1
    guess_trans(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))));
end

lb=[];
ub=[];



%%% Test Constraints and Objective Function %%%
f =  objective(guess_trans, T, params, list, Ftarget, indic);
[c, ceq] = constraints(guess_trans, T, params, init201519, list, Ems, indic);

%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

objf=@(x)objective(x, T, params, list, Ftarget, indic);
constf=@(x)constraints(x, T, params, init201519, list, Ems, indic);
options = optimset('algorithm','sqp','Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);

% options = optimset('algorithm','sqp','TolCon',1e-2,'Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
% options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
% without target the SP and OPTIMAL POL are the best

if x==guess_trans
    fprintf('In version target=%d, the initial guess and the OPtimal pol are the same. With target=0 this is the LF and with target=1 it is the SP one.', indic.target);
end

% options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% [x2,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
% 
% x=savebfp;%  output.bestfeasible.x;
% options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);

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
 
end