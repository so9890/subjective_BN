% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl real
symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl];
list.sp  = string(symms.sp); 
nn= length(list.sp); 

% %-- indexx
% indexxsp.lab = boolean(zeros(length(list.sp)*T,1));
% indexxsp.exp = boolean(zeros(length(list.sp)*T,1));
% indexxsp.sqr = boolean(zeros(length(list.sp)*T,1));
% indexxsp.oneab = boolean(zeros(length(list.sp)*T,1));
% 
% indexxsp.exp=ones(length(list.sp)*T,1);
% 
% indexx('sp')=indexxsp;
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
if ~isfile('OPT_LF_SIM_NOTARGET.mat')
    [LF_SIM, pol] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init, indexx);   
    save('OPT_LF_SIM_NOTARGET','LF_SIM');
else
    help=load('OPT_LF_SIM_NOTARGET.mat');
    LF_SIM=help.LF_SIM;
end


x0 = zeros(nn*T,1);
x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =LF_SIM(list.allvars=='hhf',2:T+1); % hhf; first period in LF is baseline
x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =LF_SIM(list.allvars=='hhg',2:T+1); % hhg
x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =LF_SIM(list.allvars=='hhn',2:T+1); % hhg
x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =LF_SIM(list.allvars=='hlf',2:T+1); % hlf
x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',2:T+1); % hlg 
x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =LF_SIM(list.allvars=='hln',2:T+1); % hlg 
x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =LF_SIM(list.allvars=='xf',2:T+1); % hlf
x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =LF_SIM(list.allvars=='xg',2:T+1); % hlg 
x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =LF_SIM(list.allvars=='xn',2:T+1); % hlg 
x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',2:T+1);  % Af
x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',2:T+1);  % Ag
x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',2:T+1);  % An
x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',2:T+1);  % hl
x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',2:T+1);  % hh
x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',2:T+1);  % C

% Rescale if model with emission target %%
if indic.target==1  
    Ftarget = Ems 
    kappaa = ; % ratio of targeted F to non-emission
    x0 = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution; do not take init which refers to 2010-2014!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ag0 = LF_SIM(list.allvars=='Ag',1);
Af0 = LF_SIM(list.allvars=='Af',1);
An0 = LF_SIM(list.allvars=='An',1);

initOPT= eval(symms.init);

%%% Transform variables to unbounded vars => requires less constraints! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guess_trans=log(x0);
% guess_trans(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))))./...
%     x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))));
% guess_trans(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))))./...
%     x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))));

lb = zeros((nn*T),1);
ub = Inf*ones((nn*T),1);

ub((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = params(list.params=='upbarH');
ub((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = params(list.params=='upbarH');

% if indic.target==1
%   ub((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T) = (Ems'+params(list.params=='deltaa'))./params(list.params=='omegaa');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f =  objectiveSP(x0,T,params, list);
[c, ceq] = constraintsSP(x0, T, params, initOPT, list, Ems, indic);

objfSP=@(x)objectiveSP(x,T,params, list);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic);

options = optimset('algorithm','sqp','Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
% options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,x0,[],[],[],[],lb,ub,constfSP,options);
