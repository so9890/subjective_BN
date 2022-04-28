function [symms, list, opt_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

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
    kappaa=(1-3e-4); % to scale sp result
    x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =kappaa*sp_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =kappaa*0.2*sp_all(:,list.allvars=='hhg'); % hhg
    x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =kappaa*sp_all(:,list.allvars=='hlf'); % hlf
    x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =kappaa*1.2*sp_all(:,list.allvars=='hlg'); % hlg 
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =kappaa*0.4*sp_all(:,list.allvars=='C');   % C
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =kappaa*0.4*sp_all(:,list.allvars=='F');

    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =kappaa*0.6*sp_all(:,list.allvars=='G');   % G
    x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =sp_all(:,list.allvars=='Af');  % Af
    x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =sp_all(:,list.allvars=='Ag');  % Ag
    x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =sp_all(:,list.allvars=='An');  % An
    x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =kappaa*0.4*sp_all(:,list.allvars=='hl');  % hl
    x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =kappaa*0.3*sp_all(:,list.allvars=='hh');  % hh
   
    
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
            [LF_SIM] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init201519, indexx);   
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
f =  objective(guess_trans, T, params, list, Ftarget, indic)
[c, ceq] = constraints(guess_trans, T, params, init201519, list, Ems, indic)
for i=1:length(symms.opt)
model(i)=jacobian(f, symms.opt(i)); % should give derivative=0 if not present
end
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

% transform
out_trans=exp(x);
out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

% save results
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, wsn, wsg, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf]= OPT_aux_vars(out_trans, list, params, T, init201519);
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
 
 opt_all=eval(symms.allvars);

if indic.target==1
    save('OPT_target', 'opt_all')
else
    save('OPT_notarget', 'opt_all')
end

end