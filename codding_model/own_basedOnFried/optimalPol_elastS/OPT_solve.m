function [symms, list, opt_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;
% indic.taus=0; % == 0 if no taus allowed
 
% symbilic variables and lists
syms hhf hhg hlf hlg C F G Af Ag An hl hh S sg real
% if indic.target==0
symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh S sg];
list.opt  = string(symms.opt); 
% else
%     symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh sn sff sg taus];
%     list.opt  = string(symms.opt); 
% end
nn= length(list.opt); % number of variables

%%% Initial Guess %%%
% uses social planner of LF allocation 
%%%%%%%%%%%%%%%%%%%%%
if indic.target==1

    helper=load(sprintf('SP_target_active_set_0505_spillover%d.mat', indic.spillovers));
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
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =sp_all(:,list.allvars=='S'); 
    x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =sp_all(:,list.allvars=='sg'); 

    % initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
  
elseif indic.target==0
        % for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
        % to get initial guess
        taus=0;
        tauf=0;
        taul=0;
        lambdaa=1; % balances budget with tauf= taul=0
        pol=eval(symms.pol);

        if ~isfile(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers))
            [LF_SIM, polLF, FVAL] =solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx);
            save(sprintf('FB_LF_SIM_NOTARGET_spillover%d', indic.spillovers),'LF_SIM');
            helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers));
        %      LF_SIM=help.LF_SIM;
             [LF_SIM]=solve_LF_VECT(T, list, pol, params,symms, init201519, helper)
            if pol~=polLF
                error('LF not solved under fb policy');
            end
        else
             helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers));
             [LF_SIM]=solve_LF_VECT(T, list, pol, params,symms, init201519, helper);
             LF_SIM=LF_SIM';
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
     x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =LF_SIM(list.allvars=='S',1:T);  % hh
     x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % hh
     
    helper=load(sprintf('SP_notarget_active_set_0505_spillover%d.mat', indic.spillovers));
    sp_all=helper.sp_all;
    
    if ~isvarname('sp_all')
        error('did not load sp solution')
    end 
% 
   
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
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =(1-1e-10)*sp_all(:,list.allvars=='S'); 
    x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =sp_all(:,list.allvars=='sg');

end

%%% Transform to unbounded variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- most of variables bounded by zero below
guess_trans=log(x0);

% guess_trans(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus')))=x0(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus'))) ;
%- exceptions with upper bound; hl, hh, F in case of target
guess_trans(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))))./...
    x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))));
guess_trans(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))))./...
    x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))));
if indic.target==1
    guess_trans(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))=sqrt(x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
else
    guess_trans(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))=log((params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))))./...
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))); 
end
guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=sqrt(x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))));
if indic.target==1
    guess_trans(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))));
end

lb=[];
ub=[];


%%
% Test Constraints and Objective Function %%%
f =  objective(guess_trans, T, params, list, Ftarget, indic);
[c, ceq] = constraints(guess_trans, T, params, init201519, list, Ems, indic);

% to examine stuff
ind=1:length(ceq);
ss=ind(abs(ceq)>1e-4);
tt=floor(ss/T); 
%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

objf=@(x)objective(x, T, params, list, Ftarget, indic);
constf=@(x)constraints(x, T, params, init201519, list, Ems, indic);

%if indic.target==0
options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-10,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%else
% options = optimset('algorithm','sqp','TolCon',1e-2,'Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
%options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%end

 [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
% without target the SP and OPTIMAL POL are the best
     save(sprintf('sqp_solu_targetOPT_505_spillover%d', indic.spillovers))

     % WHAT I did; 1) with LF as initial point started from active set with
     % TolCOn= 1e-6, then increased to 1e-8 with previous result, then with 1e-10 
 options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%end
% FINAL RUN WITH RESULT FROM ACTIVE SET NOW IN SQP
[xas,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
 save(sprintf('active_set_solu_notargetOPT_505_spillover%d_possible', indic.spillovers))
x=xas;
% if sum(x==xas)
%     fprintf('In version target=%d, the initial guess and the OPtimal pol are the same. With target=0 this is the LF and with target=1 it is the SP one.', indic.target);
% end

% transform
out_trans=exp(x);

out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
if indic.target == 0
    out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
else
    out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = (x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)).^2;
end
out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;

if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

% save results
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, wsgtil, S]= OPT_aux_vars(out_trans, list, params, T, init201519, indic);
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
 
opt_all=eval(symms.allvars);

if indic.target==1
    save(sprintf('OPT_target_active_set_0505_spillover%d', indic.spillovers), 'opt_all')
else
    save(sprintf('OPT_notarget_active_set_0505_spillover%d', indic.spillovers), 'opt_all')
end

end