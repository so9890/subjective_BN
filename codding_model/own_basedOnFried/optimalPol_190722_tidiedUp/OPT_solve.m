function [symms, list, opt_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;
 
% symbilic variables and lists
syms hhf hhg hlf hlg C F G Af Ag An hl hh S sg Lf Lg h real

if indic.taus ==1
    symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh S sg];    
else
    if indic.noskill==0
        symms.opt = [hhf hhg hlf hlg C F G Af Ag An hl hh S];
    else
        symms.opt = [Lf Lg C F G Af Ag An h S];
    end
end
list.opt  = string(symms.opt); 

nn= length(list.opt); % number of variables
%%% Initial Guess %%%
% uses social planner of LF allocation 
%%%%%%%%%%%%%%%%%%%%%
if indic.target==1
    if isfile(sprintf('OPT_target_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul))
     helper=load(sprintf('OPT_target_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul));
     opt_all=helper.opt_all;
    
    x0 = zeros(nn*T,1);
    if indic.noskill==0
        x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =opt_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =opt_all(:,list.allvars=='hhg'); % hhg
        x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =opt_all(:,list.allvars=='hlf'); % hlf
        x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =opt_all(:,list.allvars=='hlg'); % hlg 
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =opt_all(:,list.allvars=='hl');  % hl
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =opt_all(:,list.allvars=='hh');  % hh
    else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
       % x0(T*(find(list.opt=='lambdaa')-1)+1:T*(find(list.opt=='lambdaa'))) =opt_all(:,list.allvars=='lambdaa'); % hlf
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))) =opt_all(:,list.allvars=='hh'); % as starting value use hh
    end
    
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =opt_all(:,list.allvars=='F');
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
    x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =opt_all(:,list.allvars=='Af');  % Af
    x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =opt_all(:,list.allvars=='Ag');  % Ag
    x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =opt_all(:,list.allvars=='An');  % An
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =opt_all(:,list.allvars=='S'); 
    if indic.taus==1
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg'); 
    end
    else

     helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, Sparams.etaa));
     helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_gammac.mat', indic.spillovers, indic.noskill, indic.sep));

     sp_all=helper.sp_all;

     if ~isvarname('sp_all')
         error('did not load sp solution')
     end 

     x0 = zeros(nn*T,1);
     if indic.noskill==0
         x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =sp_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
         x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =sp_all(:,list.allvars=='hhg'); % hhg
         x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =sp_all(:,list.allvars=='hlf'); % hlf
         x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =sp_all(:,list.allvars=='hlg'); % hlg 
         x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =sp_all(:,list.allvars=='hl');  % hl
         x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =sp_all(:,list.allvars=='hh');  % hh
     else
         x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =sp_all(:,list.allvars=='Lf'); % hlf
         x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =sp_all(:,list.allvars=='Lg'); % hlg 
        % x0(T*(find(list.opt=='lambdaa')-1)+1:T*(find(list.opt=='lambdaa')))   =sp_all(:,list.allvars=='lambdaa');  % hl
         x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h')))   =sp_all(:,list.allvars=='hh');  % hh    
     end
     x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =sp_all(:,list.allvars=='C');   % C
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =sp_all(:,list.allvars=='F');
     x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =sp_all(:,list.allvars=='G');   % G
     x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =sp_all(:,list.allvars=='Af');  % Af
     x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =sp_all(:,list.allvars=='Ag');  % Ag
     x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =sp_all(:,list.allvars=='An');  % An
     x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =(0.999).*sp_all(:,list.allvars=='S'); 
     if indic.taus==1
         x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =sp_all(:,list.allvars=='sg'); 
     end
  end  
elseif indic.target==0
       
    if isfile(sprintf('OPT_notarget_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul))
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, 0));
         
    opt_all=helper.opt_all;

    x0 = zeros(nn*T,1);
    if indic.noskill==0
        x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =opt_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =opt_all(:,list.allvars=='hhg'); % hhg
        x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =opt_all(:,list.allvars=='hlf'); % hlf
        x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =opt_all(:,list.allvars=='hlg'); % hlg 
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =opt_all(:,list.allvars=='hl');  % hl
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =opt_all(:,list.allvars=='hh');  % hh
    else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
     %   x0(T*(find(list.opt=='lambdaa')-1)+1:T*(find(list.opt=='lambdaa'))) =opt_all(:,list.allvars=='lambdaa'); % hlf
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))) =opt_all(:,list.allvars=='hh'); % as starting value use hh
    end
    
    x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =opt_all(:,list.allvars=='F');
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
    x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =opt_all(:,list.allvars=='Af');  % Af
    x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =opt_all(:,list.allvars=='Ag');  % Ag
    x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =opt_all(:,list.allvars=='An');  % An
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =opt_all(:,list.allvars=='S'); 
    if indic.taus==1
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg'); 
    end
    
    else    

    %- use competitive equilibrium with policy (taus=0; tauf=0; taul=0)
if etaa ~=1 
helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d.mat', indic.spillovers, indic.noskill, indic.sep));
  LF_SIM=helper.LF_SIM';
else
    helper= load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')));
  LF_SIM=helper.LF_BAU';
end


     x0 = zeros(nn*T,1);

     if indic.noskill==0
         x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
         x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
         x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
         x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
         x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
         x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh

     else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =LF_SIM(list.allvars=='Lf',1:T); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =LF_SIM(list.allvars=='Lg',1:T); % hhg
     %   x0(T*(find(list.opt=='lambdaa')-1)+1:T*(find(list.opt=='lambdaa'))) =LF_SIM(list.allvars=='lambdaa',1:T);
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))) =LF_SIM(list.allvars=='hh',1:T);
     end
     
     x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =LF_SIM(list.allvars=='C',1:T);   % C
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =LF_SIM(list.allvars=='F',1:T);
     x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =LF_SIM(list.allvars=='G',1:T);   % G
     x0(T*(find(list.opt=='Af')-1)+1:T*(find(list.opt=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
     x0(T*(find(list.opt=='Ag')-1)+1:T*(find(list.opt=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
     x0(T*(find(list.opt=='An')-1)+1:T*(find(list.opt=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An

     x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))     =LF_SIM(list.allvars=='S',1:T);  % hh
     
     if indic.taus==1
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % hh
     end
  
    end
end

%%

% Transform to unbounded variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- most of variables bounded by zero below
guess_trans=log(x0);

% guess_trans(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus')))=x0(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus'))) ;
%- exceptions with upper bound; hl, hh, S, and F in case of target

if indic.noskill==0
    guess_trans(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))))./...
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))));
    guess_trans(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))))./...
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))));
else
    guess_trans(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))))./...
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))));
end
% initial value of scientists might be zero (in version with target; other transformation required)
if indic.target==1
    guess_trans(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))=sqrt(x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
else
    guess_trans(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))))./...
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
f =  objective(guess_trans, T, params, list, Ftarget, indic, init201519);
[c, ceq] = constraints(guess_trans, T, params, init201519, list, Ems, indic);

% to examine stuff
% ind=1:length(ceq);
% ss=ind(abs(ceq)>1e-4);
% tt=floor(ss/T); 

%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519);
constf=@(x)constraints(x, T, params, init201519, list, Ems, indic);

%if indic.target==0
%else
% options = optimset('algorithm','sqp','TolCon',1e-6,'Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
% options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%end
if indic.target==1
%         options = optimset('algorithm','sqp','TolCon',1e-6,'Tolfun',1e-10,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
        if indic.spillovers==1
            options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        else
            options = optimset('algorithm','sqp','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        end
        [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
%           ll=load(sprintf('as_solu_notargetOPT_1205_spillover%d_taus%d_noskill%d_notaul0', indic.spillovers, indic.taus, indic.noskill))
%         [xsqp,fval,exitflag,output,lambda] = fmincon(ss.objf,xsqp,[],[],[],[],ss.lb,ss.ub,ss.constf,options);
% 
% ss=load(sprintf('sqp_solu_notargetOPT_1005_spillover%d_taus%d_noskill%d.mat', indic.spillovers, indic.taus, indic.noskill))
%         [x,fval,exitflag,output,lambda] = fmincon(ss.objf,ss.x,[],[],[],[],ss.lb,ss.ub,ss.constf,options);
elseif indic.target==0
       options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);

    %    options = optimset('algorithm','active-set','TolCon',1e-6,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
       [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
%         options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        %[x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
        options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
%         ss=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus), 'x')
        if ~ismember(exitflag, [1,4,5])
            error('active set did not solve sufficiently well')
        end
end
%f =  objective(x, T, params, list, Ftarget, indic);

%% transform
% helper=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus))
% x=xsqp;
out_trans=exp(x);
if indic.noskill==0
    out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
    out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
else
    out_trans((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)=upbarH./(1+exp(x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
end

if indic.target == 0
    out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarH./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
else
    out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = (x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)).^2;
end
if indic.taus==1
    out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
end
if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

%% save results
if indic.taus==1 % with taus
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, wsgtil, S]= OPT_aux_vars(out_trans, list, params, T, init201519, indic);
else
    if indic.noskill==0
    [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S]= OPT_aux_vars_notaus(out_trans, list, params, T, init201519, indic);
     else
        [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom(out_trans, list, params, T, init201519, indic);
            wlf=w; wlg =w; wln=w; wh=w; wl=w; hh=h; hl=h; 
            hhf=zeros(size(pn));hhg=zeros(size(pn));hhn=zeros(size(pn));hlf=zeros(size(pn));hlg=zeros(size(pn));hln=zeros(size(pn));
    end
    taus = zeros(size(pn));

end
gammall = zeros(size(pn));
gammalh = zeros(size(pn));

%%
opt_all=eval(symms.allvars);

% test if opt solution is a laissez faire solution 
% function throws error if solution is not a solution to LF
helper.LF_SIM=opt_all';
test_LF_VECT(T, list,  params,symms, init201519, helper, indic);
%%
if indic.target==1
    save(sprintf('OPT_target_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d', indic.spillovers, indic.taus, indic.noskill, indic.notaul), 'opt_all', 'Sparams')
else
    save(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d', indic.spillovers, indic.taus, indic.noskill, indic.notaul), 'opt_all', 'Sparams')
end

end