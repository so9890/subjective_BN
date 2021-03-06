function [symms, list, sp_all]=SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems)

% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn Lg Ln Lf h real
if indic.noskill==0
    symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn];
else
    symms.sp = [Lf Lg Ln xn xf xg An Af Ag C h F sff sg sn];
end
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

%- use competitive equilibrium with policy (taus=0; tauf=0; taul=0)
  helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d.mat', indic.spillovers, indic.noskill));
  LF_SIM=helper.LF_SIM;


if indic.target==0
    x0 = zeros(nn*T,1);
    Ftarget = 0; % placeholder

    if ~isfile(sprintf('SP_notarget_active_set_0505_spillover%d_noskill%d.mat', indic.spillovers, indic.noskill))
        CONTINUE
        if indic.noskill==0
            x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
            x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
            x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =LF_SIM(list.allvars=='hhn',1:T); % hhg
            x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
            x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
            x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =LF_SIM(list.allvars=='hln',1:T); % hlg 
            x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
            x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
        else
            x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
            x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='hln'))) =LF_SIM(list.allvars=='hln',1:T); % hlg 
          %  x0(T*(find(list.sp=='Lf')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
            x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
        end
        
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =LF_SIM(list.allvars=='xf',1:T); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =LF_SIM(list.allvars=='xg',1:T); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =LF_SIM(list.allvars=='xn',1:T); % hlg 
        x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
        x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
        x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An

        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',1:T);  % C
        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =LF_SIM(list.allvars=='F',1:T);  % C
        x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =LF_SIM(list.allvars=='sn',1:T);  % C
        x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % C
        x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =LF_SIM(list.allvars=='sff',1:T);  % C
    else
        helper=load(sprintf('SP_notarget_active_set_0505_spillover%d_noskill%d.mat', indic.spillovers, indic.noskill));
        sp_all=helper.sp_all;
            x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
            x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
            x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =sp_all(1:T, list.allvars=='hhn'); % hhg
            x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
            x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
            x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =sp_all(1:T, list.allvars=='hln'); % hlg 
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg 
            x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
            x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
            x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
            x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
            x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
            x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
    end
elseif indic.target==1 
    
    Ftarget = (Ems+deltaa)/omegaa;
    x0 = zeros(nn*T,1);
    
    if ~isfile(sprintf('SP_target_active_set_0505_spillover%d_noskill%d.mat', indic.spillovers, indic.noskill))
        kappaa = Ftarget./LF_SIM(list.allvars=='F',1:T); % ratio of targeted F to non-emission
        kappaa = kappaa*(1-1e-10);
        
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
    else
        helper= load(sprintf('SP_target_active_set_0505_spillover%d_noskill%d', indic.spillovers, indic.noskill), 'sp_all');
        sp_all=helper.sp_all;
            x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
            x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
            x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =sp_all(1:T, list.allvars=='hhn'); % hhg
            x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
            x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
            x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =sp_all(1:T, list.allvars=='hln'); % hlg 
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg 
            x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
            x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
            x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
            x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
            x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
            x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution; do not take init which refers to 2010-2014!
% this version here under assumption of first best policy in initial policy
% that is: tauf=0, taul=0, taus=0, lambdaa to balage budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initOPT= init201519; % as calibrated under BAU policy

%%% Transform variables to unbounded vars => requires less constraints! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 guess_trans=log(x0);
 guess_trans(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))))./...
     x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))));
 guess_trans(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))))./...
     x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))));
 guess_trans(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff')))=sqrt(x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))));
 guess_trans(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))=sqrt(x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn'))));
 guess_trans(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))=sqrt(x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg'))));
if indic.target==1
    guess_trans(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))=log((Ftarget'-x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F'))))./...
     x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F'))));
end
lb=[];
ub=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- reload parameters! to be correct!
% [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx]=get_params( T, indic, lengthh);

f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic, initOPT);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic);

% ind=1:length(ceq);
% ss=ind(abs(ceq)>1e-9);
% whe=floor(ss/T);

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic);

%  options = optimoptions('Algorithm','sqp','TolStep',1e-10,'TolFun',1e-16,'MaxFunEvals',500000,'MaxIter',6200,'Display','Iter','MaxSQPIter',10000);

options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
options = optimset('algorithm','active-set','TolCon',1e-12,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);

% print results
% fprintf('\"solved\" with precision %f', output.constrviolation);


out_trans=exp(x);
out_trans((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)=upbarH./(1+exp(x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
out_trans((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)=upbarH./(1+exp(x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
out_trans((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)=(x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
out_trans((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)=(x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
out_trans((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)=(x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;

if indic.target==1
    out_trans((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)=Ftarget'./(1+exp(x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)));
end
% additional variables 
 [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsn, wsf,  taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF]= SP_aux_vars_2S(out_trans, list, params, T, init201519);
        ws=wsn;
gammall = zeros(size(pn));
gammalh = zeros(size(pn));

sp_all=eval(symms.allvars);


if indic.target==1
    save(sprintf('SP_target_active_set_0505_spillover%d_noskill%d', indic.spillovers, indic.noskill), 'sp_all')
else
    save(sprintf('SP_notarget_active_set_0505_spillover%d_noskill%d', indic.spillovers, indic.noskill), 'sp_all')
end
end