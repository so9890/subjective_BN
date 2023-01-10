function [LF_SIM]=solve_LF_VECT(T, list, params,symms, init201519, helper, indic, Ems, MOM)
%test OPT policy result without target in competitive equilibrium
read_in_params;


if indic.xgrowth==1
    list.choice=list.choice_xgrowth;
    symms.choice=symms.choice_xgrowth;
end
varrs=helper.LF_SIM; % used to get policy
y=log(varrs);
z=sqrt(varrs);

% create new list
syms HL HH H w Lf Lg Ln tauf real
    if indic.noskill==0
        symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'), HL, HH];
    else

        symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                    &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                    &list.choice~='wl'&list.choice~='gammall'), H, w, Lf, Lg, Ln];
    end
if indic.limit_LF==1
    symms.test=[symms.test, tauf];
end
list.test=string(symms.test);

if indic.sep==2
      symms.test= symms.test(list.test~='wsf'&list.test~='wsg'& list.test~='wsn'& list.test~='ws' & list.test~='S'...
                     &list.test~='gammasf'&list.test~='gammasg'& list.test~='gammasn'& list.test~='gammas');
      list.test=string(symms.test);
elseif indic.sep==3
    syms se real
      symms.test= [symms.test(list.test~='gammasf'), se];
      list.test=string(symms.test);
end

% read in results from social planner
if indic.noskill==0
    HL= log((params(list.params=='upbarH')-varrs(list.allvars=='hl', :))./(varrs(list.allvars=='hl', :)))';
    HH= log((params(list.params=='upbarH')-varrs(list.allvars=='hh', :))./(varrs(list.allvars=='hh', :)))';
    hhf=y(list.allvars=='hhf', :)';
    hhg=y(list.allvars=='hhg', :)';
    hlg=y(list.allvars=='hlg', :)';
    hlf=y(list.allvars=='hlf', :)';
    wh =y(list.allvars=='wh', :)';
    wl=y(list.allvars=='wl', :)';
    hhn=y(list.allvars=='hhn', :)'; 
    hln =y(list.allvars=='hln', :)';
    gammall =sqrt(varrs(list.allvars=='gammall', :))';

else
    H= log((params(list.params=='upbarH')-varrs(list.allvars=='hh', :))./(varrs(list.allvars=='hh', :)))';
    w=y(list.allvars=='wh', :)';
    Lf=y(list.allvars=='Lf', :)';
    Ln=y(list.allvars=='Ln', :)';
    Lg=y(list.allvars=='Lg', :)';
    
end
C=y(list.allvars=='C', :)';
 
F=y(list.allvars=='F', :)';
G=y(list.allvars=='G', :)';

%- research sector if endogenous growth
se=(varrs(list.allvars=='sff',:)'+varrs(list.allvars=='sg',:)');

if indic.xgrowth==0
    Af = y(list.allvars=='Af', :)';
    Ag = y(list.allvars=='Ag', :)';
    if indic.xgrowth==0 && indic.noskill==1
        sff =log((params(list.params=='upbarH')-varrs(list.allvars=='sff', :))./(varrs(list.allvars=='sff', :)))';
        sg =log((params(list.params=='upbarH')-varrs(list.allvars=='sg', :))./(varrs(list.allvars=='sg', :)))';
        sn =log((params(list.params=='upbarH')-varrs(list.allvars=='sn', :))./(varrs(list.allvars=='sn', :)))';
        se = log((params(list.params=='upbarH')-se)./se);
    else
        sff =z(list.allvars=='sff', :)';
        sg =z(list.allvars=='sg', :)';
        sn =z(list.allvars=='sn', :)';
        se=sqrt(se);
    end
    An =y(list.allvars=='An', :)';

    if indic.sep==0
        gammas =sqrt(zeros(size(An)));
        S =log((params(list.params=='upbarS')-varrs(list.allvars=='S', :))./(varrs(list.allvars=='S', :)))';
        ws=z(list.allvars=='ws', :)';
    else
        gammasg =z(list.allvars=='gammasg', :)';
        gammasn =z(list.allvars=='gammasn', :)';
        gammasf =z(list.allvars=='gammasf', :)';
        wsn=z(list.allvars=='wsn', :)';
        wsf=z(list.allvars=='wsf', :)';
        wsg=z(list.allvars=='wsg', :)';
    end
end
gammalh =z(list.allvars=='gammalh', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';

if indic.notaul~=6
    lambdaa=varrs(list.allvars=='lambdaa', :)';
else
    lambdaa=varrs(list.allvars=='lambdaa', :)';     % lambdaa is in fact taul
end

if indic.limit_LF==1
    tauf=varrs(list.allvars=='tauf',:)';
end

x0=eval(symms.test);
x0=x0(:);

f=laissez_faireVECT_sep_NoRed(x0, params, list, varrs, init201519,T, indic, Ems);
% 
 if max(abs(f))>1e-9
%    % equations where results are off
     ind=1:length(f);
     pos=ind(abs(f)>1e-9);
     eqset= floor(pos./T);
%     error('optimal policy does not solve laissez faire. For the relevant equations see eqset');
%     % to evaluate stuff
 end

modFF = @(x)laissez_faireVECT_sep_NoRed(x, params, list, varrs, init201519,T, indic, Ems);
options = optimoptions('fsolve', 'TolFun', 10e-6,'Display','iter', 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[x0, fval, exitf] = fsolve(modFF, x0, options);

% options = optimoptions('fsolve', 'TolFun', 10e-6,'Display','iter', 'MaxFunEvals',8e3, 'MaxIter', 3e5); %,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% [x0, fval, exitf] = fsolve(modFF, x0, options);

lb=[];
ub=[];

objf=@(x)objectiveCALIBSCI(x);
    if indic.xgrowth==0
%         if indic.sep>=1
            constLF=@(x)laissez_faireVECT_sep_fmincon(x, params, list, varrs, init201519,T, indic, Ems);
%         else
%             constLF=@(x)laissez_faireVECT_fmincon(x, params, list, varrs, init201519,T, indic);
%         end
    else
        constLF=@(x)laissez_faireVECT_xgrowth_fmincon(x, params, list, varrs, init201519, T, indic, MOM);
    end

    if (indic.xgrowth==0 && indic.noskill==1) 
        options = optimset('algorithm','sqp','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',5e10,'MaxIter',6e10,'Display','iter');
    else
        options = optimset('algorithm','active-set','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    end
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constLF,options);

%  count=0;
%  while exitflag==-2 && count<4
%      count=count+1
[x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constLF,options);
% end
% test solution to 
if indic.xgrowth==0
%     if indic.sep==0
%         f=laissez_faireVECT(x, params, list, varrs, init201519,T, indic);
%     else
        f=laissez_faireVECT_sep_NoRed(x, params, list, varrs, init201519,T, indic, Ems);
%     end
else
    f=laissez_faireVECT_xgrowth(x, params, list, varrs, init201519, T, indic);
end

if max(abs(f))>1e-7
    error('LF function does not solve')
end

% save results
if indic.xgrowth==0
% if indic.sep==0
%     LF_SIM=aux_solutionLF_VECT(x, list, symms, varrs, params, T, indic);
% else
    LF_SIM=aux_solutionLF_VECT_sep(x, list, symms, varrs, params, T, indic);
% end
else
    LF_SIM=aux_solutionLF_VECT_xgrowth(x, list, symms,varrs, params, T , indic, init201519);
end
end
