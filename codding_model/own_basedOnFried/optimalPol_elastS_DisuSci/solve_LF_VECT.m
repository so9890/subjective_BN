function [LF_SIM]=solve_LF_VECT(T, list, params,symms, init201519, helper, indic)
%test OPT policy result without target in competitive equilibrium
read_in_params;

if indic.sep==1
    %- new set of choice variables
    if indic.ineq==0
        list.choice=list.sepchoice;
        symms.choice=symms.sepchoice;
        list.allvars=list.sepallvars;
        symms.allvars=symms.sepallvars;
    else
        list.choice=list.sepchoice_ineq;
        symms.choice=symms.sepchoice_ineq;
        list.allvars=list.sepallvars_ineq;
        symms.allvars=symms.sepallvars_ineq;
    end
end

if indic.xgrowth==1
    list.choice=list.choice_xgrowth;
    symms.choice=symms.choice_xgrowth;
end
%helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d.mat', indic.spillovers));
varrs=helper.LF_SIM;
y=log(varrs);
z=sqrt(varrs);

% create new list
syms HL HH H w Lf Lg Ln real
if indic.noskill==0
    symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'), HL, HH];
else
    if indic.noneutral==0
    symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'), H, w, Lf, Lg, Ln];
    else
            symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'...
                &list.choice~='An' &list.choice~='sn'&list.choice~='pee' &list.choice~='pn'&list.choice~='gammasn' &list.choice~='wsn'), H, w, Lf, Lg];
    end
end

list.test=string(symms.test);

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
if indic.ineq==0
    if indic.BN==0
        C=y(list.allvars=='C', :)';
    else
        C=log((B-varrs(list.allvars=='C', :))./(varrs(list.allvars=='C', :)))';
    end
else
    if indic.BN==0
        Ch=y(list.allvars=='Ch', :)';
        Cl=y(list.allvars=='Cl', :)';
    else
        Ch=log((Bh-varrs(list.allvars=='Ch', :))./(varrs(list.allvars=='Ch', :)))';
        Cl=log((Bl-varrs(list.allvars=='Cl', :))./(varrs(list.allvars=='Cl', :)))';
    end
end
F=y(list.allvars=='F', :)';
G=y(list.allvars=='G', :)';

%- research sector if endogenous growth
if indic.xgrowth==0
    Af=y(list.allvars=='Af', :)';
    Ag =y(list.allvars=='Ag', :)';
    sff =z(list.allvars=='sff', :)';
    sg =z(list.allvars=='sg', :)';
    sn =z(list.allvars=='sn', :)';
    An =y(list.allvars=='An', :)';

    if indic.sep==0
        gammas =sqrt(zeros(size(lambdaa)));
        S =z(list.allvars=='S', :)';
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
lambdaa=y(list.allvars=='lambdaa', :)';


x0=eval(symms.test);
x0=x0(:);

% f=laissez_faireVECT(sol, params, list, varrs, init201519,T);
% 
% if max(abs(f))>1e-9
%    % equations where results are off
%     ind=1:length(f);
%     pos=ind(abs(f)>1e-9);
%     eqset= floor(pos./T);
%     error('optimal policy does not solve laissez faire. For the relevant equations see eqset');
%     % to evaluate stuff
% end

lb=[];
ub=[];

objf=@(x)objectiveCALIBSCI(x);
    if indic.xgrowth==0
        if indic.sep==1
            constLF=@(x)laissez_faireVECT_sep_fmincon(x, params, list, varrs, init201519,T, indic);
        else
            constLF=@(x)laissez_faireVECT_fmincon(x, params, list, varrs, init201519,T, indic);
        end
    else
        constLF=@(x)laissez_faireVECT_xgrowth_fmincon(x, params, list, varrs, init201519, T, indic);
    end

options = optimset('algorithm','active-set','TolCon', 1e-7,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constLF,options);

% test solution to 
if indic.noneutral==0
    if indic.xgrowth==0
        if indic.sep==0
            f=laissez_faireVECT(x, params, list, varrs, init201519,T, indic);
        else
            f=laissez_faireVECT_sep(x, params, list, varrs, init201519,T, indic);
        end
    else
        f=laissez_faireVECT_xgrowth(x, params, list, varrs, init201519, T, indic);
    end
else
    f=laissez_faireVECT_sep_non(x, params, list, varrs, init201519, T, indic);
end
if max(abs(f))>1e-7
    error('LF function does not solve')
end

% save results
if indic.noneutral==0
    if indic.xgrowth==0
    if indic.sep==0
        LF_SIM=aux_solutionLF_VECT(x, list, symms, varrs, params, T, indic);
    else
        LF_SIM=aux_solutionLF_VECT_sep(x, list, symms, varrs, params, T, indic);
    end
    else
        LF_SIM=aux_solutionLF_VECT_xgrowth(x, list, symms,varrs, params, T , indic, init201519);
    end
else
    LF_SIM=aux_solutionLF_VECT_sep_non(x, list, symms,varrs, params, T , indic);

end
end
