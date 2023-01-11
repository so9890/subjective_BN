function [f]=test_LF_VECT(T, list,  params,symms, init201519, helper, indic)
%test OPT policy result without target in competitive equilibrium
read_in_params;

% if indic.xgrowth==1
%     list.choice=list.choice_xgrowth;
%     symms.choice=symms.choice_xgrowth;
% end

varrs=helper.LF_SIM;
y=log(varrs);
z=sqrt(varrs);

% create new list
syms HL HH H w Lf Lg Ln real
if indic.noskill==0
    symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'), HL, HH];
else
    symms.test= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'), H, w, Lf, Lg, Ln];
end
list.test=string(symms.test);

% read in results from social planner
if indic.noskill==0
    HL= log((params(list.params=='upbarH')-varrs(list.allvars=='hl', :))./(varrs(list.allvars=='hl', :)))';
    HH= log((params(list.params=='upbarH')-varrs(list.allvars=='hh', :))./(varrs(list.allvars=='hh', :)))';
    hhf=y(list.allvars=='hhf', :)';
    hhn=y(list.allvars=='hhn', :)'; 
    hhg=y(list.allvars=='hhg', :)';
    hln =y(list.allvars=='hln', :)';
    hlg=y(list.allvars=='hlg', :)';
    hlf=y(list.allvars=='hlf', :)';
    wh =y(list.allvars=='wh', :)';
    wl=y(list.allvars=='wl', :)';
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
% if indic.xgrowth==0
    Af =y(list.allvars=='Af', :)';
    Ag =y(list.allvars=='Ag', :)';
    An =y(list.allvars=='An', :)';
      if indic.xgrowth==0 && indic.noskill==1
        sff =log((params(list.params=='upbarH')-varrs(list.allvars=='sff', :))./(varrs(list.allvars=='sff', :)))';
        sg =log((params(list.params=='upbarH')-varrs(list.allvars=='sg', :))./(varrs(list.allvars=='sg', :)))';
        sn =log((params(list.params=='upbarH')-varrs(list.allvars=='sn', :))./(varrs(list.allvars=='sn', :)))';
    else
        sff =z(list.allvars=='sff', :)';
        sg =z(list.allvars=='sg', :)';
        sn =z(list.allvars=='sn', :)';
    end
    if indic.sep==0
        ws=z(list.allvars=='ws', :)';
        S =log((params(list.params=='upbarS')-varrs(list.allvars=='S', :))./(varrs(list.allvars=='S', :)))';
        gammas =zeros(size(S));
    else
        wsf=z(list.allvars=='wsf', :)';
        wsg=z(list.allvars=='wsg', :)';
        wsn=z(list.allvars=='wsn', :)';
        gammasg=z(list.allvars=='gammasg', :)';
        gammasn=z(list.allvars=='gammasn', :)';
        gammasf=z(list.allvars=='gammasf', :)';    
    end
%end

gammalh =z(list.allvars=='gammalh', :)';
pg=y(list.allvars=='pg', :)';
pn=y(list.allvars=='pn', :)';
pee=y(list.allvars=='pee', :)';
pf=y(list.allvars=='pf', :)';
lambdaa=varrs(list.allvars=='lambdaa', :)';


x0=eval(symms.test);
x0=x0(:);

% test solution to 
% if indic.xgrowth==0
    if indic.sep<=1
       f= laissez_faireVECT_sep_NoRed(x0, params, list, varrs, init201519, T, indic);
%        f=laissez_faireVECT_sep(x0, params, list, varrs, init201519,T, indic);
    else
        f=laissez_faireVECT(x0, params, list, varrs, init201519,T, indic);
    end
% else
%     f=laissez_faireVECT_xgrowth(x0, params, list, varrs, init201519,T, indic);
% end
% to examine stuff
if ((indic.notaul==0 || indic.notaul==4) && indic.noskill==1 && indic.xgrowth==0)|| indic.notaul>7
    comp=1e-7;
else
    comp=1e-8;
end
 if max(abs(f))>comp
     error('LF function does not solve at 1e-8')
 else
     fprintf('Solution solves LF problem')
 end

end
