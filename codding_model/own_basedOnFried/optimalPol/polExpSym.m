function [tauf, taul, SWF] = polExpSym(pf, params, pol, list)

% symbolic solution to problem with emission target
syms pf real
symms.polnoLam= symms.pol(list.pol~='lambdaa');
list.polnoLam= string(symms.polnoLam);
list.symall = [list.polnoLam, string(pf)];
vecs=sym('a',[T,length([list.symall])]);
    for s = [list.symall]  % loop over list entries
        vecs(:,[list.symall]==s)=sym(sprintf('%s%d',s),  [T,1]); 
    end

% in this 
syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn lambdaa...
    wh wl ws pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A real

symms.targprod=[pn pg];
list.targprod = string(symms.targprod);
%indic.tauffixed=1;
symms.targlab=[hhn hhg hhf hh hl gammalh gammall wh wl lambdaa];
list.targlab=string(symms.targlab); 

% if indic.tauffixed== 1
%     symms.targlab=[hhn hhg hhf hh hl gammalh gammall wh wl taul];
%     symms.expi=[symms.choice, taul];
% else
%     symms.targlab=[hhn hhg hhf hh hl gammalh gammall wh wl tauf];
%     symms.expi=[symms.choice(list.choice~='F'), tauf];
% end

% list.expi=string(symms.expi);
% 
% indexxLF=indexx('LF');
% indexxEXPI.lab=indexxLF.lab(list.choice~='F');
% indexxEXPI.sqr=indexxLF.sqr(list.choice~='F');
% indexxEXPI.oneab=indexxLF.oneab(list.choice~='F');
% indexxEXPI.exp=indexxLF.exp(list.choice~='F');
% 
% indexx('EXPI')=indexxEXPI;
%% - initial guess

x0= zeros(T*length(list.targprod),1);
x0((find(list.targprod=='pn')-1)*T+1:find(list.targprod=='pn')*T)= pn;
x0((find(list.targprod=='pg')-1)*T+1:find(list.targprod=='pg')*T)= pg;
guess_trans=log(x0);

f= target_prod(guess_trans, pf, Ems, params, list,T);

prodf = @(x)target_prod(x, pf, Ems, params, list,T);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd] = fsolve(prodf, guess_trans, options);
trProd=exp(solProd);

pn = trProd((find(list.targprod=='pn')-1)*T+1:find(list.targprod=='pn')*T);
pg = trProd((find(list.targprod=='pg')-1)*T+1:find(list.targprod=='pg')*T);

%- read in auxiliary variables following from pn,pf,pg
[C, Lnwln, Lgwlg, Lfwlf, pf, F, pee, E, Y, N, G, xn, xg, xf, ...
  AfLf, AgLg, AnLn]=resProdTarget(list, pn, pf, pg, params, tauf, Ems'); 

%% - Household and labour side
x0= zeros(T*length(list.targlab),1);
x0((find(list.targlab=='hhg')-1)*T+1:find(list.targlab=='hhg')*T)=log(sp_t(:,list.allvars=='hhg'));
x0((find(list.targlab=='hhf')-1)*T+1:find(list.targlab=='hhf')*T)= log(sp_t(:,list.allvars=='hhf')); 
x0((find(list.targlab=='gammalh')-1)*T+1:find(list.targlab=='gammalh')*T)= sqrt(sp_t(:,list.allvars=='gammalh')); 
x0((find(list.targlab=='gammall')-1)*T+1:find(list.targlab=='gammall')*T)= sqrt(sp_t(:,list.allvars=='gammall')); 
x0((find(list.targlab=='hl')-1)*T+1:find(list.targlab=='hl')*T)= log((params(list.params=='upbarH')-sp_t(:,list.allvars=='hl'))./sp_t(:,list.allvars=='hl')); 
x0((find(list.targlab=='hh')-1)*T+1:find(list.targlab=='hh')*T)= log((params(list.params=='upbarH')-sp_t(:,list.allvars=='hh'))./sp_t(:,list.allvars=='hh')); 
x0((find(list.targlab=='wh')-1)*T+1:find(list.targlab=='wh')*T)= exp(sp_t(:,list.allvars=='wh')); 
x0((find(list.targlab=='wl')-1)*T+1:find(list.targlab=='wl')*T)= exp(sp_t(:,list.allvars=='wl')); 
x0((find(list.targlab=='lambdaa')-1)*T+1:find(list.targlab=='lambdaa')*T)= exp(sp_t(:,list.allvars=='lambdaa')); 


f = target_lab(x0,T, C, Lnwln, Lgwlg, Lfwlf, AfLf, AgLg, AnLn, params, list, taul, pf, F, tauf);
labf = @(x)target_lab(x,T, C, Lnwln, Lgwlg, Lfwlf, AfLf, AgLg, AnLn, params, list, taul, pf, F, tauf);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sollab, fval] = fsolve(labf, x0, options);
trProd=exp(solProd);
end

