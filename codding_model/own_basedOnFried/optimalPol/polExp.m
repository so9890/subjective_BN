function [tauf taul SWF] = polExp(pf, params, list)
% function to find necessary taul or tauf to meet emission target
% leaving the other one fixed at the optimal value absent emission target

syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn ...
    wh wl ws pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A real

symms.targprod=[pn pg];
list.targprod = string(symms.targprod);

if indic.tauffixed== 1
    symms.targlab=[hhn hhg hhf hh hl gammalh gammall wh wl taul];
else
    symms.targlab=[hhn hhg hhf hh hl gammalh gammall wh wl tauf];
end
list.targlab=string(symms.targlab); 

%- initial guess
pf=ones(T,1);
pn=ones(T,1);
pg=ones(T,1);

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
  AfLf, AgLg, AnLn]=resProdTarget(list, pn, pf, pg, params); 

%- Household and labour side

end

