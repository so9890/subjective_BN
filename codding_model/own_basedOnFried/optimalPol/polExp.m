function [pn, pf, pg, F, G, N, Af, Ag, An, tauf] = polExp(pf, params, list, taul, T, Ems, indexx)


% function calculates exess supply of fossil as a function of pf

helper=load('SP_target.mat');
sp_t=helper.sp_all;
% polsym=sym('tauf',[T, length([tauf])]);

% in this 

syms tauf muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An HL HH sff sg sn lambdaa...
    wh wl ws pg pn pee gammalh gammall wlg wln wlf xn xg xf SGov Emnet A real

symms.targprod=[pn pg];
list.targprod = string(symms.targprod);
%indic.tauffixed=1;
symms.targlab=[hhn hhg hhf HH HL gammalh gammall wh wl lambdaa tauf];
list.targlab=string(symms.targlab); 

indexxEXPI.lab = boolean(zeros(size(list.targlab)));
indexxEXPI.exp = boolean(zeros(size(list.targlab)));
indexxEXPI.sqr = boolean(zeros(size(list.targlab)));
indexxEXPI.oneab = boolean(zeros(size(list.targlab)));

indexxEXPI.lab(list.targlab=='HL'| list.targlab=='HH')=1;
indexxEXPI.exp(list.targlab~='HL'& list.targlab~='HH' & list.targlab~='tauf'& list.targlab~='gammall'& list.targlab~='gammalh' )=1;
indexxEXPI.sqr(list.targlab=='gammall'| list.targlab=='gammalh')=1;
indexxEXPI.oneab(list.targlab=='tauf')=1;
 
indexx('targLab')=indexxEXPI;

%% - initial guess
pn=ones(T,1);
pg=ones(T,1);

x0= zeros(T*length(list.targprod),1);
x0((find(list.targprod=='pn')-1)*T+1:find(list.targprod=='pn')*T)= pn;
x0((find(list.targprod=='pg')-1)*T+1:find(list.targprod=='pg')*T)= pg;
guess_trans=log(x0);

f= target_prod(guess_trans, pf, Ems(1:T), params, list,T);

prodf = @(x)target_prod(x, pf, Ems(1:T), params, list,T);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd] = fsolve(prodf, guess_trans, options);
trProd=exp(solProd);

pn = trProd((find(list.targprod=='pn')-1)*T+1:find(list.targprod=='pn')*T);
pg = trProd((find(list.targprod=='pg')-1)*T+1:find(list.targprod=='pg')*T);

%- read in auxiliary variables following from pn,pf,pg
[Lnwln, Lgwlg, pf, F, pee, E, Y, N, G, xn, xg,  ...
            AgLg, AnLn]=resProdTarget(list, pn, pf, pg, params, Ems(1:T)');

%% - Household and labour side
x0= zeros(T*length(list.targlab),1);
x0((find(list.targlab=='hhn')-1)*T+1:find(list.targlab=='hhn')*T)=sp_t(1:T,list.allvars=='hhn');

x0((find(list.targlab=='hhg')-1)*T+1:find(list.targlab=='hhg')*T)=sp_t(1:T,list.allvars=='hhg');
x0((find(list.targlab=='hhf')-1)*T+1:find(list.targlab=='hhf')*T)=sp_t(1:T,list.allvars=='hhf'); 
x0((find(list.targlab=='gammalh')-1)*T+1:find(list.targlab=='gammalh')*T)= sp_t(1:T,list.allvars=='gammalh'); 
x0((find(list.targlab=='gammall')-1)*T+1:find(list.targlab=='gammall')*T)= sp_t(1:T,list.allvars=='gammall'); 
x0((find(list.targlab=='HL')-1)*T+1:find(list.targlab=='HL')*T)= sp_t(1:T,list.allvars=='hl'); 
x0((find(list.targlab=='HH')-1)*T+1:find(list.targlab=='HH')*T)= sp_t(1:T,list.allvars=='hh'); 
x0((find(list.targlab=='wh')-1)*T+1:find(list.targlab=='wh')*T)= (sp_t(1:T,list.allvars=='wh')); 
x0((find(list.targlab=='wl')-1)*T+1:find(list.targlab=='wl')*T)= (sp_t(1:T,list.allvars=='wl')); 
x0((find(list.targlab=='lambdaa')-1)*T+1:find(list.targlab=='lambdaa')*T)= (sp_t(1:T,list.allvars=='lambdaa')); 
x0((find(list.targlab=='tauf')-1)*T+1:find(list.targlab=='tauf')*T)= sp_t(1:T,list.allvars=='tauf'); 

guess_trans=trans_guess(indexx('targLab'), x0, params, list.params);

f = target_lab(x0, T, Lnwln, Lgwlg, xn, xg, AgLg, AnLn, params, list, taul,pf, F, Y);
labf = @(x)target_lab(x, T, Lnwln, Lgwlg, xn, xg, AgLg, AnLn, params, list, taul,pf, F, Y);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sollab, fval] = fsolve(labf, x0, options);
% 
helper=trans_allo_out(indexx('targLab'), sollab, params, list.params);
[ Ln, Lg, Lf, hln, hlg, hlf, xf, Ag, An, Af, A, wln, wlf, wlg, ...
        hhn, hhg, hhf, gammalh, gammall, hl, hh, wh, wl, tauf, lambdaa]=resLabTarget(helper,  list, pf, params,  Lnwln, Lgwlg, F, Y, AgLg, AnLn);


end

