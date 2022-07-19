function [x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
    = calibration_matching_tech(MOM, symms, list, parsHelp, polhelp)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters
% to this end, the model is solved 
% 1) from an aggregate pproduction side
% 2) the household probelm and labour producing firms
% 3) research side. 

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear
% zh: rationalises wage premium together with thetag, thetaf, thetan
% deltay: weight on energy in production function 
% omegaa: emissions per fossil 

% output
% indexx: map of indices indicating variable transformation for code.
%         contains laissez-faire index and calibration index. 
%% initaliase stuff needed for diverse solution 

%- calibration research side
syms Af_lag Ag_lag An_lag ...
     sff sg sn  wsf wsg wsn gammaa chiis phii etaa real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sff sg sn wsf wsg wsn gammaa chiis phii etaa];
list.calib3 =string(symms.calib3);

%- all variables: to save base year variables!
syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn ...
    wh wl pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A...
     tauf taus taul lambdaa Ln Lg Lf SWF S gammasn gammasf gammasg real
symms.allvars= [muu, hhf, hhg, hhn, hln, hlf, hlg, C, F, G, N, Y, E, Af, Ag, An, ...
    hl, hh,  sff, sg, sn, wh, wl, wsf, wsg, wsn, pg, pn, pee, pf,  wlg, wln, wlf, xn, xg, xf, ...
    gammalh, gammall, SGov, Emnet, A, tauf, taus, taul, lambdaa, Ln, Lg, Lf, SWF, S ];
list.allvars  = string(symms.allvars);

%- variables and index for laissez faire
symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sff, sg, sn,...
             wh, wl, pg, pn, pee, pf, gammalh, gammall, lambdaa, wsf, wsg, wsn, gammasf, gammasn, gammasg];
list.choice  = string(symms.choice);

indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));
indexxLF.oneab = boolean(zeros(size(list.choice)));

indexxLF.lab(list.choice=='hl'| list.choice=='hh'| list.choice=='sff'| list.choice=='sn'| list.choice=='sg')=1;
indexxLF.exp(list.choice~='hl'& list.choice~='hh'& list.choice~='sff'& list.choice~='sg'& list.choice~='sn' & list.choice~='gammasf'& list.choice~='gammasg'& list.choice~='gammasn'&list.choice~='gammall'& list.choice~='gammalh')=1;
indexxLF.sqr(list.choice=='gammall'| list.choice=='gammalh'| list.choice=='gammasf' | list.choice=='gammasg' | list.choice=='gammasn' )=1;

%- calibration productivity
syms omegaa deltay real
symms.prod= [pn, pg, omegaa, deltay];
list.prod = string(symms.prod);


%- calibration labour firms and HH side
syms thetan thetaf thetag lambdaa chii zh real
symms.calib = [hhn hhg hhf hh hl gammalh gammall wh wl thetan thetaf thetag ...
     lambdaa zh chii];
 list.calib = string(symms.calib);
%-- indexx
indexxcalib.lab = boolean(zeros(size(list.calib)));
indexxcalib.exp = boolean(zeros(size(list.calib)));
indexxcalib.sqr = boolean(zeros(size(list.calib)));
indexxcalib.oneab = boolean(zeros(size(list.calib)));
indexxcalib.sci = boolean(zeros(size(list.choice)));

indexxcalib.lab(list.calib=='hl'| list.calib=='hh')=1;
indexxcalib.exp(list.calib~='thetan' & list.calib~='thetaf' & list.calib~='thetag'& list.calib~='zh' ...
    & list.calib~='hl' & list.calib~='hh'& list.calib~='gammalh' & list.calib~='gammall')=1;
indexxcalib.sqr(list.calib=='gammall'| list.calib=='gammalh')=1;
indexxcalib.oneab(list.calib=='thetag'| list.calib=='thetan'| ...
    list.calib=='thetaf' | list.calib=='zh') = 1;

% save indices to map 
indexx = containers.Map({'LF', 'calib'}, {indexxLF, indexxcalib});
 %% First calibration reduced model: ONLY producers' side

pn=log(1);
pg=log(1);
omegaa=log(30);
deltay =  log((1-0.4)/0.4);

x0=eval(symms.prod);

%- solve
f = calibProd_GE(x0, MOM, list, parsHelp);
prodf = @(x)calibProd_GE(x,  MOM, list, parsHelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,'Display', 'Iter', 'Algorithm', 'levenberg-marquardt');%, );%, );%,  );
[solProd, fval] = fsolve(prodf, x0, options);
trProd=exp(solProd);
trProd(list.prod=='deltay')=1/(1+trProd(list.prod=='deltay'));

% required for next functions
[C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd_GE(list, trProd, MOM, parsHelp, polhelp, 'calib');

%% Labour side

%- initial guess
hhn = 0.02;
hhg = 0.01;
hhf = 0.02;
zh=0.2;
wl=1;
wh=MOM.whwl*wl;
thetan=0.4;
thetaf=0.4;
thetag=0.5;
    hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
    hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
    hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
hl = (hln+hlf+hlg)/((1-zh));
hh = (hhn+hhf+hhg)/zh;
gammall = 0;
gammalh = 0;
chii =10;
lambdaa=1;

x0=eval(symms.calib);

%- transform
guess_trans=trans_guess(indexxcalib, x0, parsHelp, list.paramsdir);

%test:
f =calibLabour(guess_trans,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, F, parsHelp, list, polhelp);

%- solve
Labf = @(x)calibLabour(x,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, F, parsHelp, list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solLab, fval] = fsolve(Labf, guess_trans, options);
trLab=trans_allo_out(indexxcalib, solLab, parsHelp, list.paramsdir);


%% Research side
%- from previous
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
Sparams=cell2struct(num2cell(trLab), cell_par, 2);

hhn = trLab(list.calib=='hhn');
hhg = trLab(list.calib=='hhg');
hhf = trLab(list.calib=='hhf');
hln = hhn*(1-Sparams.thetan)/(Sparams.thetan)*MOM.whwl; % hln
hlf = hhf*(1-Sparams.thetaf)/(Sparams.thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-Sparams.thetag)/(Sparams.thetag)*MOM.whwl; % hlg 

Lg = hhg.^Sparams.thetag.*hlg.^(1-Sparams.thetag);
Ln = hhn.^Sparams.thetan.*hln.^(1-Sparams.thetan);
Lf = hhf.^Sparams.thetaf.*hlf.^(1-Sparams.thetaf);

% initial guess
% ss=load('calib_initRes');
%%
 sg  = log(0.33);
 wsf  = log(3);
 wsg  = log(3);
 wsn  = log(3);
 gammaa =log(1);
 sff = log(0.33);
LB_etaa = 0.9999; 
  etaa = log(1.2-LB_etaa); %log(1.3);
 phii =  log((1-0.4)/0.4);
sn=sg;
%  sn  = log(MOM.targethour-exp(sg)-exp(sff)); 
    Af = AfLf/Lf;
    An = AnLn/Ln;
    Ag = AgLg/Lg;
 Af_lag = log(Af/(1.02));
 An_lag = log(An/(1.02));
 Ag_lag = log(Ag/(1.02));
 chiis =log(0.04);
% rhon=log(1);
% rhog=log(1);
% rhof=log(1);

x0 = eval(symms.calib3);

% x0=ss.x1;
%%

%- test
%  f= calibRem_nows_fsolve(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag);
% 
% % solving model
modF3 = @(x)calibRem_nows_fsolve(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, LB_etaa); 
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter','Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

% savex=x;
 sol3=load('calib_initResNEW1605_test.mat');
  [x, fval, exitf] = fsolve(modF3, sol3.x, options);
  [x, fval, exitf] = fsolve(modF3, x, options);
  save('calib_initResNEW1605', 'x');
  options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter'); %,'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

 [x2, fval, exitf] = fsolve(modF3, x, options);

%   lb=[];
%   ub=[];
%   objf=@(x)objectiveCALIBSCI(x);
% % % 
% % % 
%   handf=@(x)calibRem_nows(x, MOM, list, trProd, parsHelp, polhelp,  Af, An, Ag);
% 
%    options = optimset('algorithm','sqp','TolCon', 1e-10,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% % 
%   options = optimset('algorithm','active-set','TolCon', 1e-10,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%   [xsqp,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,handf,options);
 x=x2;
%%
 Af_lag  = exp(x(list.calib3=='Af_lag'));
  Ag_lag  = exp(x(list.calib3=='Ag_lag'));
  An_lag  = exp(x(list.calib3=='An_lag'));
etaa    = (exp(x(list.calib3=='etaa'))+LB_etaa);
chiis   = exp(x(list.calib3=='chiis'));
gammaa  = exp(x(list.calib3=='gammaa'));
phii    = 1/(1+exp(x(list.calib3=='phii')));
  sff     = exp(x(list.calib3=='sff'));
  sg      = exp(x(list.calib3=='sg'));
  sn      = exp(x(list.calib3=='sn'));
    wsf     = exp(x(list.calib3=='wsf'));
  wsg      = exp(x(list.calib3=='wsg'));
  wsn      = exp(x(list.calib3=='wsn'));
%% rhof     = exp(x(list.calib3=='rhof'));
% rhog      = exp(x(list.calib3=='rhog'));
% rhon      = exp(x(list.calib3=='rhon'));
 symms.calib4=symms.calib3(list.calib3~='chiis'&list.calib3~='etaa'&list.calib3~='gammaa'&list.calib3~='phii');
 list.calib4=string(symms.calib4);
 x0=eval(symms.calib4);
% % 
% % x02=x(list.calib3~='gammaa' & list.calib3~='chiis');
% % f=calibRem_nows_fixed(x02, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, gammaa, chiis);
% 
% sol4=load('initcalib4NEW');
% 
 modF4=@(x)calibRem_nows_fixed(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, gammaa, chiis, etaa, phii);
 options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e4, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
 [x4, fval, exitf] = fsolve(modF4, log(x0), options);
% 
% [x4, fval, exitf] = fsolve(modF4, sol4.x4, options);
% save('initcalib4NEW', 'x4')
%  %
  Af_lag  = exp(x4(list.calib4=='Af_lag'));
  Ag_lag  = exp(x4(list.calib4=='Ag_lag'));
  An_lag  = exp(x4(list.calib4=='An_lag'));
% chiis    = exp(x4(list.calib4=='chiis'));
  wsn       = exp(x4(list.calib4=='wsn'));
  wsg       = exp(x4(list.calib4=='wsg'));

  wsf       = exp(x4(list.calib4=='wsf'));

  sff     = exp(x4(list.calib4=='sff'));
  sg      = exp(x4(list.calib4=='sg'));
  sn      = exp(x4(list.calib4=='sn'));
%%
% syms chiis real
% symms.calib3=[symms.calib3, chiis]; 
% list.calib3=string(symms.calib3);
% chiis= ws/(C^parsHelp(list.paramsdir=='thetaa')*MOM.targethour^parsHelp(list.paramsdir=='sigmaas'));

resSci=eval(symms.calib3); 

% get calibrated parameters and policy
[Sparams, Spol, params, pol]=parsSol(symms,trProd, trLab, resSci, parsHelp, list, polhelp);

%% save all results 
[x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519]=fsolution(symms, trProd, trLab, resSci, parsHelp, list, polhelp, MOM); 

%% Test if is calibration and baseline model solve LF in baseyear
guess_transLF=trans_guess(indexxLF, x0LF, params, list.params);
indic.noskill=0;
f=laissez_faire_nows_sep(guess_transLF, params, list, pol, init201014, indic);

if max(abs(f))>1e-12
    error('calibration is not a solution to LF')
else
    fprintf('Hurray!!! LF solves at baseline calibration!!!');
end
end