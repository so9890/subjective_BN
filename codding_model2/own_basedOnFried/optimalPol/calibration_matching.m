function [x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
    = calibration_matching(MOM, symms, list, parsHelp, polhelp)
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
     sff sg sn  wsn wse real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sff sg sn wse wsn];
list.calib3 =string(symms.calib3);

% symms.calib4 = [Af_lag Ag_lag An_lag ...
%      sff sg sn];
% list.calib4 =string(symms.calib4);

%- all variables: to save base year variables!
syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn ...
    wh wl pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A...
     tauf taus taul lambdaa Ln Lg Lf SWF real
symms.allvars= [muu, hhf, hhg, hhn, hln, hlf, hlg, C, F, G, N, Y, E, Af, Ag, An, ...
    hl, hh,  sff, sg, sn, wh, wl, wse, wsn, pg, pn, pee, pf,  wlg, wln, wlf, xn, xg, xf, ...
    gammalh, gammall, SGov, Emnet, A, tauf, taus, taul, lambdaa, Ln, Lg, Lf, SWF ];
list.allvars  = string(symms.allvars);

%- variables and index for laissez faire
symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sff, sg, sn,...
             wh, wl, pg, pn, pee, pf, gammalh, gammall, lambdaa, wse, wsn];
list.choice  = string(symms.choice);

indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));
indexxLF.oneab = boolean(zeros(size(list.choice)));

indexxLF.lab(list.choice=='hl'| list.choice=='hh')=1;
indexxLF.exp(list.choice~='hl'& list.choice~='hh' & list.choice~='gammall'& list.choice~='gammalh' )=1;
indexxLF.sqr(list.choice=='gammall'| list.choice=='gammalh' )=1;

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
f = calibProd(x0, MOM, list, parsHelp);
prodf = @(x)calibProd(x,  MOM, list, parsHelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd] = fsolve(prodf, x0, options);
trProd=exp(solProd);
trProd(list.prod=='deltay')=1/(1+trProd(list.prod=='deltay'));

% required for next functions
[C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd(list, trProd, MOM, parsHelp, polhelp, 'calib');

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
Sparams.S=parsHelp(list.paramsdir=='S');

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
sg  = log(0.002);
wsn  = log(10);
wse  = log(10);
% log((params(listt=='S')-guess(indexx.sci))./guess(indexx.sci))
sff = log(0.004);
sn  = log(Sparams.S); %log((Sparams.S-Sparams.S*(1-1e-10))/(Sparams.S*(1-1e-10)));
% sg  = log((Sparams.S-0.002)/0.002);
% % log((params(listt=='S')-guess(indexx.sci))./guess(indexx.sci))
% sff = log((Sparams.S-0.004)/0.004);
% sn  = log(Sparams.S); %log((Sparams.S-Sparams.S*(1-1e-10))/(Sparams.S*(1-1e-10)));
    Af = AfLf/Lf;
    An = AnLn/Ln;
    Ag = AgLg/Lg;
Af_lag = log(Af);
An_lag = log(An);
Ag_lag = log(Ag);

phis = log(20);
% gammaa= parsHelp(list.paramsdir=='gammaa');
% etaa= parsHelp(list.paramsdir=='etaa');
% phii= parsHelp(list.paramsdir=='phii');
% alphan= parsHelp(list.paramsdir=='alphan');
% rhon= parsHelp(list.paramsdir=='rhon');
% rhof= parsHelp(list.paramsdir=='rhof');
% rhog= parsHelp(list.paramsdir=='rhog');
% A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

% gammasn =sqrt(exp(phis)- (gammaa*etaa*(A./An).^phii.*Sparams.S.^(etaa-1).*pn.*N.*(1-alphan).*An)./(rhon^etaa.*An));
% gammasf =sqrt(0);
% gammasg =sqrt(0);
x0 = eval(symms.calib3);

%%

%- test
% f= calibRem_nows(x0, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag);
% 
% % solving model
% modF3 = @(x)calibRem_nows(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag); 
% options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% 
% [x, fval, exitf] = fsolve(modF3, x0, options);
% % [x, fval, exitf] = fsolve(modF3, x, options);

lb=[];
ub=[];
objf=@(x)objectiveCALIBSCI(x);
constf=@(x)calibRem_nows(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag); 
options = optimset('algorithm','active-set','TolCon', 1e-12,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constf,options);

%%
Af_lag  = exp(x(list.calib3=='Af_lag'));
Ag_lag  = exp(x(list.calib3=='Ag_lag'));
An_lag  = exp(x(list.calib3=='An_lag'));
% phis    = exp(x(list.calib3=='phis'));
wsn      = exp(x(list.calib3=='wsn'));
wse      = exp(x(list.calib3=='wse'));

sff     = exp(x(list.calib3=='sff'));
sg      = exp(x(list.calib3=='sg'));
sn      = exp(x(list.calib3=='sn'));
% sff      = Sparams.S/(1+exp(x(list.calib3=='sff')));
% sg      = Sparams.S/(1+exp(x(list.calib3=='sg')));
% sn      = Sparams.S/(1+exp(x(list.calib3=='sn')));
% gammasg      = (x(list.calib3=='gammasg'))^2;
% gammasf      = (x(list.calib3=='gammasf'))^2;
% gammasn      = (x(list.calib3=='gammasn'))^2;

% now rerun with phis as given 

% x40=x(list.calib3~='phis');
% f= calibRem_nows_nophis(x40, MOM, list,trProd, parsHelp, polhelp, Af, An, Ag, phis);
% 
% % solving model
% modF3 = @(x)calibRem_nows_nophis(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, phis); 
% options = optimoptions('fsolve', 'TolFun', 10e-16, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% 
% [x, fval, exitf] = fsolve(modF3, x40, options);
% 
% 
% Af_lag  = exp(x(list.calib4=='Af_lag'));
% Ag_lag  = exp(x(list.calib4=='Ag_lag'));
% An_lag  = exp(x(list.calib4=='An_lag'));
% 
% sff      = Sparams.S/(1+exp(x(list.calib4=='sff')));
% sg      = Sparams.S/(1+exp(x(list.calib4=='sg')));
% sn      = Sparams.S/(1+exp(x(list.calib4=='sn')));
% gammasg      = (x(list.calib4=='gammasg'))^2;
% gammasf      = (x(list.calib4=='gammasf'))^2;
% gammasn      = (x(list.calib4=='gammasn'))^2;
% 
resSci=eval(symms.calib3); 

% get calibrated parameters and policy
[Sparams, Spol, params, pol]=parsSol(symms,trProd, trLab, resSci, parsHelp, list, polhelp);

%% save all results 
[x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519]=fsolution(symms, trProd, trLab, resSci, parsHelp, list, polhelp, MOM); 

%% Test if is calibration and baseline model solve LF in baseyear
guess_transLF=trans_guess(indexxLF, x0LF, params, list.params);
f=laissez_faire_nows(guess_transLF, params, list, pol, init201014);

if max(abs(f))>1e-10
    error('calibration is not a solution to LF')
else
    fprintf('Hurray!!! LF solves at baseline calibration!!!');
end
end