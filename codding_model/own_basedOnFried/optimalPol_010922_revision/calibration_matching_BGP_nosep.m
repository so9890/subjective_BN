function [x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519, Sparams, Spol, params, pol, symms, MOM, indexx, list]...
    = calibration_matching_BGP_nosep(MOM, symms, list, parsHelp, polhelp, indic)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters
% to this end, the  economy is assumed to be on a balanced growth path with 
% research happening in all sectors. 

% THIS code: no separate markets for researchers! 

% To be chosen: 
% An0, Af0, Ag0 (2019)
AgAn = 0.9;
AgAf = 0.4;

% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear
% zh: rationalises wage premium together with thetag, thetaf, thetan
% deltay: weight on energy in production function 

% output
% indexx: map of indices indicating variable transformation for code.
%         contains laissez-faire index and calibration index. 

% 
if polhelp(list.poldir=='tauf')~=0 || polhelp(list.poldir=='taus')~=0 
error('tauf not 0 or taus not =0; have to change calibration and pay attention how env revs are recycled!')
end
%% initaliase stuff needed for diverse solution 

%- calibration research side
syms Af_lag Ag_lag An_lag ...
     sff sg sn  wsf wsg wsn ws gammaa chiis real
symms.calib3 = [Af_lag Ag_lag An_lag ...
     sff sg sn wsf wsg wsn chiis];
 
list.calib3 =string(symms.calib3);

%- all variables: to save base year variables!
syms muu chii hhf hhg hhn hln hlf hlg C F G N Y E Af Ag An hl hh sff sg sn ...
    wh wl pg pn pee pf gammalh gammall wlg wln wlf xn xg xf SGov Emnet A Ch Cl muuh muul...
     tauf taus taul lambdaa Ln Lg Lf SWF S gammasn gammasf gammasg  real
symms.allvars= [muu, hhf, hhg, hhn, hln, hlf, hlg, C, F, G, N, Y, E, Af, Ag, An, ...
    hl, hh,  sff, sg, sn, wh, wl, pg, pn, pee, pf,  wlg, wln, wlf, xn, xg, xf, ...
    gammalh, gammall, SGov, Emnet, A, tauf, taus, taul, lambdaa, Ln, Lg, Lf, SWF, S,...
    wsg, wsn, wsf, gammasg, gammasf, gammasn ];
list.allvars  = string(symms.allvars);


symms.allvars_ineq=[muuh, muul, Ch, Cl,symms.allvars];
list.allvars_ineq=string(symms.allvars_ineq);

%- variables and index for laissez faire
symms.choice = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, Af, Ag, An, hl, hh,  sff, sg, sn,...
             wh, wl, pg, pn, pee, pf, gammalh, gammall, lambdaa, wsg, wsn, wsf, gammasg, gammasf, gammasn];
list.choice  = string(symms.choice);

indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));
indexxLF.oneab = boolean(zeros(size(list.choice)));

indexxLF.lab(list.choice=='hl'| list.choice=='hh'| list.choice=='sff'| list.choice=='sg'| list.choice=='sn')=1;
indexxLF.exp(list.choice~='hl'& list.choice~='hh' & list.choice~='gammasg'& list.choice~='gammasf'& list.choice~='gammasn'&...
    list.choice~='gammall'& list.choice~='gammalh'& list.choice~='sff'& list.choice~='sn'& list.choice~='sg' & list.choice~='lambdaa')=1;
indexxLF.sqr(list.choice=='gammall'| list.choice=='gammalh'| list.choice=='gammasg'| list.choice=='gammasf'| list.choice=='gammasn' )=1;

%- variables for version with inequality
symms.choice_ineq=[symms.choice(list.choice~='C'), Ch, Cl];
list.choice_ineq =string(symms.choice_ineq);

indexxLF_ineq.lab = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.exp = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.sqr = boolean(zeros(size(list.choice_ineq)));
indexxLF_ineq.oneab = boolean(zeros(size(list.choice_ineq)));

indexxLF_ineq.lab(list.choice_ineq=='hl'| list.choice_ineq=='hh'| list.choice_ineq=='S')=1;
indexxLF_ineq.exp(list.choice_ineq~='hl'& list.choice_ineq~='hh' & list.choice_ineq~='gammas'&list.choice_ineq~='gammall'& list.choice_ineq~='gammalh'& list.choice_ineq~='S' )=1;
indexxLF_ineq.sqr(list.choice_ineq=='gammall'| list.choice_ineq=='gammalh'| list.choice_ineq=='gammas' )=1;

%- calibration productivity
syms omegaa deltay real
symms.prod= [pn, pg, deltay, omegaa];
list.prod = string(symms.prod(symms.prod~='omegaa')); % omegaa is added later


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
indexx = containers.Map({'LF', 'calib', 'LF_ineq'}, {indexxLF, indexxcalib, indexxLF_ineq});
 %% First calibration reduced model: ONLY producers' side

pn=log(1);
pg=log(1);
deltay =  log((1-0.4)/0.4);

x0=eval(symms.prod(symms.prod~='omegaa'));

%- solve
f = calibProd(x0, MOM, list, parsHelp, polhelp);
prodf = @(x)calibProd(x,  MOM, list, parsHelp, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[solProd, fval] = fsolve(prodf, x0, options);
trProd=exp(solProd);
trProd(list.prod=='deltay')=1/(1+trProd(list.prod=='deltay'));

% required for next functions
[C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd(list, trProd, MOM, parsHelp, polhelp, 'calib');

%- save omegaa to trProd
trProd =eval(symms.prod);
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
trLab=trans_allo_out(indexxcalib, solLab, parsHelp, list.paramsdir, indic);


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

%-- initial values
if isfile('calib_308_rem_Agr_mu.mat')
    rr=load(sprintf('calib_308_rem_Agr_mu'), 'x');
else
     sg  = log(0.1);
     wsg  = log(3);
     wsf  = log(3);
     wsn  = log(3);
     gammaa =log(1);
     sff = log(0.1);

     sn  = log(MOM.targethour-exp(sg)-exp(sff)); 

     Af_lag = log(Af/(1.02));
     An_lag = log(An/(1.02));
     Ag_lag = log(Ag/(1.02));
     chiis =log(0.1);

    rr.x = eval(symms.calib3);
end
% input to function
Af = AfLf/Lf;
An = AnLn/Ln;
Ag = AgLg/Lg;


%-- find upper bound on growth rate
xup=[log(Af), log(Ag), log(An), log(0.1)];

[ceq]= calibRem_upbar_gammaa(xup, MOM, list, parsHelp,polhelp,  Af, An, Ag);
modg=@(x)calibRem_upbar_gammaa(x, MOM, list, parsHelp,polhelp,  Af, An, Ag);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter'); %,'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[x, fval, exitf] = fsolve(modg, xup, options);
gammaa=exp(x(4));

%%
%- transform gamma to have an upper bound
% hhel=exp(rr.x(list.calib3=='gammaa'));
 rr.x=rr.x(list.calib3~='gammaa'); %=log((gammaup-gammaup*0.6)/gammaup*0.6);
%-- solving model
 f=calibRem_nows_fsolve_GOOD_sep_nogam(rr.x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, gammaa); 
modF3 = @(x)calibRem_nows_fsolve_GOOD_sep_nogam(x, MOM, list, trProd, parsHelp, polhelp, Af, An, Ag, gammaa); 
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e4, 'MaxIter', 3e5, 'Display', 'Iter'); %,'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[x, fval, exitf] = fsolve(modF3, rr.x, options);

save(sprintf('calib_308_rem_nogam'), 'x');

%- read in results for final validation 
Af_lag  = exp(x(list.calib3=='Af_lag'));
Ag_lag  = exp(x(list.calib3=='Ag_lag'));
An_lag  = exp(x(list.calib3=='An_lag'));
wsg       = exp(x(list.calib3=='wsg'));
wsf       = exp(x(list.calib3=='wsf'));
wsn       = exp(x(list.calib3=='wsn'));
% gammaa   = exp(x(list.calib3=='gammaa'));
% gammaa = gammaup/(1+exp(x(list.calib3=='gammaa')));

chiis   = exp(x(list.calib3=='chiis'));
sff     = exp(x(list.calib3=='sff'));
sg      = exp(x(list.calib3=='sg'));
sn      = exp(x(list.calib3=='sn'));


resSci=[eval(symms.calib3 ), gammaa]; 
syms gammaa real 
symms.calib3=[symms.calib3, gammaa];
% get calibrated parameters and policy
[Sparams, Spol, params, pol]=parsSol_GOOD(symms,trProd, trLab, resSci, parsHelp, list, polhelp);

%% save all results 
[x0LF, SL, SP, SR, Sall, Sinit201014, init201014 , Sinit201519, init201519]=fsolution_GOOD(symms, trProd, trLab, resSci, parsHelp, list, polhelp, MOM); 

%% Test if is calibration and baseline model solve LF in baseyear
guess_transLF=trans_guess(indexxLF, x0LF, params, list.params);
indic.noskill=0;
indic.ineq=0;
f=laissez_faire_nows_sep(guess_transLF, params, list, pol, init201014, indic);

if max(abs(f))>1e-11
    error('calibration is not a solution to LF')
else
    fprintf('Hurray!!! LF solves at baseline calibration!!!');
end
end