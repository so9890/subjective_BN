%function [An0, Af0, Ag0, thetaf, thetan, thetag, lambdaa]= calibration_matching(MOM, symms, list, parsHelp, polhelp, targets)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear

% data
% Energy consumption share, GDP, Fossil to Green output share

% numeraire: Consumption Good

% balanced budget, 
% skills: match Consoli; or wage premia (then includes)


%% solve for skill 

eleh0   = 0.43;
thetag0 = 0.5; % bounded above and below
hhgHH0  = 0.1;
hlgHL0  = 0.1;
x0=[log((1-thetag0)/thetag0),log(eleh0), log((1-hhgHH0)/hhgHH0), log((1-hlgHL0)/hlgHL0)]; %log((1-zh0)/zh0)];
skillf = @(x)aux_calib_skill(x, MOM, parsHelp,list, polhelp);
options = optimoptions('fsolve', 'TolFun', 10e-16, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sol] = fsolve(skillf, x0, options);

thetag =1/(1+exp(sol(1))); 
eleh   = exp(sol(2));
hhgHH  = 1/(1+exp(sol(3))); % as targets for thetan and thetaf
hlgHL  = 1/(1+exp(sol(4)));

%% solve for intermediate and final good producers, prices

% x0=log(50); % guess pg
% prod = @(x)init_calib(x, MOM, parsHelp,list, polhelp);
% options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
% [sol, fval, exitf] = fsolve(prod, x0, options);
% pg=exp(sol);
% [pf, pe, pn, EY, NY, FY, GY, xgY, xnY, xfY]=aux_calib_Prod(MOM, pg, parsHelp,list, polhelp);
%% initial conditions

%- choiceCALIB variabes
syms hhf hhg hhn hln hlf hlg C F G Af Ag An ...
     hl hh sf sg sn wh wl ws pg pn pe pf gammalh gammall Y...
     Af0 Ag0 An0 thetan thetaf lambdaa real

symms.choiceCALIB = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, ...
    Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn,...
    pe, pf, gammalh, gammall, Y, Af0, Ag0, An0, thetan, thetaf, lambdaa];

list.choiceCALIB  = string(symms.choiceCALIB);

%initial values
Af0     = 1.877; % Fried 
Ag0     = 0.9196;
An0     = 1; 
% for init code
syms Ag_lag Af_lag An_lag real
symms.laggs = [Ag_lag, Af_lag, An_lag];
list.laggs= string(symms.laggs);
laggs=[Ag0, Af0, An0];

%- guesses parameters
thetan   = 0.5;
thetaf   = thetag*0.5;
lambdaa  = 1;
Af  = Af0*1.02; % Af
Ag  = Ag0*1.002; % Ag
An  = An0*1.02; % An
el  = 1; 
eh  = el/eleh; 

%- guesses variables
hhf =.02; % hhf
hhg =.04; % hhg
hhn =.04; % hhn
hln =hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf =hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg =hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 

hh     = (hhn+hhf+hhg)/(zl*el); % high skill market clearing
hl     = (hln+hlf+hlg)/(zh*eh); % low skill market clearing

Lg = hhg.^thetag.*hlg.^(1-thetag);
Ln = hhn.^thetan.*hln.^(1-thetan);
Lf = hhf.^thetaf.*hlf.^(1-thetaf);

% prices
pf=(1/MOM.FG)^(1/eppse)*pg;
pe=(pf^(1-eppse)+pg^(1-eppse))^(1/(1-eppse)); 
pn=((1-deltay^eppsy.*pe.^(1-eppsy))./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

% output
G = (Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));
F = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*(Af.*Lf); 
E = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N = (An.*Ln).*(pn.*alphan).^(alphan./(1-alphan));
Y = (deltay.*E.^((eppsy-1)/eppsy)+(1-deltay).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); 

% machines
xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

% consumption
C = zh*wh*eh*hh+zl*wl*el*hl;

% matching
%- unknowns: pg thetan thetaf Af Ag An el hhn hhf hhg
%1: 
f(1) = MOM.FG-F/G; %Af
f(2) = E*pe/Y - MOM.EpeY; % Ag
f(3) = An - 1; % An
f(4) = hhgHH - hhg/(zh*eh*hh); % as targets for thetan and thetaf
f(5) = hlgHL - hlg/(zl*el*hl);
f(6) = hl/hh - (MOM.whwl*eh/el)^((1-taul)/(taul+sigmaa)); % el
f(7) = Y+xn+xf+xg-C;  %=> pg
f(8) = (1-alphaf)*(1-tauf)*pf*F-(hhf)/thetaf; % labour demand
f(9) = (pn*N*(1-alphan))-hhn/thetan; 
f(10) = (pg*G*(1-alphag))-hhg/(thetag);

% sg  = 0.2;
% sf  = 0.3;
% sn  = 0.4;
gammalh = 0;
gammall = 0;

paramss = eval(symms.params);
poll    = eval(symms.pol);
%[x0,ni, checkk]=init(Af, An, Ag, hhg, hhf, hlg, hlf, hh, hl, F, G,  gammalh, gammall,paramss, list, poll, laggs, symms.choiceCALIB);


% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choiceCALIB)));
indexx.exp = boolean(zeros(size(list.choiceCALIB)));
indexx.sqr = boolean(zeros(size(list.choiceCALIB)));
indexx.oneab = boolean(zeros(size(list.choiceCALIB)));

indexx.lab(list.choiceCALIB=='hl'| list.choiceCALIB=='hh')=1;
indexx.exp(list.choiceCALIB~='hl'& list.choiceCALIB~='hh' & list.choiceCALIB~='gammall'&...
    list.choiceCALIB~='gammalh'& list.choiceCALIB~='thetan' & list.choiceCALIB~='thetag'& ...
    list.choiceCALIB~='thetaf' )=1;
indexx.sqr(list.choiceCALIB=='gammall'| list.choiceCALIB=='gammalh')=1;
indexx.oneab(list.choiceCALIB=='thetan'| list.choiceCALIB=='thetag'| list.choiceCALIB=='thetaf') = 1;

%-- transform
guess_trans=trans_guess(indexx, x0, paramss, list);

%- test
f=target_equ(guess_trans, MOM, paramss, list, poll);

%% - solving model

modFF = @(x)target_equ(x, MOM, paramss, list, poll);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

exitfl=0;
count=0;
countmax=100;
save=zeros(2,countmax);

while exitfl<=0 && count<countmax
    count=count+1;
    [sol, fval, exitf] = fsolve(modFF, guess_trans, options);
    guess_trans=sol;
    save(1,count)=max(fval);
    save(2,count)=exitf;
end


