function [LF_t, An, Ag,Af]=aux_solutionLF_xgrowth( SLF,pol, laggs, list, symms, indexx, params, indic, MOM,t,Emlim)

% output
% LF_t: column vector of simulated variables in period t

% read in variables
read_in_params;
read_in_pol;
% read in vars
gammalh=SLF.gammalh;
gammas=zeros(size(gammalh));

if indic.noskill==0
    hhf=SLF.hhf;
    hhg=SLF.hhg;
    hhn=SLF.hhn;
    hln=SLF.hln;
    hlg=SLF.hlg;
    hlf=SLF.hlf;
    hl=SLF.hl;
    hh=SLF.hh;
    gammall=SLF.gammall;
    wh=SLF.wh;
    wl=SLF.wl;
else
    
    Lg=SLF.Lg;
    Ln=SLF.Ln;
    Lf=SLF.Lf;
    w=SLF.w;
    h=SLF.h;
    wh=w; wl=w; gammall=gammalh; hh=h; hl=h; hhn=0;hhg=0; hhf=0; hln=0; hlg=0; hlf=0; 
end

F=SLF.F;
G=SLF.G;
C=SLF.C;
if indic.notaul ~=6
    lambdaa=SLF.lambdaa;
else
    taul=SLF.lambdaa;
end

if indic.limit_LF==1
    tauf=SLF.tauf;
end
%
%- initial condition
Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

% An=(1+vn)*An_lag;
% Ag=(1+vg)*Ag_lag;
% Af=(1+vf)*Af_lag;

    sn=  sn0;%  0.0034; %0.01*MOM.targethour;
    sff= sff0; % 5.5660e-10; %0.01*MOM.targethour;
    sg=  sg0; % 1.0305e-00.018; %0.01*MOM.targethour;

    S=(sn+sg+sff); %zeros(size(gammalh))
    gammasg=zeros(size(gammalh)); gammasn=zeros(size(gammalh)); gammasf=zeros(size(gammalh));
   
pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;

% auxiliary variables 
An=An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
Af=Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
Ag=Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);

if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf); 
else
    hh=h; hl=h; wh=w; wl=w; % this should suffice to have governmenta budget correct
end

if indic.notaul<2 || ...
   indic.notaul == 6 % tauf redistributed via income tax
    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +tauf.*F;
    Tls =zeros(size(F));    
    GovCon =zeros(size(F));
    
elseif indic.notaul == 2 ||...
        indic.notaul==3 %2,3,4,5,7
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon = tauf.*F; % GovCon = env tax consumed by government
    Tls =zeros(size(F)); 
    
elseif indic.notaul == 4 || indic.notaul ==5
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon =zeros(size(F));
    Tls  = tauf.*F;
    
elseif indic.notaul == 7 % earmarking
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon =zeros(size(F));
    Tls =zeros(size(F));    
    taus = tauf.*F./(pg.*G); % subsidy on green sector
end

E  = (SLF.F^((eppse-1)/eppse)+SLF.G^((eppse-1)/eppse))^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(SLF.pee./SLF.pn)^(eppsy).*E; % demand N final good producers 
Y  = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));
xn = SLF.pn*alphan*N;
xg = SLF.pg*(1+taus)*alphag*SLF.G;
xf = SLF.pf*alphaf*SLF.F;

wsf = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
wsg = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*(1+taus)*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*Ag);
wsn = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);
ws=wsn;
A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

muu      = C.^(-thetaa); % same equation in case thetaa == 1

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = (pg.*(1+taus)).^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 

Emnet     = omegaa*F-deltaa; % net emissions

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end
if indic.noskill==0
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
 Utilsci = zeros(size(gammalh));

 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
Cincome=Y-xn-xf-xg-GovCon- SGov;

if abs(C-Cincome)>1e-6
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
if indic.noskill==0

    guess_trans=trans_guess(indexx('LF_xgrowth'), xx, params, list.params);
    
else
    guess_trans=trans_guess(indexx('LF_noskill_xgrowth'), xx, params, list.params);
end
    f=laissez_faire_xgrowth(guess_trans, params, list, pol, laggs, indic, MOM,t,Emlim);

if (max(abs(f)))>1e-8
    fprintf('f only solved at less than 1e-8')
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.allvars)';
end
 