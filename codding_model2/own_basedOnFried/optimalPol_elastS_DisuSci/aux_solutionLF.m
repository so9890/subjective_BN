function LF_t=aux_solutionLF(Sparams, SLF,pol, laggs, list, symms, indexx, params, indic)

% output
% LF_t: column vector of simulated variables in period t

% read in vars
gammalh=SLF.gammalh;
gammas=SLF.gammas;

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

Af=SLF.Af;
Ag=SLF.Ag;
An=SLF.An;

sff=SLF.sff;
sg=SLF.sg;
sn=SLF.sn;
S= SLF.S;

ws=SLF.ws;
pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;
%- params
sigmaa = Sparams.sigmaa;
chii = Sparams.chii;
sigmaas = Sparams.sigmaas;
chiis = Sparams.chiis;

eppse = Sparams.eppse;
eppsy = Sparams.eppsy;
deltay= Sparams.deltay;
thetaa = Sparams.thetaa;
thetan = Sparams.thetan;
thetag = Sparams.thetag;
thetaf = Sparams.thetaf;
rhof   = Sparams.rhof;
rhon   = Sparams.rhon;
rhog   = Sparams.rhog;
zh     = Sparams.zh; 
% lambdaa = pol(list.pol=='lambdaa');
lambdaa=SLF.lambdaa;
tauf = pol(list.pol=='tauf');
taul = pol(list.pol=='taul');
taus = pol(list.pol=='taus');
alphag=Sparams.alphag;
alphaf=Sparams.alphaf;
alphan=Sparams.alphan;
deltaa =Sparams.deltaa;
omegaa =Sparams.omegaa;

% auxiliary variables 

muu = C^(-thetaa);
E  = (SLF.F^((eppse-1)/eppse)+SLF.G^((eppse-1)/eppse))^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(SLF.pee./SLF.pn)^(eppsy).*E; % demand N final good producers 
Y = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));
xn=SLF.pn*Sparams.alphan*N;
xg=SLF.pg*Sparams.alphag*SLF.G;
xf=SLF.pf*(1-pol(list.pol=='tauf'))*Sparams.alphaf*SLF.F;
Cincome=Y-xn-xf-xg;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

muu     = C.^(-thetaa);
if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl*hl).^(1-taul))...
            +tauf.*pf.*F;
else
    SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
end

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

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
 Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);

 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
if abs(C-Cincome)>1e-10
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
if indic.noskill==0
    guess_trans=trans_guess(indexx('LF'), xx, params, list.params);
else
    guess_trans=trans_guess(indexx('LF_noskill'), xx, params, list.params);
end
    f=laissez_faire_nows(guess_trans, params, list, pol, laggs, indic);

if (max(abs(f)))>1e-8
    fprintf('f only solved at less than 1e-8')
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.allvars)';
end
 