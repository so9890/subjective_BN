function LF_t=aux_solutionLF_LambdaFlex(Sparams, SLF,pol, laggs, list, symms, indexx, params)

% output
% LF_t: column vector of simulated variables in period t

% read in vars
hhf=SLF.hhf;
hhg=SLF.hhg;
hhn=SLF.hhn;
hln=SLF.hln;
hlg=SLF.hlg;
hlf=SLF.hlf;
F=SLF.F;
G=SLF.G;
C=SLF.C;
Af=SLF.Af;
Ag=SLF.Ag;
An=SLF.An;
hl=SLF.hl;
hh=SLF.hh;
sff=SLF.sff;
sg=SLF.sg;
sn=SLF.sn;
gammalh=SLF.gammalh;
gammall=SLF.gammall;
wh=SLF.wh;
wl=SLF.wl;
ws=SLF.ws;
pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;

%- params
sigmaa = Sparams.sigmaa;
chii = Sparams.chii;

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
lambdaa = SLF.lambdaa;

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
xf=SLF.pf*(1-tauf)*Sparams.alphaf*SLF.F;
Cincome=Y-xn-xf-xg;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
muu     = C.^(-thetaa);
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl*hl).^(1-taul))...
        +tauf.*pf.*F;
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
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
 SWF = Utilcon-Utillab;

% test market clearing
if abs(C-Cincome)>1e-10
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
guess_trans=trans_guess(indexx, xx, params, list.params);
f=laissez_faire_LambdaFlex(guess_trans, params, list, pol, laggs);

if (max(abs(f)))>1e-8
    fprintf('f only solved at less than 1e-8')
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.allvars)';
end
 
 