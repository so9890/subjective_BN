function LF_t=aux_solutionLF_sep( SLF,pol, laggs, list, symms, indexx, params, indic, Emlim, t)

% output
% LF_t: column vector of simulated variables in period t

%- params
read_in_params;
read_in_pol;

% read in vars
gammalh=SLF.gammalh;
if indic.sep~=2
    if indic.sep==1
        
        gammasg=SLF.gammasg;
        gammasn=SLF.gammasn;
        gammasf=SLF.gammasf;
    elseif indic.sep==3
        
        gammasg=SLF.gammasg;
        gammasn=SLF.gammasn;
        gammasf=0;
    elseif indic.sep==0
        gammas=SLF.gammas;
    end
elseif indic.sep==2
    gammas =0;
end

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
if indic.ineq==0
    C=SLF.C;
    Cl=C;Ch=C;
else
    Cl=SLF.Cl;
    Ch=SLF.Ch;
C=Ch;
end
if indic.limit_LF==1
    tauf=SLF.tauf;
end
Af=SLF.Af;
Ag=SLF.Ag;
An=SLF.An;

sff=SLF.sff;
sg=SLF.sg;
sn=SLF.sn;
if indic.sep==0
    S=SLF.S;
else
    S= (sff+sg+sn)./zs;
end

if indic.sep==3
    se=SLF.se;
end

pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;
if indic.notaul ~=6
    lambdaa=SLF.lambdaa;
else
    taul=SLF.lambdaa;
end

% auxiliary variables 


if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf); 
else
    hh=h; hl=h; wh=w; wl=w; % this should suffice to have governmenta budget correct
end
if indic.sep ~=2
    if indic.sep~=0
        wsg=SLF.wsg;
        wsn=SLF.wsn;
        wsf=SLF.wsf;
    else
        ws=SLF.ws;
    end
else
    ws     = wspar;
end

if indic.Sun~=2
    helpS=zeros(size(F)); % helper for gov revenues from taxing scientists
elseif indic.Sun==2
    helpS=zs*(ws.*S-lambdaa.*(ws.*S).^(1-taul)); % version with lump sum transfers to finance subsidies on machine producers
end
if indic.notaul<2 || ...
   indic.notaul == 6 % tauf redistributed via income tax

    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +helpS...
        +tauf.*F;
    Tls =zeros(size(F));    
    GovCon =zeros(size(F));
elseif indic.notaul == 2 ||...
        indic.notaul==3 %2,3,4,5,7
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
           +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
           +helpS;
    GovCon = tauf.*F; % GovCon = env tax consumed by government
    Tls =zeros(size(F)); 
elseif indic.notaul == 4 || indic.notaul ==5
    if indic.noskill==0
        SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +helpS;
    else
        SGov = w.*h-lambdaa.*(w.*h).^(1-taul) +helpS;
    end
    GovCon =zeros(size(F));
    Tls  = tauf.*F;
elseif indic.notaul >= 7 % earmarking
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
         +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
         +helpS;
    GovCon =zeros(size(F));
    Tls = zeros(size(F));
end
if indic.xgrowth==0
    Tlsall = Tls-zs.*ws.*S;
else
    Tlsall=Tls;
end

E  = (SLF.F^((eppse-1)/eppse)+SLF.G^((eppse-1)/eppse))^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(SLF.pee./SLF.pn)^(eppsy).*E; % demand N final good producers 
Y  = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));
xn = SLF.pn*alphan*N;
xg = SLF.pg*(1+taus)*alphag*SLF.G;
xf = SLF.pf*alphaf*SLF.F;

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
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+((1-zh)).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
if indic.sep~=0
    Utilsci = chiis*sn.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sff.^(1+sigmaas)./(1+sigmaas);
else
    Utilsci=zs.*chiis*S.^(1+sigmaas)./(1+sigmaas);
end
 SWF = Utilcon-Utillab-Utilsci;

% test market clearing

Cincome=Y-xn-xf-xg-GovCon- SGov;

diff=C-Cincome;

if max(abs(diff))>1e-6
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
if indic.sep<=1

    if indic.noskill==0
        if indic.ineq==0
               guess_trans=trans_guess(indexx('LF'), xx, params, list.params);
        else
               guess_trans=trans_guess(indexx('LF_ineq'), xx, params, list.params);
        end
    else
        guess_trans=trans_guess(indexx(sprintf('LF_noskill_sep%d', indic.sep)), xx, params, list.params);
    end
    
    f=laissez_faire_nows_sep(guess_trans, params, list, pol, laggs, indic, Emlim, t);

else
   guess_trans=trans_guess(indexx(sprintf('LF_sep%d',indic.sep)), xx, params, list.params);
   if indic.sep==2
       f=laissez_faire_nows_partialS(guess_trans, params, list, pol, laggs, indic, Emlim, t);
   else
       f=laissez_faire_nows_sepSe(guess_trans, params, list, pol, laggs, indic, Emlim, t);
   end
end


if (max(abs(f)))>1e-7
    fprintf('f only solved at less than 1e-7, max abs f = %f',max(abs(f)) )
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.allvars)';

end
 