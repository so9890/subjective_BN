function LF_t=aux_solutionLF_sep(Sparams, SLF,pol, laggs, list, symms, indexx, params, indic)

% output
% LF_t: column vector of simulated variables in period t

% read in vars
gammalh=SLF.gammalh;
gammasg=SLF.gammasg;
gammasn=SLF.gammasn;
gammasf=SLF.gammasf;

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
Af=SLF.Af;
Ag=SLF.Ag;
An=SLF.An;

sff=SLF.sff;
sg=SLF.sg;
sn=SLF.sn;
S= sff+sg+sn;

wsg=SLF.wsg;
wsn=SLF.wsn;
wsf=SLF.wsf;
pg=SLF.pg;
pn=SLF.pn;
pee=SLF.pee;
pf=SLF.pf;
%- params
sigmaa = Sparams.sigmaa;
chii = Sparams.chii;
sigmaas = Sparams.sigmaas;
chiis = Sparams.chiis;
if indic.BN_red==0
    B= Sparams.B;
    Bl= Sparams.Bl;
    Bh= Sparams.Bh;
else
    B= 0.75*Sparams.B;
    Bl= Sparams.Bl;
    Bh= 0.75*Sparams.Bh;
    
end
zetaa = Sparams.zetaa;

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

E  = (SLF.F^((eppse-1)/eppse)+SLF.G^((eppse-1)/eppse))^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(SLF.pee./SLF.pn)^(eppsy).*E; % demand N final good producers 
Y = (deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));
xn=SLF.pn*Sparams.alphan*N;
xg=SLF.pg*Sparams.alphag*SLF.G;
xf=SLF.pf*(1-pol(list.pol=='tauf'))*Sparams.alphaf*SLF.F;
Cincome=Y-xn-xf-xg;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

if indic.ineq==0
    if indic.BN==0
        muu      = C.^(-thetaa); % same equation in case thetaa == 1
    else
        muu      = -(C-B).^(zetaa-1);
    end
muhh=muu;
muul=muu;
else
    if indic.BN==0
        muuh      = Ch.^(-thetaa); % same equation in case thetaa == 1
        muul      = Cl.^(-thetaa); % same equation in case thetaa == 1

    else
        muuh      = -(Ch-Bh).^(zetaa-1);
        muul      = -(Cl-Bl).^(zetaa-1);
    end
muu=muuh;
end
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
if indic.ineq==0
    if indic.BN==0
        if thetaa~=1
            Utilcon = (C.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = log(C);
        end
    else
        Utilcon=-(C-B).^(zetaa)./(zetaa); 
    end
else
    if indic.BN==0
        if thetaa~=1
            Utilcon = zh.*(Ch.^(1-thetaa))./(1-thetaa)+(1-zh).*(Cl.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = zh.*log(Ch)+(1-zh).*log(Cl);
        end
    else
        Utilcon=zh.*(-(Ch-Bh).^(zetaa)./(zetaa))+(1-zh).*(-(Cl-Bl).^(zetaa)./(zetaa)); 
    end

end

if indic.noskill==0
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
 Utilsci = chiis*sn.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sff.^(1+sigmaas)./(1+sigmaas);

 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
if indic.ineq==0
    diff=C-Cincome;
else
    diff=zh.*Ch+(1-zh).*Cl-Cincome;
end
if max(abs(diff))>1e-10
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice);
if indic.noskill==0
if indic.ineq==0

    if indic.sep==1
        if indic.BN==0
           guess_trans=trans_guess(indexx('LF_sep'), xx, params, list.params);
        else
           guess_trans=trans_guess(indexx('LF_sep_BN'), xx, params, list.params);
        end
    else
       guess_trans=trans_guess(indexx('LF'), xx, params, list.params);
    end
else

    if indic.sep==1
        if indic.BN==0
           guess_trans=trans_guess(indexx('LF_sep_ineq'), xx, params, list.params);
        else
           guess_trans=trans_guess(indexx('LF_sep_BN_ineq'), xx, params, list.params);
        end
    else
       guess_trans=trans_guess(indexx('LF_ineq'), xx, params, list.params);
    end
end
else
    guess_trans=trans_guess(indexx(sprintf('LF_noskill_sep%d', indic.sep)), xx, params, list.params);
end
f=laissez_faire_nows_sep(guess_trans, params, list, pol, laggs, indic);

if (max(abs(f)))>1e-10
    fprintf('f only solved at less than 1e-10')
else
    fprintf('saved variables are correct!')
end

% save stuff
if indic.ineq==0
    LF_t= eval(symms.allvars)';
else
    LF_t= eval(symms.allvars_ineq)';
end
end
 