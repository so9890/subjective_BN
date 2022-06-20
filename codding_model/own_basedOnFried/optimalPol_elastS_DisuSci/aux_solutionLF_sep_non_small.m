function LF_t=aux_solutionLF_sep_non_small(Sparams, SLF,pol, laggs, list, symms, indexx, params, indic)

% output
% LF_t: column vector of simulated variables in period t
minn=indic.minn;
% read in vars
gammalh=SLF.gammalh;
F=SLF.F;
pf=SLF.pf;

if indic.xgrowth==0
    gammasg=SLF.gammasg;
    gammasf=SLF.gammasf;

    Af=SLF.Af;
    Ag=SLF.Ag;
    sff=SLF.sff;
    sg=SLF.sg;
    wsg=SLF.wsg;
    wsf=SLF.wsf;
    S= sff+sg;

else

    %- initial condition
    read_in_params;
    Af_lag=laggs(list.init=='Af0');
    Ag_lag=laggs(list.init=='Ag0');

    Ag=(1+vg)*Ag_lag;
    Af=(1+vf)*Af_lag;
    sff=zeros(size(gammalh));
    sg=zeros(size(gammalh));
    S=zeros(size(gammalh));
    ws=zeros(size(gammalh));
    wsg=ws; wsf=ws; 
    gammasg=ws;  gammasf=ws;
end
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
tauf = pol(list.pol=='tauf');
taul = pol(list.pol=='taul');
taus = pol(list.pol=='taus');
alphag=Sparams.alphag;
alphaf=Sparams.alphaf;
alphan=Sparams.alphan;
deltaa =Sparams.deltaa;
omegaa =Sparams.omegaa;

% auxiliary variables 
h       = ((1-taul)/chii)^(1/(1+sigmaa)); 
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil

lambdaa =(w.*h+tauf.*pf.*F)./(w.*h).^(1-taul); 
SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
         +tauf.*pf.*F; 
C       = lambdaa.*(w.*h).^(1-taul)+SGov;

Lf      = F./(((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af); 
pg      = (w/((1-alphag)*alphag^(alphag/(1-alphag))*Ag))^(1-alphag);
G       = F.*(pf./pg).^eppse; % green demand
Lg      = G./(Ag*(pg.*alphag).^(alphag./(1-alphag)));

muu   = C.^(-thetaa); % same equation in case thetaa == 1
E     = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
pee = ones(size(E));
Y   = E; %(deltay^(1/eppsy)*E^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N^((eppsy-1)/eppsy))^(eppsy/(eppsy-1));
xg  = pg*Sparams.alphag*G;
xf  = pf*(1-pol(list.pol=='tauf'))*Sparams.alphaf*F;
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 
%- non neutral placeholders
N       = zeros(size(E));
pn      = zeros(size(E));
xn      = zeros(size(E));
sn      = zeros(size(E));
gammall = zeros(size(E));
gammasn = zeros(size(E));
Ln      = zeros(size(E));
wln     = zeros(size(E));
hhn     = zeros(size(E));
hln     = zeros(size(E));
wsn     = zeros(size(E));
An      = zeros(size(E));

%- het skill stuff
hhf     = zeros(size(E));
hhg     = zeros(size(E));
hlf     = zeros(size(E));
hlg     = zeros(size(E));

hl     = h;
hh     = h;
wh      = w;
wl      = w;

Cincome=Y-xf-xg;

A     = (rhof*Af+rhog*Ag)/(rhof+rhog);
Emnet = omegaa*F-deltaa; % net emissions

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
if max(abs(diff))>1e-7
    fprintf('market clearing does not hold, diff=%f', diff)
else
    fprintf('goods market cleared!')
end

% test variables read in properly
xx=eval(symms.choice_small);

guess_trans=trans_guess(indexx('LF_noneutral_sep_noskill'), xx, params, list.params, minn);
f=laissez_faire_nows_sep_non_noskillSmall(guess_trans, params, list, pol, laggs, indic);

if (max(abs(f)))>1e-7
    fprintf('f only solved at less than 1e-7')
else
    fprintf('saved variables are correct!')
end

% save stuff
LF_t= eval(symms.sepallvars)';
end
 