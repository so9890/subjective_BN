%% read in stuff
% inputs: x is a symbolic vector of choice variables for T periods
%params=symms.params; % by passing a symbolic vector read_in_params gives symbolic pars
read_in_params;

% separate to derive auxiliary variables
 hhf    = x((find(list.optsym=='hhf')-1)*T+1:find(list.optsym=='hhf')*T);
 hhg    = x((find(list.optsym=='hhg')-1)*T+1:(find(list.optsym=='hhg'))*T);
 hlf    = x((find(list.optsym=='hlf')-1)*T+1:find(list.optsym=='hlf')*T);
 hlg    = x((find(list.optsym=='hlg')-1)*T+1:find(list.optsym=='hlg')*T);
 C      = x((find(list.optsym=='C')-1)*T+1:find(list.optsym=='C')*T);
 F      = x((find(list.optsym=='F')-1)*T+1:find(list.optsym=='F')*T);
 G      = x((find(list.optsym=='G')-1)*T+1:find(list.optsym=='G')*T);
 Af     = x((find(list.optsym=='Af')-1)*T+1:find(list.optsym=='Af')*T);
 Ag     = x((find(list.optsym=='Ag')-1)*T+1:find(list.optsym=='Ag')*T);
 An     = x((find(list.optsym=='An')-1)*T+1:find(list.optsym=='An')*T);
 hl     = x((find(list.optsym=='HL')-1)*T+1:find(list.optsym=='HL')*T);
 hh     = x((find(list.optsym=='HH')-1)*T+1:find(list.optsym=='HH')*T);
 sn      = x((find(list.optsym=='sn')-1)*T+1:find(list.optsym=='sn')*T);
 sff      = x((find(list.optsym=='sff')-1)*T+1:find(list.optsym=='sff')*T);
 sg      = x((find(list.optsym=='sg')-1)*T+1:find(list.optsym=='sg')*T);
 
 % lagrange multiplier
mu_IMP     = x((find(list.optsym=='mu_IMP')-1)*T+1:find(list.optsym=='mu_IMP')*T);
% mu_MarketS = x((find(list.optsym=='mu_MarketS')-1)*T+1:find(list.optsym=='mu_MarketS')*T);

mu_NProd   = x((find(list.optsym=='mu_NProd')-1)*T+1:find(list.optsym=='mu_NProd')*T);
mu_OPThhg  = x((find(list.optsym=='mu_OPThhg')-1)*T+1:find(list.optsym=='mu_OPThhg')*T);
mu_OPThhn  = x((find(list.optsym=='mu_OPThhn')-1)*T+1:find(list.optsym=='mu_OPThhn')*T);
mu_OPThlg  = x((find(list.optsym=='mu_OPThlg')-1)*T+1:find(list.optsym=='mu_OPThlg')*T);
mu_OPThln  = x((find(list.optsym=='mu_OPThln')-1)*T+1:find(list.optsym=='mu_OPThln')*T);
mu_OPTLg  = x((find(list.optsym=='mu_OPTLg')-1)*T+1:find(list.optsym=='mu_OPTLg')*T);
mu_OPTLn  = x((find(list.optsym=='mu_OPTLn')-1)*T+1:find(list.optsym=='mu_OPTLn')*T);
mu_OPThh  = x((find(list.optsym=='mu_OPThh')-1)*T+1:find(list.optsym=='mu_OPThh')*T);
mu_LOMAf  = x((find(list.optsym=='mu_LOMAf')-1)*T+1:find(list.optsym=='mu_LOMAf')*T);
mu_LOMAg  = x((find(list.optsym=='mu_LOMAg')-1)*T+1:find(list.optsym=='mu_LOMAg')*T);
mu_LOMAn  = x((find(list.optsym=='mu_LOMAn')-1)*T+1:find(list.optsym=='mu_LOMAn')*T);
mu_demSff  = x((find(list.optsym=='mu_demSff')-1)*T+1:find(list.optsym=='mu_demSff')*T);
mu_demSg  = x((find(list.optsym=='mu_demSg')-1)*T+1:find(list.optsym=='mu_demSg')*T);
mu_demSn  = x((find(list.optsym=='mu_demSn')-1)*T+1:find(list.optsym=='mu_demSn')*T);

if indic.target==1
    mu_target  = x((find(list.optsym=='mu_target')-1)*T+1:find(list.optsym=='mu_target')*T);
end    
KT_hl       = x((find(list.optsym=='KT_hl')-1)*T+1:find(list.optsym=='KT_hl')*T);
KT_hh       = x((find(list.optsym=='KT_hh')-1)*T+1:find(list.optsym=='KT_hh')*T);
KT_S       = x((find(list.optsym=='KT_S')-1)*T+1:find(list.optsym=='KT_S')*T);

% initial values: CALIBRATED dont change
An0=init201519(list.init=='An0');
Ag0=init201519(list.init=='Ag0');
Af0=init201519(list.init=='Af0');

%% auxiliary variables

hhn     = zh*hh-(hhf+hhg);
hln     = (1-zh)*hl-(hlf+hlg); 
hhhl    = hh./hl;
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
%A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0;Ag(1:T-1)];
An_lag  = [An0;An(1:T-1)];
%A_lag   = [max([Af0,Ag0,An0]);A(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);


muu = C.^(-thetaa); % same equation in case thetaa == 1
% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% wages and policy elements

tauf    = 1-((F./(Af.*Lf)).^((1-alphaf)/alphaf))./(alphaf*pf); % production fossil

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 


wh      = thetaf*Lf./hhf.*wlf; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*Lf./hlf.*wlf;

    
wsf      = (chiis*sff.^sigmaas); % assuming inner solution
wsg      = chiis*sg.^sigmaas; % assuming inner solution
wsn      = chiis*sn.^sigmaas; % assuming inner solution


% assuming interior solution households
if indic.notaul==0
    taul   = (log(wh./wl)-sigmaa*log(hhhl))./(log(hhhl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
else
    taul   = zeros(size(sn));
end
% lambdaa so that gov budget is balanced
lambdaa = (zh*(wh.*hh)+(1-zh)*(wl.*hl)+tauf.*pf.*F)./...
            (zh*(wh.*hh).^(1-taul)+(1-zh)*(wl.*hl).^(1-taul)); 
        % subsidies, profits and wages scientists cancel


xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
S= sn+sff+sg;        
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
%- vector of utilities
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
end
 SWF = Utilcon-Utillab-Utilsci;