function [c, ceq] = constraints(x, T, targets, params, list, Ems)
% function to read in constraints on government problem

% there is no order of constraints
%% read in stuff
% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% C, F, G : intermediate outputs
% Af, Ag, An : technology
% hl, hh

 hhf    = x(1:T);
 hhg    = x(  T+1:2*T);
 hlf    = x(2*T+1:3*T);
 hlg    = x(3*T+1:4*T);
 C      = x(4*T+1:5*T);
 F      = x(5*T+1:6*T);
 G      = x(6*T+1:7*T);
 Af     = x(7*T+1:8*T);
 Ag     = x(8*T+1:9*T);
 An     = x(9*T+1:10*T);
 hl     = x(10*T+1:11*T);
 hh     = x(11*T+1:12*T);
 
 
% parameters
thetaa = params(list.params=='thetaa');
sigmaa = params(list.params=='sigmaa');
% barHl  = params(list.params=='barHl');
% barHh  = params(list.params=='barHh');
zl  = params(list.params=='zl');
zh  = params(list.params=='zh');

%betaa = params(list.params=='betaa');
S      = params(list.params=='S');

thetaf = params(list.params=='thetaf');
thetan = params(list.params=='thetan');
thetag = params(list.params=='thetag');

alphag = params(list.params=='alphag');
alphaf = params(list.params=='alphaf');
alphan = params(list.params=='alphan');

eppsy = params(list.params=='eppsy');
eppse = params(list.params=='eppse');
deltay = params(list.params=='deltay');
gammaa = params(list.params=='gammaa');
etaa = params(list.params=='etaa');
rhof = params(list.params=='rhof');
rhog = params(list.params=='rhog');
rhon = params(list.params=='rhon');
phii = params(list.params=='phii');

Af0 = params(list.params=='Af0');
Ag0 = params(list.params=='Ag0');
An0 = params(list.params=='An0');

omegaa = targets(list.targets=='omegaa'); % carbon content of fossil energy
deltaa = targets(list.targets=='deltaa'); % natural sink
Ems    = Ems;   % vector of emission targets

%% auxiliary variables

hhn     = zh*hh-(hhf+hhg);
hln     = zl*hl-(hlf+hlg); 
hhhl    = hh./hl;
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
A       = max([Af', Ag', An'])'; 
Af_lag  = [Af0,Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0,Ag(1:T-1)];
An_lag  = [An0,An(1:T-1)];
A_lag   = [A0,A(1:T-1)];

mu      = C.^(-thetaa); % same equation in case thetaa == 1
% Muhh    = -zh*hh.^(-sigmaa);
% Muhl    = -zl*hl.^(-sigmaa);
% dIdhh   = lambdaa
% dIdhl 
sf      = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
sn      =  S-(sf+sg); 

E       = (F^((eppse-1)/eppse)+G^((eppse-1)/eppse)).^(eppse/(eppse-1));
% prices and policy elements
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); 
wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
ws      = (gammaa*etaa*(A_lag./Af_lag).^phii.*rhof^etaa.*sf.^(etaa-1).*pf.*F*(1-alphaf))./(Af./Af_lag); 


pe      = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay)^eppsy.*pe.^(1-eppsy)./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

%taus    = 1-((gammaa*etaa*(A_lag./Ag_lag).^phii*rhog^etaa.*sg^(etaa-1).*pg.*G.*(1-alphag))./(Ag./Ag_lag.*ws));
taul    = (exp(wh./wl)-sigmaa*(hhhl))./(exp(hhhl)+exp(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
lambdaa = hl.^(sigmaa+taul)./(mu.*(1-taul).*wl.^(1-taul));      % from FOC on low skill supply


% auxiliary stuff depending on prices
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +zl*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
        % subsidies and profits and wages scientists cancel
Y       = (pe.*E.^(1/eppse)/deltay).^(eppsy); % demand E final good producers 
N       = ((1-deltay)./pn.*Y.^(1/eppsy)).^eppse; % demand N final good producers 

wln     = pn.^(1/(1-alphan)).*(1-alphan)*alphan^(alphan/(1-alphan).*An); % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag)*alphag^(alphag/(1-alphag).*Ag);


xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
 c = zeros(2*T,1); % contains upper bounds for all direct optimisation 
                         % periods and 2 additional ones: 
                         %  IMP
 
                         % the constraints valid for T periods are:
                         % emission targets,
                         % resource constraint
                         
%%% 1. Emission constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1:T) = omegaa*F-(Ems-deltaa);


%%% 2. Implementability constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- create discount vector
% disc=repmat(betaa, 1,T);
% expp=0:T-1;
% vec_discount= disc.^expp;
% 
% c(2*T+1) =vec_discount*((mu.*C)-...
%     ((zl.*hl.^(sigmaa+1)+zh.*hh.^(sigmaa+1))./(1-taul)+mu.*SGov)); %

% rather without savings technology: one binding constraint per period
c(T+1:T*2)=(mu.*C)-((zl.*hl.^(sigmaa+1)+zh.*hh.^(sigmaa+1))./(1-taul)+mu.*SGov);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];
 

 ceq(1:T)       = pe.*E.^(1/eppse)-pg.*G.^(1/eppse); % green energy demand energy producers
 ceq(T+1:T*2)   = (deltay*E.^((eppsy-1)/eppsy)+(1-deltay)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
 ceq(T*2+1:T*3) = N-(An.*Ln).*(pn.*alphan)^(alphan./(1-alphan)); % from production function neutral good
 ceq(T*3+1:T*4) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
 ceq(T*4+1:T*5) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
 ceq(T*5+1:T*6) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
 ceq(T*6+1:T*7) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
 ceq(T*7+1:T*8) = ((1-alphan)*etaa*gammaa*An_lag.^(1-phii).*A_lag.^phii.*sn.^etaa.*pn.*N)-ws.*sn.*An*rhon^etaa; % wage scientists neutral
 ceq(T*8+1:T*9) = An_lag.*(1+(sn./rhon).^etaa.*(A_lag./An_lag).^phii)-An; % LOM neutral technology 
 ceq(T*9+1:T*10) = C+xf+xn+xg-Y; % final good market clearing

ceq = ceq';