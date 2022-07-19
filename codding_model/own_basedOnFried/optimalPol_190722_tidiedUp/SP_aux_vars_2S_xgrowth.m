function [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_xgrowth(x, list, params, T, init, indic)

read_in_params;

hhf    = x((find(list.sp=='hhf')-1)*T+1:find(list.sp=='hhf')*T);
hhg    = x((find(list.sp=='hhg')-1)*T+1:(find(list.sp=='hhg'))*T);
hhn    = x((find(list.sp=='hhn')-1)*T+1:(find(list.sp=='hhn'))*T);
hlf    = x((find(list.sp=='hlf')-1)*T+1:find(list.sp=='hlf')*T);
hlg    = x((find(list.sp=='hlg')-1)*T+1:find(list.sp=='hlg')*T);
hln    = x((find(list.sp=='hln')-1)*T+1:find(list.sp=='hln')*T);
hl     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
hh     = x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T);
    
xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);


C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
Ch=C; Cl=C;


F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);


% initial values
An_lag=init(list.init=='An0');
Ag_lag=init(list.init=='Ag0');
Af_lag=init(list.init=='Af0');

Af=zeros(size(xf));
Ag=zeros(size(xf));
An=zeros(size(xf));

for i=1:T
    An(i)=(1+vn)*An_lag;
    Ag(i)=(1+vg)*Ag_lag;
    Af(i)=(1+vf)*Af_lag;
    %- update laggs
    An_lag=An(i);
    Af_lag=Af(i);
    Ag_lag=Ag(i);
end

% aux variables

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
%A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
% Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
% Ag_lag  = [Ag0;Ag(1:T-1)];
% An_lag  = [An0;An(1:T-1)];
% %A_lag   = [max([Af0,Ag0,An0]);A(1:T-1)];
% A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);
% 
% sff     = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
% sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa); 
% sn      = ((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);
S       = zeros(size(xf));
sff     = zeros(size(xf));
sg      = zeros(size(xf));
sn      = zeros(size(xf));

% the absolute amount of scientists supplied is determined by scientist
% preferences, the planner has to take this as given (otherwise there is no bound on science!)
% alternatively add disutility from scientists to objective function 
 
N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 

E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% prices compatible with sp solution 
pg = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag;
pn = (N./(An.*Ln)).^((1-alphan)/alphan)./alphan;
pf = (G./F).^(1/eppse).*pg; 
tauf = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf.*pf); 
% tauf2= 1-xf./(pf.*alphaf.*F); 
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
wsf     = zeros(size(xf));
wsn     = zeros(size(xf));

% S       = (wsn/chiis)^(1/sigmaas);
wsg     = zeros(size(xf));

wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
taul    = (log(wh./wl)-sigmaa*log(hh./hl))./(log(hh./hl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
lambdaa = (zh*(wh.*hh)+(1-zh)*(wl.*hl)+tauf.*pf.*F)./...
            (zh*(wh.*hh).^(1-taul)+(1-zh)*(wl.*hl).^(1-taul)); 
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
muuh=muu; muul=muu;

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 


%- continuation value: assuming after period T constant growth rate
%  note that A and Af refer to the last direct period T so do not use lags
%  here!, Assumption that research input sff is constant after period T
gammac = gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii;

% utility

        if thetaa~=1
            Utilcon = (C.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = log(C);
        end
    


Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
Utilsci = zeros(size(xf));

 SWF = Utilcon-Utillab-Utilsci;
 
 contUtil= Utilcon(T)/(1-betaa* (1+gammac)^(1-thetaa));
 contUtillab =1/(1-betaa)*(chii.*(zh.*hh(T).^(1+sigmaa)+(1-zh).*hl(T).^(1+sigmaa))./(1+sigmaa));
 contUtilsci = 1/(1-betaa)*(chiis*S(T).^(1+sigmaas)./(1+sigmaas));
 
 PVcontUtil = contUtil-contUtillab-contUtilsci;
end