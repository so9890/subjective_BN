function [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf,  taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill(x, list, params, T, init)

read_in_params;


Lg    = x((find(list.sp=='Lg')-1)*T+1:find(list.sp=='Lg')*T);
Ln    = x((find(list.sp=='Ln')-1)*T+1:find(list.sp=='Ln')*T);
%Lf     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
h     = x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T);

xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);
Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);
An     = x((find(list.sp=='An')-1)*T+1:find(list.sp=='An')*T);

C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);
sff     = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);
sg      = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);
sn      = x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T);


% initial values
An0=init(list.init=='An0');
Ag0=init(list.init=='Ag0');
Af0=init(list.init=='Af0');

% aux variables


%A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0;Ag(1:T-1)];
An_lag  = [An0;An(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);
S       = sff+sg +sn;
% the absolute amount of scientists supplied is determined by scientist
% preferences, the planner has to take this as given (otherwise there is no bound on science!)
% alternatively add disutility from scientists to objective function 
 
N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 

E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
Lf      =  (F./xf.^alphaf).^(1/(1-alphaf))./Af;

% prices compatible with sp solution 
pg = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag;
pn = (N./(An.*Ln)).^((1-alphan)/alphan)./alphan;
pf = (G./F).^(1/eppse).*pg; 
tauf = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf.*pf); 
% tauf2= 1-xf./(pf.*alphaf.*F); 
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
% S       = (wsn/chiis)^(1/sigmaas);
wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus
taus    = 1-wsgtil./wsn;
w      = pg.*(1-alphag).*G./Lg;
muu = C.^(-thetaa);


% finding lambdaa and taul
taul0 = 0.2*ones(size(sn));
lambdaa0=ones(size(sn));
ff=@(x)[w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)+tauf.*pf.*F;% balanced budget gov.
        chii*h.^(sigmaa+x(1:T))-(muu.*x(T+1:2*T).*(1-x(1:T)).*(w).^(1-x(1:T)))];
optionsfs = optimoptions('fsolve', 'TolFun', 10e-12,'Display','none');% 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

soll = fsolve(ff, [taul0; lambdaa0], optionsfs); 
taul=soll(1:T);
lambdaa=soll(T+1:T*2);


SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

% growth rate consumption

gammac = gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii;

% utility
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end

Utillab = chii*(h.^(1+sigmaa))./(1+sigmaa);
Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
SWF = Utilcon-Utillab-Utilsci;

 
contUtil= Utilcon(T)/(1-betaa* (1+gammac)^(1-thetaa));
contUtillab =1/(1-betaa)*(chii*h(T).^(1+sigmaa)./(1+sigmaa));
contUtilsci = 1/(1-betaa)*(chiis*S(T).^(1+sigmaas)./(1+sigmaas));
PVcontUtil = contUtil-contUtillab-contUtilsci;

end