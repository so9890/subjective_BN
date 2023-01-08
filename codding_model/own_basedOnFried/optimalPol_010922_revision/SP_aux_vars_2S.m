function [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S(x, list, params, T, init, indic)

read_in_params;

hhf    = x((find(list.sp=='hhf')-1)*T+1:find(list.sp=='hhf')*T);
hhg    = x((find(list.sp=='hhg')-1)*T+1:(find(list.sp=='hhg'))*T);
hlf    = x((find(list.sp=='hlf')-1)*T+1:find(list.sp=='hlf')*T);
hlg    = x((find(list.sp=='hlg')-1)*T+1:find(list.sp=='hlg')*T);
hl     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
hh     = x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T);
    
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);
Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);

xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
An     = x((find(list.sp=='An')-1)*T+1:find(list.sp=='An')*T);
sn      = x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T);
hhn    = x((find(list.sp=='hhn')-1)*T+1:(find(list.sp=='hhn'))*T);
hln    = x((find(list.sp=='hln')-1)*T+1:find(list.sp=='hln')*T);



C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
Ch=C; Cl=C;
F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);
sff     = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);
sg      = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);


% initial values
An0=init(list.init=='An0');
Ag0=init(list.init=='Ag0');
Af0=init(list.init=='Af0');

% aux variables

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
%A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0;Ag(1:T-1)];
%A_lag   = [max([Af0,Ag0,An0]);A(1:T-1)];

% 
% sff     = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
% sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa); 
% sn      = ((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);
% the absolute amount of scientists supplied is determined by scientist
% preferences, the planner has to take this as given (otherwise there is no bound on science!)
% alternatively add disutility from scientists to objective function 
 
G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

Ln      = hhn.^thetan.*hln.^(1-thetan);
An_lag  = [An0;An(1:T-1)];
N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
pn      = (N./(An.*Ln)).^((1-alphan)/alphan)./alphan;
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);
S       = (sff+sg +sn);
A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector


% prices compatible with sp solution 
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf); % production fossil

tauf      = (G./F).^(1/eppse).*pg-pf; % optimality energy producers

pee     = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*Af_lag)./(Af.*rhof^etaa); 
wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        (pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        (pf).^(1/(1-alphaf)).*Af;
taul    = (log(wh./wl)-sigmaa*log(hh./hl))./(log(hh./hl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
lambdaa = (zh*(wh.*hh)+((1-zh))*(wl.*hl)+tauf.*F)./...
            (zh*(wh.*hh).^(1-taul)+((1-zh))*(wl.*hl).^(1-taul)); 
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*F;
Emnet     = omegaa*F-deltaa; % net emissions


        muu   = C.^(-thetaa); % same equation in case thetaa == 1

    muuh=muu; muul=muu;


wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 


% utility

if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end

Utillab = chii.*(zh.*hh.^(1+sigmaa)+((1-zh)).*hl.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
 end
 SWF = Utilcon-Utillab-Utilsci;
 
%% social welfare

    %- create discount vector
     disc=repmat(betaa, 1,T);
     expp=0:T-1;
     vec_discount= disc.^expp;
     PVSWF = vec_discount*SWF;
    
    % continuation value
    %- last period growth rate as proxy for future growth rates
      gammay = Y(T)/Y(T-1)-1;
    PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(T);
    PVwork = indic.PVwork*1/(1-betaa)*(Utillab(T)+Utilsci(T)); % this decreases last period work and science 
    PV= betaa^T*(PVconsump-PVwork);

    %Objective function value:
    %!! Dot product!!! so no dot.*
    % f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);+
    objF=(vec_discount*(SWF-indic.extern*weightext*(omegaa.*F).^extexpp)+indic.PV*PV);

end