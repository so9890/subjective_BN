function [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill(x, list, params, T, init, indic)

read_in_params;


Lg    = x((find(list.sp=='Lg')-1)*T+1:find(list.sp=='Lg')*T);
Ln    = x((find(list.sp=='Ln')-1)*T+1:find(list.sp=='Ln')*T);
%Lf     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
h     = x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T);

xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);

if indic.xgrowth==0
    Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
    Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);
    An     = x((find(list.sp=='An')-1)*T+1:find(list.sp=='An')*T);
    sff     = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);
    sg      = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);
    sn      = x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T);
end

if indic.ineq==0
    C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
    Ch=C; Cl=C;
else
    Ch      = x((find(list.sp=='Ch')-1)*T+1:find(list.sp=='Ch')*T);
    Cl      = x((find(list.sp=='Cl')-1)*T+1:find(list.sp=='Cl')*T);
    C=Ch;
end
F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);


if indic.xgrowth==0
    % initial values
    An0=init(list.init=='An0');
    Ag0=init(list.init=='Ag0');
    Af0=init(list.init=='Af0');

    % aux variables

    %A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
    Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
    Ag_lag  = [Ag0;Ag(1:T-1)];
    An_lag  = [An0;An(1:T-1)];
else
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
    sff     = zeros(size(xf));
    sg      = zeros(size(xf));
    sn      = zeros(size(xf));

end
G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Lf      =  (F./xf.^alphaf).^(1/(1-alphaf))./Af;
pg = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag;
pf = (G./F).^(1/eppse).*pg; 
tauf = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf.*pf); 
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

if indic.noneutral==0
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    An_lag  = [An0;An(1:T-1)];
    N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
    pn      = (N./(An.*Ln)).^((1-alphan)/alphan)./alphan;
    A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);
    S       = sff+sg +sn;
    A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
    Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
    wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
    wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector

else
    An      =zeros(size(E));
    An_lag  =zeros(size(E));
    N       =zeros(size(E));
    pn      =zeros(size(E));
    A_lag   = (rhof*Af_lag+rhog*Ag_lag)./(rhof+rhog);
    A       = (rhof*Af+rhog*Ag)/(rhof+rhog);
    S       = sff+sg;   
    Y       = E ;
    wsn     =zeros(size(E));
    wln     =zeros(size(E)); 
    sn      =zeros(size(E));
    xn      =zeros(size(E));
    Ln      =zeros(size(E));
end
% the absolute amount of scientists supplied is determined by scientist
% preferences, the planner has to take this as given (otherwise there is no bound on science!)
% alternatively add disutility from scientists to objective function 
 
% prices compatible with sp solution 
if indic.xgrowth==0
    wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
    if indic.noneutral==0
        wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
    end
    wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus
else
    wsf     = zeros(size(xf));
    wsn     = zeros(size(xf));
    wsg     = zeros(size(xf));
end
w      = pg.*(1-alphag).*G./Lg;
muu = C.^(-thetaa);


% finding lambdaa and taul
if thetaa~=1
    taul0 = 0.2*ones(size(E));
    lambdaa0=ones(size(E));
    ff=@(x)[w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)+tauf.*pf.*F;% balanced budget gov.
            chii*h.^(sigmaa+x(1:T))-(muu.*x(T+1:2*T).*(1-x(1:T)).*(w).^(1-x(1:T)))];
    optionsfs = optimoptions('fsolve', 'TolFun', 10e-12,'Display','none');% 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

    soll = fsolve(ff, [taul0; lambdaa0], optionsfs); 
    taul=soll(1:T);
    lambdaa=soll(T+1:T*2);
else
            taul    = 1-chii.*h.^(sigmaa+1);
    lambdaa = (w.*h+tauf.*pf.*F)./(w.*h).^(1-taul);
end


SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
Emnet     = omegaa*F-deltaa; % net emissions


% growth rate consumption

gammac = gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii;

% utility
if indic.ineq==0
    if thetaa~=1
     Utilcon = (C.^(1-thetaa))./(1-thetaa);
    elseif thetaa==1
     Utilcon = log(C);
    end
else
    error('not yet separate science and inequality coded')
end

Utillab = chii*(h.^(1+sigmaa))./(1+sigmaa);
 if indic.sep==0
        Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
 end
 SWF = Utilcon-Utillab-Utilsci;

 
contUtil= Utilcon(T)/(1-betaa* (1+gammac)^(1-thetaa));
contUtillab =1/(1-betaa)*(chii*h(T).^(1+sigmaa)./(1+sigmaa));
contUtilsci = 1/(1-betaa)*(chiis*S(T).^(1+sigmaas)./(1+sigmaas));
PVcontUtil = contUtil-contUtillab-contUtilsci;

end