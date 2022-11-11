function [ xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag, sff, sg,  ...
            F, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
             pg, pf, pee, w, wsf, wsg, tauf, taul, lambdaa,...
             wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill_noneutral(x, list, params, T, init, indic)

read_in_params;


Lg    = x((find(list.sp=='Lg')-1)*T+1:find(list.sp=='Lg')*T);
h     = x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);

if indic.xgrowth==0
    Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
    Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);
    sff     = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);
    sg      = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);
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
    Ag0=init(list.init=='Ag0');
    Af0=init(list.init=='Af0');

    % aux variables

    %A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
    Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
    Ag_lag  = [Ag0;Ag(1:T-1)];
else
    % initial values
    Ag_lag=init(list.init=='Ag0');
    Af_lag=init(list.init=='Af0');

    Af=zeros(size(xf));
    Ag=zeros(size(xf));

    for i=1:T
        Ag(i)=(1+vg)*Ag_lag;
        Af(i)=(1+vf)*Af_lag;
        %- update laggs
        Af_lag=Af(i);
        Ag_lag=Ag(i);
    end
    sff     = zeros(size(xf));
    sg      = zeros(size(xf));
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

A_lag   = (rhof*Af_lag+rhog*Ag_lag)./(rhof+rhog);
A       = (rhof*Af+rhog*Ag)/(rhof+rhog);
S       = sff+sg;   
Y       = E ;

% the absolute amount of scientists supplied is determined by scientist
% preferences, the planner has to take this as given (otherwise there is no bound on science!)
% alternatively add disutility from scientists to objective function 
 
% prices compatible with sp solution 
if indic.xgrowth==0
    wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
    wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus
else
    wsf     = zeros(size(xf));
    wsg     = zeros(size(xf));
end
w      = pg.*(1-alphag).*G./Lg;
muu = C.^(-thetaa);


% finding lambdaa and taul
if thetaa==1
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
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas);
 end
 SWF = Utilcon-Utillab-Utilsci;

 
contUtil= Utilcon(T)/(1-betaa* (1+gammac)^(1-thetaa));
contUtillab =1/(1-betaa)*(chii*h(T).^(1+sigmaa)./(1+sigmaa));
contUtilsci = 1/(1-betaa)*(chiis*S(T).^(1+sigmaas)./(1+sigmaas));
PVcontUtil = contUtil-contUtillab-contUtilsci;

end