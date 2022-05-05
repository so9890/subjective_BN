function [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, wsgtil, S]= OPT_aux_vars(x, list, params, T, init201519, indic)

read_in_params;

hhf    = x((find(list.opt=='hhf')-1)*T+1:find(list.opt=='hhf')*T);
hhg    = x((find(list.opt=='hhg')-1)*T+1:(find(list.opt=='hhg'))*T);
hlf    = x((find(list.opt=='hlf')-1)*T+1:find(list.opt=='hlf')*T);
hlg    = x((find(list.opt=='hlg')-1)*T+1:find(list.opt=='hlg')*T);
C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
Af     = x((find(list.opt=='Af')-1)*T+1:find(list.opt=='Af')*T);
Ag     = x((find(list.opt=='Ag')-1)*T+1:find(list.opt=='Ag')*T);
An     = x((find(list.opt=='An')-1)*T+1:find(list.opt=='An')*T);
hl     = x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T);
hh     = x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T);  
S      = x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T);
sg     = x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T);    

 
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

tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); % production fossil
wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;

    
ws      = chiis*S.^sigmaas./muu; 
sff     = ((gammaa*etaa*(A_lag./Af_lag).^phii.*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(ws.*Af.*rhof^etaa)).^(1/(1-etaa));
%sg      = ((gammaa*etaa*(A_lag./Ag_lag).^phii.*pg.*G*(1-alphag).*Ag_lag)./(ws.*Ag.*rhog^etaa.*(1-taus))).^(1/(1-etaa));
sn      = ((gammaa*etaa*(A_lag./An_lag).^phii.*pn.*N*(1-alphan).*An_lag)./(ws.*An.*rhon^etaa)).^(1/(1-etaa));
% sg      = S -(sff+sn);

wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

taus    = 1-wsgtil./ws; % since (wsgtilde/(1-taus)=ws)


% assuming interior solution households
taul    = (log(wh./wl)-sigmaa*log(hhhl))./(log(hhhl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul

% lambdaa so that gov budget is balanced
lambdaa = (zh*(wh.*hh)+(1-zh)*(wl.*hl)+tauf.*pf.*F)./...
            (zh*(wh.*hh).^(1-taul)+(1-zh)*(wl.*hl).^(1-taul)); 
        % subsidies, profits and wages scientists cancel

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
        
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

% wh2 = thetag*(hhg./hlg).^(thetag-1).*wlg;
% utility
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
 Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);

 SWF = Utilcon-Utillab-Utilsci;

end