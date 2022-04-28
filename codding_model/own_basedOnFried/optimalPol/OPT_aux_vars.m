function [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsn, wsg, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf]= OPT_aux_vars(x, list, params, T, init)

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

% initial values: CALIBRATED dont change
An0=init(list.init=='An0');
Ag0=init(list.init=='Ag0');
Af0=init(list.init=='Af0');

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


muu      = C.^(-thetaa); % same equation in case thetaa == 1

% scientists follow from LOW fossil and green 
sff     = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
sn      = ((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);

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
wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa);
wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

taus    = 1-wsgtil./wsf;
wsg     = wsgtil./(1-taus); % this ensures wsg=ws

% assuming interior solution households
taul    = (exp(wh./wl)-sigmaa*exp(hhhl))./(exp(hhhl)+exp(wh./wl)); % from equating FOCs wrt skill supply, solve for taul

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

end