function [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S, gammac]= OPT_aux_vars_notaus_flex(x, list, params, T, init201519, indic)

read_in_params;

% if indic.noskill==0
    hhf    = x((find(list.opt=='hhf')-1)*T+1:find(list.opt=='hhf')*T);
    hhg    = x((find(list.opt=='hhg')-1)*T+1:(find(list.opt=='hhg'))*T);
    hlf    = x((find(list.opt=='hlf')-1)*T+1:find(list.opt=='hlf')*T);
    hlg    = x((find(list.opt=='hlg')-1)*T+1:find(list.opt=='hlg')*T);
    hhn    = x((find(list.opt=='hhn')-1)*T+1:find(list.opt=='hhn')*T);
    hln    = x((find(list.opt=='hln')-1)*T+1:find(list.opt=='hln')*T);
    hl     = x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T);
    hh     = x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T);  
% else
%     
%     Lf    = x((find(list.opt=='Lf')-1)*T+1:find(list.opt=='Lf')*T);
%     Lg     = x((find(list.opt=='Lg')-1)*T+1:find(list.opt=='Lg')*T);
%     h     = x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T); 
% end

if indic.ineq==0
    C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
    Ch=C; Cl=C;
else
    Cl      = x((find(list.opt=='Cl')-1)*T+1:find(list.opt=='Cl')*T);
    Ch      = x((find(list.opt=='Ch')-1)*T+1:find(list.opt=='Ch')*T);
    C=Ch;
end
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);

if indic.xgrowth==0
    sff     = x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T);
    sg     = x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T);
    sn     = x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T);
else
    sff=zeros(size(F));
    sn=zeros(size(F));
    sg=zeros(size(F));
end
 

%% auxiliary variables

% hhn     = zh*hh-(hhf+hhg);
% hln     = (1-zh)*hl-(hlf+hlg); 
hhhl    = hh./hl;
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
% loop over technology
if indic.xgrowth==0
    Af=zeros(T,1);
    Af_lag=[init201519(list.init=='Af0'); Af(1:T)]; % drop last value later
    Ag=zeros(T,1);
    Ag_lag=[init201519(list.init=='Ag0'); Ag(1:T)];
    An=zeros(T,1);
    An_lag=[init201519(list.init=='An0'); An(1:T)];
    A_lag=zeros(T,1);


    for i=1:T
        A_lag(i)   = (rhof*Af_lag(i)+rhon*An_lag(i)+rhog*Ag_lag(i))./(rhof+rhon+rhog);

        Af(i)=Af_lag(i).*(1+gammaa*(sff(i)/rhof).^etaa.*(A_lag(i)/Af_lag(i))^phii);
        Ag(i)=Ag_lag(i).*(1+gammaa*(sg(i)/rhog).^etaa.*(A_lag(i)/Ag_lag(i))^phii);
        An(i)=An_lag(i).*(1+gammaa*(sn(i)/rhon).^etaa.*(A_lag(i)/An_lag(i))^phii);

        %-update lags

        Af_lag(i+1)=Af(i);
        Ag_lag(i+1)=Ag(i);
        An_lag(i+1)=An(i);

    end

    Af_lag=Af_lag(1:end-1);
    An_lag=An_lag(1:end-1);
    Ag_lag=Ag_lag(1:end-1);
else
    An_lag=init201519(list.init=='An0');
    Ag_lag=init201519(list.init=='Ag0');
    Af_lag=init201519(list.init=='Af0');

    Af=zeros(size(F));
    Ag=zeros(size(F));
    An=zeros(size(F));

    for i=1:T
        An(i)=(1+vn)*An_lag;
        Ag(i)=(1+vg)*Ag_lag;
        Af(i)=(1+vf)*Af_lag;
        %- update laggs
        An_lag=An(i);
        Af_lag=Af(i);
        Ag_lag=Ag(i);
    end
    
    An_lag=An;
    Af_lag=Af;
    Ag_lag=Ag;
    A_lag=Ag_lag;
end

if indic.ineq==0
    if indic.BN==0
        muu      = C.^(-thetaa); % same equation in case thetaa == 1
    else
        muu =-(C-B).^(zetaa-1);
    end
    muuh=muu; muul=muu;
else
    if indic.BN==0
        muuh      = Ch.^(-thetaa); % same equation in case thetaa == 1
        muul      = Cl.^(-thetaa); % same equation in case thetaa == 1
    else
        muul =-(Cl-Bl).^(zetaa-1);
        muuh =-(Ch-Bh).^(zetaa-1);
    end
    muu=muuh; 
end
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

%- wages scientists  
if indic.xgrowth==0
    wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
    wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
    wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

    %- relevant for code without separate markets
    S    = sn+sg+sff;
    ws   = chiis*S.^sigmaas; 
else
    ws=zeros(size(F));
    S=zeros(size(F));
    wsf =zeros(size(F));
    wsg=zeros(size(F));
    wsn=zeros(size(F));
end
% assuming interior solution households
if indic.notaul==0
    if indic.ineq==0
        taul   = (log(wh./wl)-sigmaa*log(hhhl))./(log(hhhl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
    else
        taul   = (log(wh./wl)+log(muuh./muul)-sigmaa*log(hhhl))./(log(hhhl)+log(wh./wl)); % from equating FOCs wrt skill supply, solve for taul        
    end
else
    taul   = zeros(size(sn));
end
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

gammac =(1+gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii)-1;

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

 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
end

 SWF = Utilcon-Utillab-Utilsci;

 contUtil= Utilcon(T)/(1-betaa* (1+gammac)^(1-thetaa));
 contUtillab =1/(1-betaa)*(chii.*(zh.*hh(T).^(1+sigmaa)+(1-zh).*hl(T).^(1+sigmaa))./(1+sigmaa));
 contUtilsci = 1/(1-betaa)*(chiis*S(T).^(1+sigmaas)./(1+sigmaas));
 
 PVcontUtil = contUtil-contUtillab-contUtilsci;
end