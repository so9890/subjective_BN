function f=laissez_faire_xgrowth(x, params, list, pol, laggs, indic, MOM,t, Emlim)
% Model
% equilibrium for one period!
% takes policy as given

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
read_in_pol;

%- initial condition
Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');
% choice variables
%- transform variables directly instead of in code
if indic.noskill==0
     hhf    = exp(x(list.choice=='hhf'));
     hhg    = exp(x(list.choice=='hhg'));
     hhn    = exp(x(list.choice=='hhn'));
     hln    = exp(x(list.choice=='hln'));
     hlf    = exp(x(list.choice=='hlf'));
     hlg    = exp(x(list.choice=='hlg'));
     hl     = upbarH/(1+exp(x(list.choice=='hl')));
     hh     = upbarH/(1+exp(x(list.choice=='hh')));
     gammall = x(list.choice=='gammall')^2;
     wh     = exp(x(list.choice=='wh'));
     wl     = exp(x(list.choice=='wl'));
else
    h      = upbarH/(1+exp(x(list.choice=='h')));
    w      = exp(x(list.choice=='w'));
    Ln     = exp(x(list.choice=='Ln'));
    Lg     = exp(x(list.choice=='Lg'));
    Lf     = exp(x(list.choice=='Lf'));
end

 C      = exp(x(list.choice=='C'));

 F      = exp(x(list.choice=='F'));
 G      = exp(x(list.choice=='G'));
 gammalh = x(list.choice=='gammalh')^2;

 pg     = exp(x(list.choice=='pg'));
 pn     = exp(x(list.choice=='pn'));
 pee     = exp(x(list.choice=='pee'));
 pf     = exp(x(list.choice=='pf'));
if indic.notaul ==6
 taul = x(list.choice=='lambdaa');
else
 lambdaa  = (x(list.choice=='lambdaa')); % in calibration chosen to match GovRev
end

if indic.limit_LF==1
 tauf=x(list.choice=='tauf');
end
%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);
%- technology
% An=(1+vn)*An_lag;
% Ag=(1+vg)*Ag_lag;
% Af=(1+vf)*Af_lag;
    sn=   sn0; %x0LF(list.choice=='') 0.0034; %0.01*MOM.targethour;
    sff=  sff0; % 5.5660e-10; %0.01*MOM.targethour;
    sg=   sg0; %1.0305e-08; %0.01*MOM.targethour;


An=An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
Af=Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
Ag=Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);
 

if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf); 
else
    hh=h; hl=h; wh=w; wl=w; % this should suffice to have governmenta budget correct
end

if indic.notaul<2 || ...
   indic.notaul == 6 % tauf redistributed via income tax
    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +tauf.*F;
    Tls =zeros(size(F));    
    GovCon =zeros(size(F));
elseif indic.notaul == 2 ||...
        indic.notaul==3 %2,3,4,5,7
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon = tauf.*F; % GovCon = env tax consumed by government
    Tls =zeros(size(F));    
elseif indic.notaul == 4 || indic.notaul ==5
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon =zeros(size(F));
    Tls  = tauf.*F;
elseif indic.notaul == 7 % earmarking
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
    GovCon =zeros(size(F));
    Tls =zeros(size(F));    
    taus = tauf.*F./(pg.*G); % subsidy on green sector
end
muu   = C.^(-thetaa); % same equation in case thetaa == 1
   
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

N       =  (1-deltay)/deltay.*(pee./pn)^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = (pg.*(1+taus)).^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
if indic.noskill==0
    
%    if indic.ineq==0
        q=q+1;
        f(q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f(q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl
        %3- budget
        q=q+1;
        f(q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+Tls-C; %=> determines C

%     else  
%         error('version with inequality not yet updated')
%         q=q+1;
%         f(q)= chii*hh.^(sigmaa+taul)- ((muuh.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh.*hh.^taul); %=> determines hh
%         %2
%         q=q+1;
%         f(q)= chii*hl.^(sigmaa+taul) - ((muul.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall.*hl.^taul); %=> determines hl
%     
%         %3- budget
%         q=q+1;
%         f(q) = lambdaa.*(wl.*hl).^(1-taul)+Tls-Cl; %=> determines C
%         q=q+1;
%         f(q) = lambdaa.*(wh.*hh).^(1-taul)+Tls-Ch;
%     end
else
    q=q+1;
    f(q)= chii*h.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(w).^(1-taul))-gammalh.*h.^taul); %=> determines hh
   %3- budget
    q=q+1;
    f(q) = lambdaa.*(w.*h).^(1-taul)+Tls-C; %=> determines C
end


%4- output fossil
q=q+1;
f(q) = (alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%5- output neutral
q=q+1;
f(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f(q)=  G-Ag.*Lg.*(pg.*(1+taus).*alphag).^(alphag./(1-alphag));

%13- optimality labour input producers
if indic.noskill==0

    q=q+1;
    f(q) = thetan*Ln.*wln-wh.*hhn;
    %14
    q=q+1;
    f(q)= thetag*Lg.*wlg-wh.*hhg;
    %15
    q=q+1;
    f(q)=(1-thetan)*Ln.*wln-wl.*hln;
    %16
    q=q+1;
    f(q)=(1-thetag)*Lg.*wlg-wl.*hlg;
    %18- demand skill
    q=q+1;
    f(q) = wh - thetaf*Lf./hhf.*wlf; % from optimality labour input producers fossil, and demand labour fossil
    %19
    q=q+1;
    f(q) = wl-(1-thetaf)*Lf./hlf.*wlf;

else
    q=q+1; % labour demand fossil
    f(q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    q=q+1;
    f(q) = Ln -pn.*(1-alphan).*N./w; % labour demand neutral 
    q=q+1;
    f(q) = Lg - pg.*(1+taus)*(1-alphag).*G./w;
end

% prices and wages
%17- optimality energy producers
q=q+1;
f(q) = (pf+tauf).*F.^(1/eppse)- (G).^(1/eppse).*pg; 


%- definitions prices
%20
q=q+1;
f(q) = pee - ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
f(q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%22- market clearing (consumption good=> numeraire)
if indic.noskill==0

    q=q+1;
    f(q) = zh.*hh-(hhn + hhf+hhg); % high skill market clearing
    %23
    q=q+1;
    f(q) = (1-zh).*hl-(hln + hlf+hlg); % low skill market clearing

    %13- Kuhn Tucker Labour supply and scientists
    %25
    q=q+1;
    f(q)= gammalh.*(upbarH-hh);
    %26
    q=q+1;
    f(q)= gammall.*(upbarH-hl);
else
    q=q+1;
    f(q) = Lf+Lg+Ln-h;
    q=q+1;
    f(q)= gammalh.*(h-upbarH);
end


% balanced budget
q=q+1;
f(q)= SGov-GovRev*Y;
% if indic.noskill==0
%     f(q)= SGov-GovRev*(wh.*hh+wl.*hl);
% else
%     f(q)= SGov-GovRev*(w.*h);
% end

%- if emission limit determines tauf
if indic.limit_LF==1
q=q+1;
if t==1 % base year period, tauf =0
    f(q)= tauf ;
else
    f(q)=omegaa*F-deltaa-Emlim;
end
end

% fprintf('number equations: %d; number variables %d', q, length(list.choice));
end