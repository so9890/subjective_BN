function f=laissez_faire_nows_sep_bgp(x, params, list, pol, laggs, indic)

% to simulate economy assuming constant growth rates
% the solution will be growth rates given initial values; hence only one
% period! 
% and assumung that emissions do not grow further! 
% code assumes an interior solution

% Model
% equilibrium for one period!
% takes policy as given

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
read_in_pol;

%- initial condition: Technology level in last direct optimization period
Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');

% choice variables
%- transform variables directly instead of in code
if indic.noskill==0
     ghhf    = exp(x(list.choice=='hhf'));
     ghhg    = exp(x(list.choice=='hhg'));
     ghhn    = exp(x(list.choice=='hhn'));
     ghln    = exp(x(list.choice=='hln'));
     ghlf    = exp(x(list.choice=='hlf'));
     ghlg    = exp(x(list.choice=='hlg'));
     ghl     = upbarH/(1+exp(x(list.choice=='hl')));
     ghh     = upbarH/(1+exp(x(list.choice=='hh')));
     gwh     = exp(x(list.choice=='wh'));
     gwl     = exp(x(list.choice=='wl'));
else
    gh      = upbarH/(1+exp(x(list.choice=='h')));
    gw      = exp(x(list.choice=='w'));
    gLn     = exp(x(list.choice=='Ln'));
    gLg     = exp(x(list.choice=='Lg'));
    gLf     = exp(x(list.choice=='Lf'));
end


 gC      = exp(x(list.choice=='C'));

%  gF      = exp(x(list.choice=='F'));
 gG      = exp(x(list.choice=='G'));
 gAf     = exp(x(list.choice=='Af'));
 gAg     = exp(x(list.choice=='Ag'));
 gAn     = exp(x(list.choice=='An'));
 gsff     = upbarH/(1+exp(x(list.choice=='sff')));%exp(x(list.choice=='S')); % total labour supply
 gsg      = upbarH/(1+exp(x(list.choice=='sg')));%exp(x(list.choice=='S')); % total labour supply
 gsn      = upbarH/(1+exp(x(list.choice=='sn')));%exp(x(list.choice=='S')); % total labour supply

%  S      = upbarH/(1+exp(x(list.choice=='S')));%exp(x(list.choice=='S')); % total labour supply
 gwsf     = exp(x(list.choice=='wsf'));
 gwsn     = exp(x(list.choice=='wsn'));
 gwsg     = exp(x(list.choice=='wsg'));


 gpg     = exp(x(list.choice=='pg'));
 gpn     = exp(x(list.choice=='pn'));
 gpee     = exp(x(list.choice=='pee'));
 gpf     = exp(x(list.choice=='pf'));
 glambdaa  = exp(x(list.choice=='lambdaa'));

%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

if indic.noskill==0
    gLg      = ghhg.^thetag.*ghlg.^(1-thetag);
    gLn      = ghhn.^thetan.*ghln.^(1-thetan);
    gLf      = ghhf.^thetaf.*ghlf.^(1-thetaf); 
    
    if indic.notaul<2 % tauf redistributed via income tax
        SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
    else
        SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));
        if indic.notaul <4
            GovCon = tauf.*pf.*F;
        else
            GovCon =zeros(size(F));
        end
    end
else
    if indic.notaul<2
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
    else
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul));
        if indic.notaul <4
            GovCon = tauf.*pf.*F;
        else
            GovCon =zeros(size(F));
        end
    end
end
% lump sum transfers
if indic.notaul >=4
    Tls =tauf.*pf.*F;
else
    Tls =zeros(size(F));
end

muu   = C.^(-thetaa); % same equation in case thetaa == 1
   
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

N       =  (1-deltay)/deltay.*(pee./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
if indic.noskill==0
    
    if indic.ineq==0
        q=q+1;
        f(q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f(q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl
        %3- budget
        q=q+1;
        f(q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov+Tls-C; %=> determines C

    else  
        error('version with inequality not yet updated')
        q=q+1;
        f(q)= chii*hh.^(sigmaa+taul)- ((muuh.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f(q)= chii*hl.^(sigmaa+taul) - ((muul.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall.*hl.^taul); %=> determines hl
    
        %3- budget
        q=q+1;
        f(q) = lambdaa.*(wl.*hl).^(1-taul)+SGov-Cl; %=> determines C
        q=q+1;
        f(q) = lambdaa.*(wh.*hh).^(1-taul)+SGov-Ch;
    end
else
    q=q+1;
    f(q)= chii*h.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(w).^(1-taul))-gammalh.*h.^taul); %=> determines hh
   %3- budget
    q=q+1;
    f(q) = lambdaa.*(w.*h).^(1-taul)+SGov-C; %=> determines C
end

%4- output fossil
q=q+1;
f(q) = ((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%5- output neutral
q=q+1;
f(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f(q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%7- demand green scientists
q=q+1;
f(q)= wsf - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 

%8
q=q+1;
f(q)= wsg - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);

%9
q=q+1;
f(q)= wsn - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

%10- LOM technology
q=q+1;
f(q) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
%11
q=q+1;
f(q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
%12
q=q+1;
f(q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);
 
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
    f(q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    q=q+1;
    f(q) = Ln -pn.*(1-alphan).*N./w; % labour demand neutral 
    q=q+1;
    f(q) = Lg - pg.*(1-alphag).*G./w;
end

% prices and wages
%17- optimality energy producers
q=q+1;
f(q) = pf.*F.^(1/eppse)- (G).^(1/eppse).*pg; 


%- definitions prices
%20
q=q+1;
f(q) = pee - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
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

% optimality scientists
q=q+1;
f(q)= (chiis)*sff^sigmaas-(wsf-gammasf); % scientist hours supply
q=q+1;
f(q)= (chiis)*sg^sigmaas-((wsg-gammasg));
q=q+1;
f(q)= (chiis)*sn^sigmaas-((wsn-gammasn));
% Kuhn tucker scientists
q=q+1;
f(q)= gammasf.*(sff-upbarH);
q=q+1;
f(q)= gammasg.*(sg-upbarH);
q=q+1;
f(q)= gammasn.*(sn-upbarH);
%market for scientists
% q=q+1;
% f(q)= sff+sg+sn-S; % determines wage in neutral sector

% balanced budget
q=q+1;
f(q)= SGov;
%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end