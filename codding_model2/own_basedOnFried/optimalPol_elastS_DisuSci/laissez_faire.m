function f=laissez_faire(x, params, list, pol, laggs)
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
 hhf    = exp(x(list.choice=='hhf'));
 hhg    = exp(x(list.choice=='hhg'));
 hhn    = exp(x(list.choice=='hhn'));
 hln    = exp(x(list.choice=='hln'));
 hlf    = exp(x(list.choice=='hlf'));
 hlg    = exp(x(list.choice=='hlg'));
 C      = exp(x(list.choice=='C'));
 F      = exp(x(list.choice=='F'));
 G      = exp(x(list.choice=='G'));
 Af     = exp(x(list.choice=='Af'));
 Ag     = exp(x(list.choice=='Ag'));
 An     = exp(x(list.choice=='An'));
 hl     = upbarH/(1+exp(x(list.choice=='hl')));
 hh     = upbarH/(1+exp(x(list.choice=='hh')));
 sff     = exp(x(list.choice=='sff'));
 sg     = exp(x(list.choice=='sg'));
 sn     = exp(x(list.choice=='sn'));
 gammalh = x(list.choice=='gammalh')^2;
 gammall = x(list.choice=='gammall')^2;
 wh     = exp(x(list.choice=='wh'));
 wl     = exp(x(list.choice=='wl'));
 ws     = exp(x(list.choice=='ws'));
 pg     = exp(x(list.choice=='pg'));
 pn     = exp(x(list.choice=='pn'));
 pee     = exp(x(list.choice=='pee'));
 pf     = exp(x(list.choice=='pf'));

%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl*hl).^(1-taul))...
            +tauf.*pf.*F;
            % subsidies and profits and wages scientists cancel
N       =  (1-deltay)/deltay.*(pee./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
% wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
q=q+1;
f(q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
%2
q=q+1;
f(q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl

%3- budget
q=q+1;
f(q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov-C; %=> determines C

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
f(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
%9
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

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

% prices and wages
%17- optimality energy producers
q=q+1;
f(q) = pf.*F.^(1/eppse)- (G).^(1/eppse).*pg; 

%18- demand skill
q=q+1;
f(q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
%19
q=q+1;
f(q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;

%- definitions prices
%20
q=q+1;
f(q) = pee - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
f(q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%22- market clearing (consumption good=> numeraire)
q=q+1;
f(q) = zh.*hh-(hhn + hhf+hhg); % high skill market clearing
%23
q=q+1;
f(q) = (1-zh).*hl-(hln + hlf+hlg); % low skill market clearing
%24
q=q+1;
f(q) = S-(sn+sff+sg);

%13- Kuhn Tucker Labour supply
%25
q=q+1;
f(q)= gammalh.*(hh-upbarH);
%26
q=q+1;
f(q)= gammall.*(hl-upbarH);

%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end
