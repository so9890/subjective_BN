function [c, ceq]=laissez_faire_NoSkill_nows_fmincon(x, params, list, pol, laggs)
% Model without skill

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
%  hhf    = exp(x(list.choice=='hhf'));
%  hhg    = exp(x(list.choice=='hhg'));
%  hhn    = exp(x(list.choice=='hhn'));
%  hln    = exp(x(list.choice=='hln'));
%  hlf    = exp(x(list.choice=='hlf'));
%  hlg    = exp(x(list.choice=='hlg'));
 C      = exp(x(list.choice=='C'));
 F      = exp(x(list.choice=='F'));
 G      = exp(x(list.choice=='G'));
 Af     = exp(x(list.choice=='Af'));
 Ag     = exp(x(list.choice=='Ag'));
 An     = exp(x(list.choice=='An'));
 hl     = upbarH/(1+exp(x(list.choice=='hl')));
 hh     = upbarH/(1+exp(x(list.choice=='hh')));
 
 Lf    = exp(x(list.choice=='Lf'));
 Lg     = exp(x(list.choice=='Lg'));
 
 sff    = exp(x(list.choice=='sff'));
 sg     = exp(x(list.choice=='sg'));
 sn     = exp(x(list.choice=='sn'));
 S      = exp(x(list.choice=='S')); % total labour supply
%  wse     = exp(x(list.choice=='wse'));
 ws     = exp(x(list.choice=='ws'));
%  gammasf = x(list.choice=='gammasf')^2;
%  gammasn = x(list.choice=='gammasn')^2;
%  gammasg = x(list.choice=='gammasg')^2;
 gammalh = x(list.choice=='gammalh')^2;
 gammall = x(list.choice=='gammall')^2;
 wh     = exp(x(list.choice=='wh'));
 wl     = exp(x(list.choice=='wl'));
 pg     = exp(x(list.choice=='pg'));
 pn     = exp(x(list.choice=='pn'));
 pee     = exp(x(list.choice=='pee'));
 pf     = exp(x(list.choice=='pf'));
 lambdaa  = exp(x(list.choice=='lambdaa'));

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
c=[];
q=0;

%1- household optimality (muu auxiliary variable determined above)
q=q+1;
ceq(q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
%2
q=q+1;
ceq(q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl

%3- budget
q=q+1;
ceq(q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov-C; %=> determines C

%4- output fossil
q=q+1;
ceq(q) = ((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%5- output neutral
q=q+1;
ceq(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
ceq(q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%7- demand green scientists
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 

%8
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);

%9
q=q+1;
ceq(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

%10- LOM technology
q=q+1;
ceq(q) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
%11
q=q+1;
ceq(q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
%12
q=q+1;
ceq(q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);
 
%13- optimality labour input producers
q=q+1;
ceq(q) = thetan*Ln.*wln-wh.*hhn;
%14
q=q+1;
ceq(q)= thetag*Lg.*wlg-wh.*hhg;
%15
q=q+1;
ceq(q)=(1-thetan)*Ln.*wln-wl.*hln;
%16
q=q+1;
ceq(q)=(1-thetag)*Lg.*wlg-wl.*hlg;

% prices and wages
%17- optimality energy producers
q=q+1;
ceq(q) = pf.*F.^(1/eppse)- (G).^(1/eppse).*pg; 

%18- demand skill
q=q+1;
ceq(q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
%19
q=q+1;
ceq(q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;

%- definitions prices
%20
q=q+1;
ceq(q) = pee - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
ceq(q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%22- market clearing (consumption good=> numeraire)
q=q+1;
ceq(q) = zh.*hh-(hhn + hhf+hhg); % high skill market clearing
%23
q=q+1;
ceq(q) = (1-zh).*hl-(hln + hlf+hlg); % low skill market clearing

%13- Kuhn Tucker Labour supply and scientists
%25
q=q+1;
ceq(q)= gammalh.*(upbarH-hh);
%26
q=q+1;
ceq(q)= gammall.*(upbarH-hl);

q=q+1;
ceq(q)= S-(ws*muu/chiis)^(1/sigmaa); % scientist hours supply
q=q+1;
ceq(q)= sff+sg+sn-S; % determines wage in neutral sector

%26
% q=q+1;
% ceq(q)= gammasf.*(S-sff);
% 
% q=q+1;
% ceq(q)= gammasn.*(S-sn);
% 
% q=q+1;
% ceq(q)= gammasg.*(S-sg);

% balanced budget
q=q+1;
ceq(q)= SGov;
%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end