function f=target_equ(x, MOM, paramss, list, poll, targets)
% Model
% equilibrium for one period!
% takes policy as given

%- read in policy and parameters
read_in_pars_calib;

% choice variables

%- transform variables directly instead of in code
 hhf    = exp(x(list.choiceCALIB=='hhf'));
 hhg    = exp(x(list.choiceCALIB=='hhg'));
 hhn    = exp(x(list.choiceCALIB=='hhn'));
 hln    = exp(x(list.choiceCALIB=='hln'));
 hlf    = exp(x(list.choiceCALIB=='hlf'));
 hlg    = exp(x(list.choiceCALIB=='hlg'));
 C      = exp(x(list.choiceCALIB=='C'));
 F      = exp(x(list.choiceCALIB=='F'));
 G      = exp(x(list.choiceCALIB=='G'));
 Af     = exp(x(list.choiceCALIB=='Af')); % technology 2015-2019
 Ag     = exp(x(list.choiceCALIB=='Ag'));
 An     = exp(x(list.choiceCALIB=='An'));
 hl     = upbarH/(1+exp(x(list.choiceCALIB=='hl')));
 hh     = upbarH/(1+exp(x(list.choiceCALIB=='hh')));
 sf     = exp(x(list.choiceCALIB=='sf'));
 sg     = exp(x(list.choiceCALIB=='sg'));
 sn     = exp(x(list.choiceCALIB=='sn'));
 gammalh = x(list.choiceCALIB=='gammalh')^2;
 gammall = x(list.choiceCALIB=='gammall')^2;
 wh     = exp(x(list.choiceCALIB=='wh'));
 wl     = exp(x(list.choiceCALIB=='wl'));
 ws     = exp(x(list.choiceCALIB=='ws'));
 pg     = exp(x(list.choiceCALIB=='pg'));
 pn     = exp(x(list.choiceCALIB=='pn'));
 pe     = exp(x(list.choiceCALIB=='pe'));
 pf     = exp(x(list.choiceCALIB=='pf'));

% parameters
thetan = 1/(1+exp(x(list.choiceCALIB=='thetan')));
thetaf = 1/(1+exp(x(list.choiceCALIB=='thetaf')));
thetag = 1/(1+exp(x(list.choiceCALIB=='thetag')));
Af_lag = exp(x(list.choiceCALIB=='Af_lag'));
Ag_lag = exp(x(list.choiceCALIB=='Ag_lag'));
An_lag = exp(x(list.choiceCALIB=='An_lag'));
el     = exp(x(list.choiceCALIB=='el'));
eh     = exp(x(list.choiceCALIB=='eh'));
omegaa = exp(x(list.choiceCALIB=='omegaa'));
lambdaa = exp(x(list.choiceCALIB=='lambdaa'));
deltay = 1/(1+exp(x(list.choiceCALIB=='deltay')));

%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

SGov    = zh*(wh.*hh*eh-lambdaa.*(wh.*hh*eh).^(1-taul))...
            +zl*(wl.*hl*el-lambdaa.*(wl.*hl*el).^(1-taul))...
            +tauf.*omegaa*pf.*F;
            % subsidies and profits and wages scientists cancel
N       =  (1-deltay)/deltay.*(pe./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan).*An); % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag).*Ag);
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

%% equations
q=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- calibration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) perceived as if Af, Ag, An in 2015-2019 are parameters but from here back
% out Af0, Ag0, An0
q=q+1;
f(q) = MOM.FG-F/G; %Af => Af0 from LOM

q=q+1;
f(q) = E*pe/Y - MOM.EpeY; % market share Epe = determines deltay

q=q+1;
f(q) = Y - MOM.Y; % scales model!

q=q+1;
f(q) = wh/wl-MOM.whwl; %=> determines Af as fcn of Ag

%2) emissions
q=q+1;
f(q) = omegaa - MOM.emissionsUS2019/F;

%3) government
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))+tauf.*pf.*omegaa*F;

%4) skill shares
q=q+1;
f(q) = MOM.hhg_hhghlg-(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag)));

q=q+1;
f(q) = thetaf-thetan;%Y-xg-xn-xf-C; % => pg

q=q+1;
f(q) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;

%5) effective skill productivity
q=q+1;
f(q) = MOM.hhehzh_total-1/(1+zh/(1-zh)*hh/hl/el*eh); % => determines eleh

q=q+1;
f(q) = el-1; % => determines el

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1- household optimality (muu auxiliary variable determined above)
q=q+1;
f(q) = hh^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wh*eh)^(1-taul))+gammalh; %=> determines hh

q=q+1;
f(q) = hl^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wl*el)^(1-taul))+gammall; %=> determines hl

%3- budget
q=q+1;
f(q) = zh*lambdaa*(wh*hh*eh)^(1-taul)+zl*lambdaa*(wl*hl*el)^(1-taul)+SGov-C; %=> determines C

%4- output fossil
q=q+1;
f(q) = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*(Af.*Lf) -(F); 

%5- output neutral
q=q+1;
f(q) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f(q)=  G-(Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));

%6- demand green scientists
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag_lag);
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An_lag);

%8- LOM technology
q=q+1;
f(q) = An-An_lag.*(1+gammaa*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
q=q+1;
f(q) = Af-Af_lag*(1+gammaa*(sf/rhof)^etaa*(A_lag/Af_lag)^phii);
q=q+1;
f(q) = Ag-Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
 
%9- optimality labour input producers
q=q+1;
f(q) = thetan*Ln.*wln-wh.*hhn;
q=q+1;
f(q)= thetag*Lg.*wlg-wh.*hhg;
q=q+1;
f(q)=(1-thetan)*Ln.*wln-wl.*hln;
q=q+1;
f(q)=(1-thetag)*Lg.*wlg-wl.*hlg;

% prices and wages
%- optimality energy producers
q=q+1;
f(q) = pf *F.^(1/eppse)- (G).^(1/eppse).*pg; 

%- demand skill

q=q+1;
f(q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*wlf; % from optimality labour input producers fossil, and demand labour fossil
q=q+1;
f(q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*wlf;

%- definitions prices
q=q+1;
f(q) = pe - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
q=q+1;
f(q) = pn - ((1-deltay.*pe.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire


%- market clearing (consumption good=> numeraire)
q=q+1;
f(q) = zh*hh*eh-(hhn + hhf+hhg); % high skill market clearing
q=q+1;
f(q) = zl*hl*el-(hln + hlf+hlg); % low skill market clearing
q=q+1;
f(q) = S-(sn+sf+sg);

%13- Kuhn Tucker Labour supply
q=q+1;
f(q)= gammalh*(hh-upbarH);
q=q+1;
f(q)= gammall*(hl-upbarH);

%fprintf('number equations: %d; number variables %d', q, length(list.choiceCALIB));
end
