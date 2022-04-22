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

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% equations
q=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- calibration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) perceived as if Af, Ag, An in 2015-2019 are parameters but from here back
% out Af0, Ag0, An0
%1
q=q+1;
f(q) = MOM.FG-F/G; %Af => Af0 from LOM
%f(q) = F/E-0.8;

%2
q=q+1;
f(q) = E*pe/Y - MOM.EpeY; % market share Epe = determines deltay

%3
q=q+1;
%q=q+1;
f(q) = Y - MOM.Y; % scales model!

%4
q=q+1;
f(q) = wh/wl-MOM.whwl; %=> determines ehel

%2) emissions
%5
q=q+1;
f(q) = omegaa - MOM.emissionsUS2019/F;

%3) government
%6
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))...
             +tauf.*pf.*omegaa*F;

%4) skill shares
%7
q=q+1;
f(q) = MOM.hhg_hhghlg-hhg/(hhg+hlg);%(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag)));

%8
q=q+1;
f(q) = Y-xg-xn-xf-C; % => pg 

%9
q=q+1;
f(q) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;

%5) effective skill productivity
%10
q=q+1;
f(q) =  thetaf-thetan;%
%f(q) = MOM.hhehzh_total-1/(1+zh/(1-zh)*hh/hl*eh/el); % => determines eleh;
%function in conflict with other shares

%11 % income share labour versus machines 
q=q+1;
f(q) = Ag/An-MOM.AgAn;%MOM.el; % => determines el
%f(q) = wlg/wln-MOM.wlgwln;%MOM.el; % => determines el

% f(q) = (zh*wh*eh*hh+zl*wl*el*hl+ws*(sf+sg+sn))/Y-MOM.labourshare;% => 0.66; el
%f(q) = C/Y-MOM.labourshare;% => 0.66; el

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1- household optimality (muu auxiliary variable determined above)
%12
q=q+1;
f(q) = hh^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wh*eh)^(1-taul))+gammalh; %=> determines hh

%13
q=q+1;
f(q) = hl^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wl*el)^(1-taul))+gammall; %=> determines hl

%3- budget
%14
q=q+1;
f(q) = zh*lambdaa*(wh*hh*eh)^(1-taul)+zl*lambdaa*(wl*hl*el)^(1-taul)+SGov-C; %=> determines C

%4- output fossil
%15
q=q+1;
f(q) = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*(Af.*Lf) -(F); 

%5- output neutral
%16
q=q+1;
f(q) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
%17
q=q+1;
f(q)=  G-(Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));

%6- demand green scientists
%18
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 
%19
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag_lag);
%20
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An_lag);

%8- LOM technology
%21
q=q+1;
f(q) = An-An_lag.*(1+gammaa*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
%22
q=q+1;
f(q) = Af-Af_lag*(1+gammaa*(sf/rhof)^etaa*(A_lag/Af_lag)^phii);
%23
q=q+1;
f(q) = Ag-Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
 
%9- optimality labour input producers

%24
q=q+1;
f(q) = thetan*Ln.*wln-wh.*hhn;
%25
q=q+1;
f(q)= thetag*Lg.*wlg-wh.*hhg;
%26
q=q+1;
f(q)=(1-thetan)*Ln.*wln-wl.*hln;
%27
q=q+1;
f(q)=(1-thetag)*Lg.*wlg-wl.*hlg;

% prices and wages
%- optimality energy producers

%28
q=q+1;
f(q) = pf *F.^(1/eppse)- (G).^(1/eppse).*pg; 

%- demand skill
%29
q=q+1;
f(q) = wh - thetaf*Lf/hhf.*wlf; % from optimality labour input producers fossil, and demand labour fossil
%30
q=q+1;
f(q) = wl-(1-thetaf)*Lf/hlf*wlf;

%- definitions prices
%31
q=q+1;
f(q) = pe - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition

%32
q=q+1;
f(q) = pn - ((1-deltay.*pe.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire


%- market clearing (consumption good=> numeraire)
%33
q=q+1;
f(q) = zh*hh*eh-(hhn+hhf+hhg); % high skill market clearing
%34
q=q+1;
f(q) = zl*hl*el-(hln+hlf+hlg); % low skill market clearing
%35
q=q+1;
f(q) = S-(sn+sf+sg);

%13- Kuhn Tucker Labour supply
%36
q=q+1;
f(q)= gammalh*(upbarH-hh);
%37
q=q+1;
f(q)= gammall*(upbarH-hl);

%fprintf('number equations: %d; number variables %d', q, length(list.choiceCALIB));
end