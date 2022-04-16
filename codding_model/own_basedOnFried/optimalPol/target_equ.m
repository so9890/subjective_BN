function f=target_equ(x, MOM, paramss, list, poll)

% read in already calibrated parameters and policy
read_in_pars_calib;

% variables
hhf    = exp(x(list.choiceCALIB=='hhf'));
hhg    = exp(x(list.choiceCALIB=='hhg'));
hhn    = exp(x(list.choiceCALIB=='hhn'));
hln    = exp(x(list.choiceCALIB=='hln'));
hlf    = exp(x(list.choiceCALIB=='hlf'));
hlg    = exp(x(list.choiceCALIB=='hlg'));
C      = exp(x(list.choiceCALIB=='C'));
F      = exp(x(list.choiceCALIB=='F'));
G      = exp(x(list.choiceCALIB=='G'));
Af     = exp(x(list.choiceCALIB=='Af'));
Ag     = exp(x(list.choiceCALIB=='Ag'));
An     = exp(x(list.choiceCALIB=='An'));
hl     = upbarH-exp(x(list.choiceCALIB=='hl'));
hh     = upbarH-exp(x(list.choiceCALIB=='hh'));
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
Y     = exp(x(list.choiceCALIB=='Y'));
% parameters to be matched
thetaf = 1-exp(x(list.choiceCALIB=='thetaf'));
thetan = 1-exp(x(list.choiceCALIB=='thetan'));
thetag = 1-exp(x(list.choiceCALIB=='thetag'));

lambdaa = exp(x(list.choiceCALIB=='lambdaa'));
Af0     = exp(x(list.choiceCALIB=='Af0'));
Ag0     = exp(x(list.choiceCALIB=='Ag0'));
An0     = exp(x(list.choiceCALIB=='An0'));
A_lag   = max([Af0, An0, Ag0]);
% auxiliary equations
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
muu      = C.^(-thetaa);
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +zl*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*omegaa*F;
N       =  ((1-deltay)/deltay.*pe./pn)^(eppsy).*E; % demand N final good producers 
%
wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan).*An); % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag).*Ag);


%% matching equations: Baseyear: 2015-2019
% note: the data refers to 2015-2019
% => Af0, Ag0, An0 refer to 2010-2015 
% => take Af, Ag, An as initial values! 
%    to start simulation in 2020

q=0;
%1 final output
q=q+1;
f(q) = Y-MOM.Y;

%2 energy expenditure share
q=q+1;
f(q) = E*pe/Y-MOM.EpeY;

%3 fossil to renewable energy consumption
q=q+1;
f(q) = F/G -MOM.FG;

%4 Government budget;
q=q+1;
f(q) = MOM.Debt-SGov;

%5 wage premium skill
q=q+1;
f(q) = wh/wl-MOM.whwl;

%6 hours worked
q=q+1;
f(q) = 

%7 relative skill input neutral sector
q=q+1;
f(q) = hhn/hln - MOM.hhnhln;

%% MODEL equations laissez faire
%8- household optimality (muu auxiliary variable determined above)
q=q+1;
f(q) = hh^(sigmaa+taul)-(muu*lambdaa*(1-taul)*wh^(1-taul))+gammalh; %=> determines hh

%9
q=q+1;
f(q) = hl^(sigmaa+taul)-(muu*lambdaa*(1-taul)*wl^(1-taul))+gammall; %=> determines hl

%10- budget
q=q+1;
f(q) = zh*lambdaa*(wh*hh)^(1-taul)+zl*lambdaa*(wl*hl)^(1-taul)+SGov-C; %=> determines C

%11- output fossil
q=q+1;
f(q) = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*(Af.*Lf) -(F); 

%12- output neutral
q=q+1;
f(q) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f(q)=  G-(Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));

%6- demand green scientists
q=q+1;
f(q)= (1-taus)*(Ag./Ag0.*ws)*rhog^etaa*sg-((gammaa*etaa*(A_lag./Ag0).^phii.*sg^(etaa).*pg.*G.*(1-alphag)));

%7- wage scientists neutral
q=q+1;
f(q)= ws.*sn.*An*rhon^etaa- (etaa*gammaa*An0.^(1-phii).*A_lag.^phii.*sn.^etaa.*pn.*(1-alphan)*N); 

% scientists fossil
q=q+1;
f(q) =  ws*(rhof^etaa.*Af./Af0)*sf- (gammaa*etaa*(A_lag./Af0).^phii.*sf.^(etaa).*pf.*(1-tauf).*F*(1-alphaf)); 

%8- LOM technology
q=q+1;
f(q) = An-An0.*(1+gammaa*(sn./rhon).^etaa.*(A_lag./An0).^phii);
q=q+1;
f(q) = Af-Af0*(1+gammaa*(sf/rhof)^etaa*(A_lag/Af0)^phii);
q=q+1;
f(q) = Ag-Ag0*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag0)^phii);
 
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
f(q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
q=q+1;
f(q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;

%- definitions prices
q=q+1;
f(q) = pe - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
q=q+1;
f(q) = pn - ((1-deltay^eppsy.*pe.^(1-eppsy))./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

%- market clearing (consumption good=> numeraire)
q=q+1;
f(q) = zh*hh-(hhn+ hhf+hhg); % high skill market clearing
q=q+1;
f(q) = zl*hl-(hln + hlf+hlg); % low skill market clearing
q=q+1;
f(q) = S-(sn+sf+sg);

% production function Y 
q=q+1;
f(q) = Y - (deltay.*E.^((eppsy-1)/eppsy)+(1-deltay).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); 

%13- Kuhn Tucker Labour supply
q=q+1;
f(q)= gammalh*(hh-upbarH);
q=q+1;
f(q)= gammall*(hl-upbarH);
end
