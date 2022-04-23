function f=laissez_faire(x, params, list, pol, laggs, targets)
% Model
% equilibrium for one period!
% takes policy as given
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
read_in_pol;

%- initial condition
Af_lag=laggs(list.laggs=='Af_lag');
An_lag=laggs(list.laggs=='An_lag');
Ag_lag=laggs(list.laggs=='Ag_lag');

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
 sf     = exp(x(list.choice=='sf'));
 sg     = exp(x(list.choice=='sg'));
 sn     = exp(x(list.choice=='sn'));
 gammalh = x(list.choice=='gammalh')^2;
 gammall = x(list.choice=='gammall')^2;
 wh     = exp(x(list.choice=='wh'));
 wl     = exp(x(list.choice=='wl'));
 ws     = exp(x(list.choice=='ws'));
 pg     = exp(x(list.choice=='pg'));
 pn     = exp(x(list.choice=='pn'));
 pe     = exp(x(list.choice=='pe'));
 pf     = exp(x(list.choice=='pf'));

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

%% model equations
q=0;

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
f(q) = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*Af.*Lf -F; 

%5- output neutral
q=q+1;
f(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f(q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%7- demand green scientists
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 
%8
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(rhog^etaa.*(1-taus)*Ag./Ag_lag);
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(rhon^etaa.*An./An_lag);

% q=q+1;
% f(q)= (1-taus)*(Ag./Ag_lag.*ws)*rhog^etaa*sg-((gammaa*etaa*(A_lag./Ag_lag).^phii.*sg^(etaa).*pg.*G.*(1-alphag)));
% 
% %7- wage scientists neutral
% q=q+1;
% f(q)= ws.*sn.*An*rhon^etaa- (etaa*gammaa*An_lag.^(1-phii).*A_lag.^phii.*sn.^etaa.*pn.*(1-alphan)*N); 
% 
% % scientists fossil
% q=q+1;
% f(q) =  ws*(rhof^etaa.*Af./Af_lag)*sf- (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa).*pf.*(1-tauf).*F*(1-alphaf)); 

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
f(q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
q=q+1;
f(q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;

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

%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end
