function f=target_equ(x, params, list, pol, laggs, targets)
% Model
% equilibrium for one period!
% takes policy as given

%- read in policy and parameters
read_in_params;
read_in_pol;


% choice variables

%- transform variables directly instead of in code
 hhf    = exp(x(list.calibJoint=='hhf'));
 hhg    = exp(x(list.calibJoint=='hhg'));
 hhn    = exp(x(list.calibJoint=='hhn'));
 hln    = exp(x(list.calibJoint=='hln'));
 hlf    = exp(x(list.calibJoint=='hlf'));
 hlg    = exp(x(list.calibJoint=='hlg'));
 C      = exp(x(list.calibJoint=='C'));
 F      = exp(x(list.calibJoint=='F'));
 G      = exp(x(list.calibJoint=='G'));
 Af     = exp(x(list.calibJoint=='Af')); % technology 2015-2019
 Ag     = exp(x(list.calibJoint=='Ag'));
 An     = exp(x(list.calibJoint=='An'));
 hl     = upbarH/(1+exp(x(list.calibJoint=='hl')));
 hh     = upbarH/(1+exp(x(list.calibJoint=='hh')));
 sf     = exp(x(list.calibJoint=='sf'));
 sg     = exp(x(list.calibJoint=='sg'));
 sn     = exp(x(list.calibJoint=='sn'));
 gammalh = x(list.calibJoint=='gammalh')^2;
 gammall = x(list.calibJoint=='gammall')^2;
 wh     = exp(x(list.calibJoint=='wh'));
 wl     = exp(x(list.calibJoint=='wl'));
 ws     = exp(x(list.calibJoint=='ws'));
 pg     = exp(x(list.calibJoint=='pg'));
 pn     = exp(x(list.calibJoint=='pn'));
 pe     = exp(x(list.calibJoint=='pe'));
 pf     = exp(x(list.calibJoint=='pf'));

% parameters
thetan = 1/(1+exp(x(list.calibJoint=='thetan')));
thetaf = 1/(1+exp(x(list.calibJoint=='thetaf')));
HAVE TO THINK ABOUT LAGEGD VARIABLES HERE
Af_lag = exp(x(list.calibJoint=='Af_lag'));
Ag_lag = exp(x(list.calibJoint=='Ag_lag'));
An_lag = exp(x(list.calibJoint=='An_lag'));
el     = exp(x(list.calibJoint=='el'));
eh     = exp(x(list.calibJoint=='eh'));
omegaa = exp(x(list.calibJoint=='omegaa'));
lambdaa = exp(x(list.calibJoint=='lambdaa'));
deltay = 1/(1+exp(x(list.calibJoint=='deltay')));

%- read in auxiliary equations
[Af_lag, Ag_lag, An_lag, A_lag, Lg, Ln, Lf, muu, E, SGov, N, Y,wln, wlg, wlf, xn, xg, xf ] ...
= auxiliary_stuff(params, list, pol, targets, laggs, C, hhg, hhf, hhn, hlg, hln, hlf, F, G, wh, hh, hl, wl, An,...
                  Ag, Af, pn, pe, pf, pg);

%% equations
q=0;

%- calibration 
q=q+1;
f(q) = MOM.FG-F/G; %Af

q=q+1;
f(q) = E*pe/Y - MOM.EpeY; % market share Epe = determines deltay

q=q+1;
f(q) = Y - MOM.Y; % scales model!

q=q+1;
f(q) = omegaa - MOM.emissionsUS2019/F;

q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))+tauf.*pf.*omegaa*F;
q=q+1;
f(q) = whg/wl-MOM.whwl; %=> determines Af as fcn of Ag

q=q+1;
f(q) = Y-xg-xn-xf-C; % => thetan

q=q+1;
f(q) = MOM.hhehzh_total-1/(1+zh/(1-zh)*hh/hl/el*eh); % => determines eleh

q=q+1;
f(q) = el-1; % => determines el

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

fprintf('number equations: %d; number variables %d', q, length(list.calibJoint));
end
