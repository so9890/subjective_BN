function f=target_equFinal(x, MOM, paramss, list, poll, targets)
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
chii = exp(x(list.choiceCALIB=='chii'));
lambdaa = exp(x(list.choiceCALIB=='lambdaa'));
deltay = 1/(1+exp(x(list.choiceCALIB=='deltay')));

%% - read in auxiliary equations
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
            
N       =  (1-deltay)/deltay.*(pe./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn*(1-alphan)*N/Ln; %pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg*(1-alphag)*G/Lg;%pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = pf*(1-tauf)*(1-alphaf)*F/Lf; %(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

% machine prices and taxes
pxf     = 1;
pxg     = 1;
pxn     = 1; 
zetaf = (1-alphaf)/alphaf;
zetag = (1-alphag)/alphag; 
zetan = (1-alphan)/alphan;

% profits
prof_f =  (1-alphaf)/alphaf*xf-ws*sf;
prof_g =  (1-alphag)/alphag*xg-ws*(1-taus)*sg;
prof_n =  (1-alphan)/alphan*xn-ws*sn;

SGov    = zh*(wh.*hh*eh-lambdaa.*(wh.*hh*eh).^(1-taul))...
            +zl*(wl.*hl*el-lambdaa.*(wl.*hl*el).^(1-taul))...
            +tauf.*omegaa*pf.*F+ws*sf+ws*sg+ws*sn+ prof_f+prof_n+prof_g...
            -ws*taus*sg-zetaf*pxf*xf-zetan*pxn*xn-zetag*pxg*xg;
            % subsidies and profits and wages scientists cancel
%% equations
q=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- calibration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1)
q=q+1;
f(q) = MOM.FG-F/G; %Af => Af0 from LOM
%f(q) = F/E-0.8;

%3
q=q+1;
%q=q+1;
f(q) = Y - MOM.Y; % scales model!

%4
q=q+1;
f(q) = wh/wl-MOM.whwl; %=> determines ehel


q=q+1;
f(q)= MOM.EpeY*Y/pe  - E; % deltay

%2) Government := > lambdaa
q=q+1;
f(q) = - MOM.Debt + SGov;
         
%- high skill share
%3) thetag
q=q+1;
f(q) =MOM.hhg_hhghlg- (hhg/(hhg+hlg)); % MOM.hhg_hhghlg-(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag))); %=> thetag
%4)
q=q+1;
f(q) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;%=> thetan, thetaf
%5)
q=q+1;
f(q) =  thetaf-thetan;%=> thetan, thetaf

% skill market clearing (bcs whwl already determined, this determines el eh)
%9)

q=q+1;
f(q) = zh*hh*eh-(hhn+hhf+hhg); % high skill market clearing
%34

%11)
q=q+1;
f(q) =  zl*wl*el*hl-MOM.lowskill;

% chii: average hours worked
q=q+1;
f(q) = hh*zh+hl*zl-MOM.targethour; 

%emissions
q=q+1;
f(q) = omegaa - MOM.emissionsUS2019/F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1- household optimality (muu auxiliary variable determined above)
%12
q=q+1;
f(q) = chii*hh^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wh*eh)^(1-taul))+gammalh/zh*hh^taul; %=> determines hh

%13
q=q+1;
f(q) = chii*hl^(sigmaa+taul)-(muu*lambdaa*(1-taul)*(wl*el)^(1-taul))+gammall/zl*hl^taul; %=> determines hl

%3- budget
%14
q=q+1;
f(q) = zh*lambdaa*(wh*hh*eh)^(1-taul)+zl*lambdaa*(wl*hl*el)^(1-taul)+SGov-C; %=> determines C

%4- output fossil
%15
q=q+1;
f(q) = ((1-tauf)*alphaf*pf).^(alphaf/(1-alphaf))*Af.*Lf -F; 

%5- output neutral
%16
q=q+1;
f(q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
%17
q=q+1;
f(q)=  G-(Ag.*Lg).*(pg.*alphag).^(alphag./(1-alphag));

%6- demand green scientists
%18
q=q+1;
f(q)= ws-(1-alphaf)*F*pf*(1-tauf)*Af_lag/Af*gammaa*etaa*sf^(etaa-1)*rhof^(-etaa)*(A_lag/Af_lag)^phii;
%19
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))*Ag_lag/Ag/(rhog^etaa.*(1-taus));
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
% 
% %28
 q=q+1;
 f(q) = wh - thetaf*Lf/hhf.*wlf; % from optimality labour input producers fossil, and demand labour fossil
% %29
 q=q+1;
 f(q) = wl-(1-thetaf)*Lf/hlf*wlf;

% alternative with one optimality and price definition
% q=q+1;
% f(q)= wh/wl-thetaf/(1-thetaf)*hlf/hhf; % optimality fossil labour good
% 
% q=q+1;
% f(q)= wh/wl-thetag/(1-thetag)*hlg/hhg; % optimality green labour good
% 
% q=q+1;
% f(q)= wh/wl-thetan/(1-thetan)*hln/hhn; % optimality neutral labour good

% wage aggregation
% q=q+1;
% f(q) = wlf-(thetaf*wh+(1-thetaf)*wl);
% 
% q=q+1;
% f(q) = wlg-(thetag*wh+(1-thetag)*wl);
% 
% q=q+1;
% f(q) = wln-(thetan*wh+(1-thetan)*wl);

% prices and wages
%- optimality energy producers

%30
q=q+1;
f(q) = pf *F.^(1/eppse)- G.^(1/eppse).*pg; 


%- definitions prices
%31
q=q+1;
f(q) = pe - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition

%32
q=q+1;
f(q) = 1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(eppsy)); %pn - ((1-deltay.*pe.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire


%- market clearing (consumption good=> numeraire)
%34
q=q+1;
f(q) = zl*hl*el-(hln+hlf+hlg); % low skill market clearing
%35
q=q+1;
f(q) = S-(sn+sf+sg);

%36
q=q+1;
f(q) = Y-xn-xg-xf-C;

%13- Kuhn Tucker Labour supply
%36
q=q+1;
f(q)= gammalh*(upbarH-hh);
%37
q=q+1;
f(q)= gammall*(upbarH-hl);

%fprintf('number equations: %d; number variables %d', q, length(list.choiceCALIB));
end