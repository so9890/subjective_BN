function f = calibLabourAnsSc(x,  MOM, trProd, paramss, list, poll )

% to find
% deltay
% zh
% thetan
% thetag
% thetaf 
% lambdaa
% omegaa=> determined as auxiliary function 

% parameters
read_in_pars_calib
% upbarH=paramss(list.paramsdir=='upbarH');
% 
% thetaa=paramss(list.paramsdir=='thetaa');
% sigmaa=paramss(list.paramsdir=='sigmaa');
% taul=poll(list.poldir=='taul');
% tauf=poll(list.poldir=='tauf'); 

% variables

hhn = exp(x(list.calibBoth=='hhn'));
hhg = exp(x(list.calibBoth=='hhg'));
hhf = exp(x(list.calibBoth=='hhf'));
gammalh = x(list.calibBoth=='gammalh')^2;
gammall = x(list.calibBoth=='gammall')^2;
hl    = upbarH/(1+exp(x(list.calibBoth=='hl')));
hh    = upbarH/(1+exp(x(list.calibBoth=='hh')));
wh = exp(x(list.calibBoth=='wh'));
wl = exp(x(list.calibBoth=='wl'));

% parameters
thetan = 1/(1+exp(x(list.calibBoth=='thetan')));
thetaf = 1/(1+exp(x(list.calibBoth=='thetaf')));
thetag = 1/(1+exp(x(list.calibBoth=='thetag')));
zh = 1/(1+exp(x(list.calibBoth=='zh')));
chii = exp(x(list.calibBoth=='chii'));
lambdaa = exp(x(list.calibBoth=='lambdaa'));

% scientists

Af_lag  = exp(x(list.calibBoth=='Af_lag'));
Ag_lag  = exp(x(list.calibBoth=='Ag_lag'));
An_lag  = exp(x(list.calibBoth=='An_lag'));
sff     = exp(x(list.calibBoth=='sff'));
sg      = exp(x(list.calibBoth=='sg'));
sn      = exp(x(list.calibBoth=='sn'));
ws      = exp(x(list.calibBoth=='ws'));
% sigmaas = exp(x(list.calibBoth=='sigmaas'));
chiis   = exp(x(list.calibBoth=='chiis'));
gammaa  = exp(x(list.calibBoth=='gammaa'));

% AUX

[ C, Lnwln, Lgwlg, Lfwlf, pf, FF, pn, pg, ~, ~, ~, N, G,...
    ~, ~, ~,AfLf, AgLg, AnLn, ~, ~]=resProd(list, trProd, MOM , paramss, poll, 'calib'); 

hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
Lg = hhg.^thetag.*hlg.^(1-thetag);
Ln = hhn.^thetan.*hln.^(1-thetan);
Lf = hhf.^thetaf.*hlf.^(1-thetaf);
muu = C^(-thetaa);
S = MOM.targethour;
Af = AfLf/Lf;
An = AnLn/Ln;
Ag = AgLg/Lg;
A =  (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
A_lag  = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

q=0;
%1)
        
%- high skill share
%3) thetag
q=q+1;
f(q) = MOM.hhg_hhghlg-hhg/(hhg+hlg); %(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag))); %=> thetag
%4)
q=q+1;
f(q) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;%=> thetan, thetaf
%5)
q=q+1;
f(q) =  thetaf-thetan;%=> thetan, thetaf

% chii: average hours worked PER FIVE YEARS
q=q+1;
f(q) = hh*zh+hl*((1-zh))-MOM.targethour;  

% zh
q=q+1;
f(q) = MOM.whwl-wh/wl; 

%2) Government := > lambdaa $ NOTE: since tauf ==0 this is fine! otherwise
%need to be specific on policy regime! 
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
             +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
             +zs*(ws.*S-lambdaa*(ws.*S)^(1-taul))...
             +tauf.*FF;
         
%budget => C
% hhn hhf hhg follow from optimality labour producers as fcn of wlh wlg wlf

% assume only=> demand hhn, hhg, hhf
q=q+1;
f(q)=(1-thetan)*Lnwln-wl.*hln;
% %27
q=q+1;
f(q)=(1-thetaf)*Lfwlf-wl.*hlf;
% %27
q=q+1;
f(q)=(1-thetag)*Lgwlg-wl.*hlg;

% skill market clearing
%9)
q=q+1;
f(q) = hh  - (hhn+hhf+hhg)/(zh); % high skill market clearing
q=q+1;
f(q) = hl  - (hln+hlf+hlg)/(((1-zh))); % low skill market clearing

% - Model

%13-Labour supply and kuhnt tucker
%- skill supply: so that each type is indifferent how much to work
q=q+1;
f(q)= chii*hh^(sigmaa+taul)- ((muu*lambdaa*(1-taul)*(wh)^(1-taul))-gammalh/zh*hh^taul); %=> determines hh
q=q+1;
f(q)= chii*hl^(sigmaa+taul) - ((muu*lambdaa*(1-taul)*(wl)^(1-taul))-gammall/((1-zh))*hl^taul); %=> determines hl

%12
q=q+1;
f(q)= gammalh*(upbarH-hh);
%13
q=q+1;
f(q)= gammall*(upbarH-hl);

% fprintf('number equations: %d; number variables %d', q, length(list.calibBoth));

%% scientists
   q=q+1;
   f(q)= (A/A_lag-1)-MOM.growth; % targeting 5 year growth rate

% LOM => lagged technology
q=q+1;
f(q) = Af- Af_lag*(1+gammaa*(sff/rhof)^etaa*(A_lag/Af_lag)^phii);
q=q+1;
f(q) = Ag- Ag_lag*(1+gammaa*(sg/rhog)^etaa*(A_lag/Ag_lag)^phii);
 q=q+1;
 f(q) = An- An_lag*(1+gammaa*(sn/rhon)^etaa*(A_lag/An_lag)^phii);
% scientist demand
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*FF.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);
%9
q=q+1;
f(q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);

q=q+1;
f(q)= sff+sg+sn-S*zs;
q=q+1;   
f(q)= chiis.*S.^(sigmaas+taul)-muu.*lambdaa.*(1-taul).*ws.^(1-taul);
            
%f(q)= MOM.targethour-(muu*ws/(chiis)).^(1/sigmaas);  % equal disutility as for other labour => pins down ws
%  q=q+1;
%  f(q)= MOM.rhon-rhon;  % equal disutility as for other labour => pins down ws

end
