function f= calibration(x, MOM, list, paramss, poll)

% to find
% deltay
% el
% eh 
% thetan
% thetag
% thetaf 
% lambdaa
% omegaa=> determined as auxiliary function 

% parameters
read_in_pars_calib;

% variables
pg = exp(x(list.calib=='pg'));
pn = exp(x(list.calib=='pn'));

hhn = exp(x(list.calib=='hhn'));
hhg = exp(x(list.calib=='hhg'));
hhf = exp(x(list.calib=='hhf'));
gammalh = x(list.calib=='gammalh')^2;
gammall = x(list.calib=='gammall')^2;
C = exp(x(list.calib=='C'));
hl    = upbarH/(1+exp(x(list.calib=='hl')));
hh    = upbarH/(1+exp(x(list.calib=='hh')));
% sf      = exp(x(list.calib=='sf'));
% sg      = exp(x(list.calib=='sg'));
% sn      = exp(x(list.calib=='sn'));
% ws      = exp(x(list.calib=='ws'));

% parameters
thetan = 1/(1+exp(x(list.calib=='thetan')));
thetaf = 1/(1+exp(x(list.calib=='thetaf')));
thetag = 1/(1+exp(x(list.calib=='thetag')));
el = exp(x(list.calib=='el'));
eh = exp(x(list.calib=='eh'));
chii = exp(x(list.calib=='chii'));

lambdaa = exp(x(list.calib=='lambdaa'));
deltay = 1/(1+exp(x(list.calib=='deltay')));
% Af0  = exp(x(list.calib=='Af0'));
% Ag0  = exp(x(list.calib=='Ag0'));
% An0  = exp(x(list.calib=='An0'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- auxiliary variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[muu, pf, Yout, pe, E, N, F, G, omegaa, Af, An, Ag, Lg, Ln, Lf, xf, xg, xn, ...
   SGov, hhD, hlD, hln, hlg, hlf, wlg, wln, wlf,...
  wh, wl] = aux_calibFinal(pn, hh, hl, deltay, eh, el, lambdaa, thetaf, thetag, thetan, C, pg, hhn, hhf, hhg, MOM, list, paramss, poll);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

q=0;
%1)
q=q+1;
f(q)= MOM.Y-Yout; % deltay

%2) Government := > lambdaa
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hhD-lambdaa.*(wh.*eh*hhD).^(1-taul))...
             +zl*(wl.*el*hlD-lambdaa.*(wl.*el*hlD).^(1-taul))...
             +tauf.*pf.*F;         
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

% hhn hhf hhg follow from optimality labour producers as fcn of wlh wlg wlf
% => ensures high skill market clearing
%6)
% q=q+1;
% f(q) = wln-(thetan*wh+(1-thetan)*wl);
% %7)
% q=q+1;
% f(q) = wlf-(thetaf*wh+(1-thetaf)*wl);
% assume only
q=q+1;
f(q)=(1-thetan)*Ln.*wln-wl.*hln;
% %27
q=q+1;
f(q)=(1-thetaf)*Lf.*wlf-wl.*hlf;

%8) already used to get wl!
%q=q+1;
%f(q) = hhg -thetag*Lg*wlg/wh; 
% instead use final good market clearing
q=q+1;
f(q) = Yout-xn-xg-xf-C;

% skill market clearing
%9)
q=q+1;
f(q) = hhD-hh; %=> eh
%10)
q=q+1;
f(q) = hlD-hl; %=> pg

%11)
q=q+1;
f(q) =  zl*wl*el*hlD-MOM.lowskill; 

% chii: average hours worked
q=q+1;
f(q) = hh*zh+hl*zl-MOM.targethour;  
%budget => C
q=q+1;
f(q) = zh*lambdaa*(wh*hh*eh)^(1-taul)+zl*lambdaa*(wl*hl*el)^(1-taul)+SGov-C;
% Y-C-xn-xf-xg; 

% - Model
% pn
q=q+1;
f(q)= 1-(deltay*pe^(1-eppsy)+(1-deltay)*pn^(1-eppsy))^(1/(1-eppsy));

%13-Labour supply and kuhnt tucker
%- skill supply: so that each type is indifferent how much to work
q=q+1;
f(q)= chii*hh^(sigmaa+taul)- ((muu*lambdaa*(1-taul)*(wh*eh)^(1-taul))-gammalh/zh*hh^taul); %=> determines hh
q=q+1;
f(q)= chii*hl^(sigmaa+taul) - ((muu*lambdaa*(1-taul)*(wl*el)^(1-taul))-gammall/zl*hl^taul); %=> determines hl

%12
q=q+1;
f(q)= gammalh*(upbarH-hh);
%13
q=q+1;
f(q)= gammall*(upbarH-hl);

% % research
% q=q+1;
% f(q) = Af- Af0*(1+gammaa*(sf/rhof)^etaa*(A0/Af0)^phii);
% q=q+1;
% f(q) = Ag- Ag0*(1+gammaa*(sg/rhog)^etaa*(A0/Ag0)^phii);
% q=q+1;
% f(q) = An- An0*(1+gammaa*(sn/rhon)^etaa*(A0/An0)^phii);
% % wages scientists
% q=q+1;
% f(q) = wsf-wsg;
% q=q+1;
% f(q) = wsn-wsg;
% % market clearing scientists
% q=q+1;
% f(q) = -S+sf+sn+sg; 
% q=q+1;
% f(q) = ws-wsg;
%fprintf('number equations: %d; number variables %d', q, length(list.calib));

end
