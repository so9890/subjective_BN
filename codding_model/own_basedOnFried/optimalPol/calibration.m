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
hhn = exp(x(list.calib=='hhn'));
hhg = exp(x(list.calib=='hhg'));
hhf = exp(x(list.calib=='hhf'));
gammalh = x(list.calib=='gammalh')^2;
gammall = x(list.calib=='gammall')^2;
C = exp(x(list.calib=='C'));
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
lambdaa = exp(x(list.calib=='lambdaa'));
deltay = 1/(1+exp(x(list.calib=='deltay')));
% Af0  = exp(x(list.calib=='Af0'));
% Ag0  = exp(x(list.calib=='Ag0'));
% An0  = exp(x(list.calib=='An0'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- auxiliary variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pf, Y, pe, pn, Etarg, N, Eopt, E, F, G, omegaa, Af, An, Ag, Lg, Ln, Lf, xf, xg, xn, ...
   SGov,  hh, hl, hhD, hlD, hln, hlg, hlf, wlg, wln, wlf,...
  wh, wl] = aux_calibFinal(deltay, eh, el, lambdaa, thetaf, thetag, thetan, C, gammall, gammalh, pg, hhn, hhf, hhg, MOM, list, paramss, poll); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

q=0;
%1)
q=q+1;
f(q)= Etarg - Eopt; % deltay

%2) Government := > lambdaa
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hhD-lambdaa.*(wh.*eh*hhD).^(1-taul))...
             +zl*(wl.*el*hlD-lambdaa.*(wl.*el*hlD).^(1-taul))...
             +tauf.*pf.*omegaa*F;         
%- high skill share
%3) thetag
q=q+1;
f(q) = MOM.hhg_hhghlg-(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag))); %=> thetag
%4)
q=q+1;
f(q) = (hhn+hhf)/(hhn+hln+hhf+hlf)-MOM.sharehighnongreen;%=> thetan, thetaf
%5)
q=q+1;
f(q) =  thetaf-thetan;%=> thetan, thetaf

% hhn hhf hhg follow from optimality labour producers as fcn of wlh wlg wlf
% => ensures high skill market clearing
%6)
q=q+1;
f(q) = hhn -thetan*Ln*wln/wh;
%7)
q=q+1;
f(q) = hhf -thetaf*Lf*wlf/wh; 
%8)
q=q+1;
f(q) = hhg -thetag*Lg*wlg/wh; 

% skill market clearing
%9)
q=q+1;
f(q) = hhD-hh; %=> eh
%10)
q=q+1;
f(q) = hlD-hl; %=> el

% goods market clearing : Holds by Walras' Law
%11)
q=q+1;
f(q) =  zl*wl*el*hlD-MOM.lowskill; 

%budget => C
q=q+1;
f(q) = zh*lambdaa*(wh*hh*eh)^(1-taul)+zl*lambdaa*(wl*hl*el)^(1-taul)+SGov-C; %Y-C-xn-xf-xg; 
%13- Kuhn Tucker Labour supply
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
end
