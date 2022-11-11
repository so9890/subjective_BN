function [ Ln, Lg, Lf, hln, hlg, hlf, xf, Ag, An, Af, A, wln, wlf, wlg, ...
        hhn, hhg, hhf, gammalh, gammall, hl, hh, wh, wl, tauf, lambdaa]=resLabTarget(x, list, pf, params,  Lnwln, Lgwlg, F, Y, AgLg, AnLn)
% params
read_in_params;

%read in vars
hhn = x(list.targlab=='hhn');
hhg = x(list.targlab=='hhg');
hhf = x(list.targlab=='hhf');
gammalh = x(list.targlab=='gammalh');
gammall = x(list.targlab=='gammall');
hl    = x(list.targlab=='HL');
hh    = x(list.targlab=='HH');
wh = (x(list.targlab=='wh'));
wl = (x(list.targlab=='wl'));
tauf = (x(list.targlab=='tauf')); % goods market clearing
lambdaa = (x(list.targlab=='lambdaa'));


xf=pf.*(1-tauf).*alphaf.*F; % includes optimal machine demand

%C= Y-xn-xf-xg; 

AfLf    = F./(alphaf*pf.*(1-tauf)).^(alphaf/(1-alphaf)); % production
Lfwlf   = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*AfLf; % optimal labour demand

hln = hhn*(1-thetan)/(thetan).*wh./wl; % hln
hlf = hhf*(1-thetaf)/(thetaf).*wh./wl; % hlf
hlg = hhg*(1-thetag)/(thetag).*wh./wl; % hlg 

% muu = C.^(-thetaa);
% SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
%             +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
%             +tauf.*pf.*F;
%         
% CHH =  zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov;% consumption determined by hh budget=> to have goods clearing

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

Ag  = AgLg/Lg; % labour good supply= labour good demand
An  = AnLn/Ln;
Af  = AfLf/Lf;
A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

wln = Lnwln/Ln;
wlg = Lgwlg/Lg;
wlf = Lfwlf/Lf;

end