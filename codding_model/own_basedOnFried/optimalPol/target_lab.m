function f = target_lab(x, T, Lnwln, Lgwlg, xn, xg, AgLg, AnLn, params, list, taul,pf, F, Y)

% second function to find allocation when emission target is achieved by
% use of either taul or tauf

% indic.tauffixed == 1 then tauf is fixed as under laissez faire/BAU

% parameters
read_in_params; 

% variables

hhn = exp(x(list.targlab=='hhn'));
hhg = exp(x(list.targlab=='hhg'));
hhf = exp(x(list.targlab=='hhf'));
gammalh = x(list.targlab=='gammalh')^2;
gammall = x(list.targlab=='gammall')^2;
hl    = upbarH/(1+exp(x(list.targlab=='HL')));
hh    = upbarH/(1+exp(x(list.targlab=='HH')));
wh = exp(x(list.targlab=='wh'));
wl = exp(x(list.targlab=='wl'));
tauf = 1/(1+exp(x(list.targlab=='tauf'))); % goods market clearing
lambdaa = exp(x(list.targlab=='lambdaa'));

% AUX
xf=pf.*(1-tauf).*alphaf.*F; % includes optimal machine demand

C= Y-xn-xf-xg; 

AfLf    = F./(alphaf*pf.*(1-tauf)).^(alphaf/(1-alphaf)); % production
Lfwlf   = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*AfLf; % optimal labour demand

hln = hhn*(1-thetan)/(thetan).*wh./wl; % hln
hlf = hhf*(1-thetaf)/(thetaf).*wh./wl; % hlf
hlg = hhg*(1-thetag)/(thetag).*wh./wl; % hlg 

muu = C.^(-thetaa);
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
        
CHH =  zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov;% consumption determined by hh budget=> to have goods clearing

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

q=0;

% assume only=> demand hhn, hhg, hhf
q=q+1;
f((q-1)*T+1:T*q)=(1-thetan)*Lnwln-wl.*hln;
% %27
q=q+1;
f((q-1)*T+1:T*q)=(1-thetaf)*Lfwlf-wl.*hlf;
% %27
q=q+1;
f((q-1)*T+1:T*q)=(1-thetag)*Lgwlg-wl.*hlg;

% skill market clearing
%9)
q=q+1;
f((q-1)*T+1:T*q) = hh  - (hhn+hhf+hhg)/(zh); % high skill market clearing
q=q+1;
f((q-1)*T+1:T*q) = hl  - (hln+hlf+hlg)/((1-zh)); % low skill market clearing

% - Model

%13-Labour supply and kuhnt tucker
%- skill supply: so that each type is indifferent how much to work

q=q+1;
f((q-1)*T+1:T*q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
q=q+1;
f((q-1)*T+1:T*q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl

%12
q=q+1;
f((q-1)*T+1:T*q)= gammalh*(upbarH-hh);
%13
q=q+1;
f((q-1)*T+1:T*q)= gammall*(upbarH-hl);

q=q+1;
f((q-1)*T+1:T*q) = SGov;


q=q+1;
f((q-1)*T+1:T*q) = CHH-C;

% optimality 

% fprintf('number equations: %d; number variables %d', q, length(list.targlab));

end
