function f = calibLabour(x,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, F, paramss, list, poll)

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
% read_in_pars_calib;
zh=paramss(list.paramsdir=='zh');
zl=paramss(list.paramsdir=='zl');
upbarH=paramss(list.paramsdir=='upbarH');
thetaa=paramss(list.paramsdir=='thetaa');
sigmaa=paramss(list.paramsdir=='sigmaa');
taul=poll(list.poldir=='taul');
tauf=poll(list.poldir=='tauf'); 

% variables

hhn = exp(x(list.calib=='hhn'));
hhg = exp(x(list.calib=='hhg'));
hhf = exp(x(list.calib=='hhf'));
gammalh = x(list.calib=='gammalh')^2;
gammall = x(list.calib=='gammall')^2;
hl    = upbarH/(1+exp(x(list.calib=='hl')));
hh    = upbarH/(1+exp(x(list.calib=='hh')));
wh = exp(x(list.calib=='wh'));
wl = exp(x(list.calib=='wl'));

% parameters
thetan = 1/(1+exp(x(list.calib=='thetan')));
thetaf = 1/(1+exp(x(list.calib=='thetaf')));
thetag = 1/(1+exp(x(list.calib=='thetag')));
el = exp(x(list.calib=='el'));
eh = exp(x(list.calib=='eh'));
chii = exp(x(list.calib=='chii'));
lambdaa = exp(x(list.calib=='lambdaa'));

% AUX
hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
muu = C^(-thetaa);
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

%11)
%q=q+1;
%f(q) =  zl*wl*el*hl-MOM.lowskill; 
q=q+1;
f(q) = el-1;
% chii: average hours worked
q=q+1;
f(q) = hh*zh+hl*zl-MOM.targethour;  

% eleh
q=q+1;
f(q) = MOM.whwl-wh/wl; 

%2) Government := > lambdaa
q=q+1;
f(q) = - MOM.Debt + zh*(wh.*eh*hh-lambdaa.*(wh.*eh*hh).^(1-taul))...
             +zl*(wl.*el*hl-lambdaa.*(wl.*el*hl).^(1-taul))...
             +tauf.*pf.*F;
         
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
f(q) = hh  - (hhn+hhf+hhg)/(zh*eh); % high skill market clearing
q=q+1;
f(q) = hl  - (hln+hlf+hlg)/(zl*el); % low skill market clearing

% - Model

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

% fprintf('number equations: %d; number variables %d', q, length(list.calib));

end
