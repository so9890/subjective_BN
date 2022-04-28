function f = target_lab(x,  MOM, C, Lnwln, Lgwlg, Lfwlf, pf, params, list, pol)

% second function to find allocation when emission target is achieved by
% use of either taul or tauf

% indic.tauffixed == 1 then tauf is fixed as under laissez faire/BAU

% parameters
read_in_params; 

if indic.tauffixed==1
    tauf=pol(list.pold=='tauf'); 
else
    taul=pol(list.pol=='taul');
end

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
pf = exp(x(list.calib=='pf')); % follows from AfLf from previous code  

if indic.tauffixed==1
    taul = 1/(1+exp(x(list.calib=='taul')));
else
    tauf = x(list.calib=='taul');
end


% AUX

hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 
muu = C^(-thetaa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

q=0;
%1)emission target
q=q+1;
f(q)= F-(Ems'+deltaa)./omegaa;


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
f(q) = hl  - (hln+hlf+hlg)/((1-zh)); % low skill market clearing

% - Model

%13-Labour supply and kuhnt tucker
%- skill supply: so that each type is indifferent how much to work
q=q+1;
f(q)= chii*hh^(sigmaa+taul)- ((muu*lambdaa*(1-taul)*(wh)^(1-taul))-gammalh/zh*hh^taul); %=> determines hh
q=q+1;
f(q)= chii*hl^(sigmaa+taul) - ((muu*lambdaa*(1-taul)*(wl)^(1-taul))-gammall/(1-zh)*hl^taul); %=> determines hl

%12
q=q+1;
f(q)= gammalh*(upbarH-hh);
%13
q=q+1;
f(q)= gammall*(upbarH-hl);

% fprintf('number equations: %d; number variables %d', q, length(list.calib));

end
