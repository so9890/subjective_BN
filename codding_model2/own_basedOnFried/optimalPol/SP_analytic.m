function f=SP_analytic(list, params, Ems)

% variables enter as vector special equations for last period T values
read_in_params;

% vars

hhf    = x((find(list.sp=='hhf')-1)*T+1:find(list.sp=='hhf')*T);
hhg    = x((find(list.sp=='hhg')-1)*T+1:(find(list.sp=='hhg'))*T);
hhn    = x((find(list.sp=='hhn')-1)*T+1:(find(list.sp=='hhn'))*T);
hlf    = x((find(list.sp=='hlf')-1)*T+1:find(list.sp=='hlf')*T);
hlg    = x((find(list.sp=='hlg')-1)*T+1:find(list.sp=='hlg')*T);
hln    = x((find(list.sp=='hln')-1)*T+1:find(list.sp=='hln')*T);
xn     = x((find(list.sp=='xn')-1)*T+1:find(list.sp=='xn')*T);
xf     = x((find(list.sp=='xf')-1)*T+1:find(list.sp=='xf')*T);
xg     = x((find(list.sp=='xg')-1)*T+1:find(list.sp=='xg')*T);
Af     = x((find(list.sp=='Af')-1)*T+1:find(list.sp=='Af')*T);
Ag     = x((find(list.sp=='Ag')-1)*T+1:find(list.sp=='Ag')*T);
An     = x((find(list.sp=='An')-1)*T+1:find(list.sp=='An')*T);
hl     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
hh     = x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T);
C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);
F      = x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T);
sg     = x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T);
sn      = x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T);
sff      = x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T);


% auxiliary variables!
F=(Ems'+deltaa)./omegaa;
xf=(F./(Af.*Lf).^(1-alphaf)).^(1/alphaf)M

