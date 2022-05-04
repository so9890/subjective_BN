function [c, ceq] = constraintsSP_Ffixed(y, T, params, init, list, Ems, indic)

% pars
read_in_params;
Ftarget=(Ems+deltaa)./omegaa; 
F=Ftarget';
% transform x: all are exponentially transformed
 x=exp(y);
% except for hours

x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));

% variables

[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            N, G, E, Y, C, hl, hh, A_lag]= ...
            SP_aux_vars_2S_Ffixed(x, list, params, T, init, F);
% inequality constraints
c=[];

% equality constraints
ceq =[];
ceq(1:T)       = zh*hh - (hhn + hhf+hhg); % high skill market clearing
ceq(T+1:2*T)   = (1-zh)*hl - (hln + hlf+hlg);
ceq(2*T+1:3*T) = C - (Y-xn-xg-xf);
% ceq(3*T+1:4*T) = sn-((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);
ceq(3*T+1:4*T) = S*2/3-(sn); 
ceq(4*T+1:5*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf); % determines xf
ceq(5*T+1:6*T) = S*1/3-(sg+sff);
ceq = ceq';
end