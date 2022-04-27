function [c, ceq] = constraintsSP(y, T, params, init, list, Ems, indic)

% pars
read_in_params;
Ftarget=(Ems+deltaa)./omegaa; 
% transform x: all are exponentially transformed
 x=exp(y);
% except for hours

x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
if indic.target==1
x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)   = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)));
end
% variables

[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg, ...
            F, N, G, E, Y, C, hl, hh, A_lag]= SP_aux_vars(x, list, params, T, init);

% inequality constraints
c=[];


%c(1:T) = sn +sg+ sff- S; % with this constraint it runs into singular
%matrix

if indic.target==1
    c(1:T)=F-Ftarget';
end

% equality constraints
ceq =[];
ceq(1:T)       = zh*hh - (hhn + hhf+hhg); % high skill market clearing
ceq(T+1:2*T)   = (1-zh)*hl - (hln + hlf+hlg);
ceq(2*T+1:3*T) = C - (Y-xn-xg-xf);
ceq(3*T+1:4*T) = sn-((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);
ceq(4*T+1:5*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf); 
end