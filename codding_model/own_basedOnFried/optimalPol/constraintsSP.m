function [c, ceq] = constraintsSP(x, T, params, init, list, Ems, indic)

% pars
read_in_params;

% transform x: all are exponentially transformed
%x=exp(y);
% except for hours

%x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)=upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
%x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)=upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));

% variables

[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg, ...
            F, N, G, E, Y, C, hl, hh, A_lag]= SP_aux_vars(x, list, params, T, init);

% inequality constraints
c=[];


c(1:T) = sn +sg+ sff- S;

if indic.target==1
c(2*T+1:3*T)=omegaa.*F-(Ems'+deltaa);
end

% equality constraints
ceq =[];
ceq(1:T)       = hh - (hhn + hhf+hhg)/zh; % high skill market clearing
ceq(T+1:2*T)   = hl - (hln + hlf+hlg)/(1-zh);
ceq(2*T+1:3*T) = C - (Y-xn-xg-xf);

end