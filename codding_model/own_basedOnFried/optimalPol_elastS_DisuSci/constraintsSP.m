function [c, ceq] = constraintsSP(y, T, params, init, list, Ems, indic)

% pars
read_in_params;
Ftarget=(Ems+deltaa)./omegaa; 
% transform x: all are exponentially transformed
 x=exp(y);
% except for hours

if indic.noskill==0
 x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
 x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
else
 x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T) = upbarH./(1+exp(y((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T)));
end

 x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T) = (y((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
 x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T) = (y((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
 x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T) = (y((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;

if indic.target==1
 x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)   = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)));
end

% variables
if indic.noskill==0
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil, gammac]= SP_aux_vars_2S(x, list, params, T, init, indic);
else
    [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf,  taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF]= SP_aux_vars_2S_noskill(x, list, params, T, init);
end
% inequality constraints
c=[];
if indic.sep==0
    c(1:T)    = S-upbarH;
    c(T+1:2*T) = Af(T)./Af_lag(T)-1-gammac; % the growth rate in T+1 should at least be as big as from T-1 to T
 if indic.target==1
      c(2*T+1:3*T)=F-Ftarget';
 end
 
else
    c(1:T)=sn-upbarH;
    c(T+1:2*T)=sff-upbarH;
    c(2*T+1:3*T)=sg-upbarH;
    c(3*T+1:4*T) = Af(T)./Af_lag(T)-1-gammac; % the growth rate in T+1 should at least be as big as from T-1 to T

 if indic.target==1
      c(4*T+1:5*T)=F-Ftarget';
  end
end


% equality constraints
ceq =[];
ceq(1:1*T)     = C - (Y-xn-xg-xf);
ceq(1*T+1:2*T) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
ceq(2*T+1:3*T) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
ceq(3*T+1:4*T) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);

if indic.noskill==0
    ceq(4*T+1:5*T) = zh*hh - (hhn + hhf+hhg); % high skill market clearing
    ceq(5*T+1:6*T) = (1-zh)*hl - (hln + hlf+hlg);
    ceq(6*T+1:7*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
else
    ceq(4*T+1:5*T)       = h - (Ln + Lf+Lg); % high skill market clearing
end

ceq = ceq';
end