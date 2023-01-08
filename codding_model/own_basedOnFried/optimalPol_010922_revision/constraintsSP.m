function [c, ceq] = constraintsSP(y, T, params, init, list, Ems, indic, percon, MOM)

% pars
read_in_params;
Ftarget=(Ems+deltaa)./omegaa; 
% Ftarg_20s=(MOM.US_Budget20_30+3*deltaa)./omegaa; 
% transform x: all are exponentially transformed
 x=exp(y);
% except for hours

if indic.noskill==0
 x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
 x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
else
 x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T) = upbarH./(1+exp(y((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T)));
end

if indic.xgrowth==0
 x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T) = (y((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
 x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T) = (y((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
 x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T) = (y((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
  x((find(list.sp=='S')-1)*T+1:find(list.sp=='S')*T) = upbarS./(1+exp(y((find(list.sp=='S')-1)*T+1:find(list.sp=='S')*T)));

end

if indic.target==1
 x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)   = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
end

% variables
if indic.noskill==0
    if indic.xgrowth==0
 [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S(x, list, params, T, init, indic);
    else
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S_xgrowth(x, list, params, T, init, indic);
    end
        
else

        [ xn,xf,xg,Ag, An, Af,...
        Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
        F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
        pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
        wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S_noskill(x, list, params, T, init, indic);

end
% inequality constraints
c=[];
% c(1)= sum(F(1:percon))-Ftarg_20s; 
if indic.xgrowth==0
    if indic.sep==0
        c(1:T)    = S-upbarS;
    else
        
        c(1:T)=sg-upbarH;
        c(T+1:2*T)=sff-upbarH;
        c(2*T+1:3*T)=sn-upbarH;

    end
end


% equality constraints
ceq =[];
if indic.xgrowth==0
    ceq(1:1*T)     = C - (Y-xn-xg-xf);

    ceq(1*T+1:2*T) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii); 
    ceq(2*T+1:3*T) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
    ceq(3*T+1:4*T) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);

    if indic.noskill==0
        ceq(4*T+1:5*T) = zh*hh - (hhn + hhf+hhg); % high skill market clearing
        ceq(5*T+1:6*T) = ((1-zh))*hl - (hln + hlf+hlg);
        ceq(6*T+1:7*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
    else
        ceq(4*T+1:5*T)       =h - (Ln + Lf+Lg); % high skill market clearing
    end

else
     ceq(1:1*T)     = C - (Y-xn-xg-xf);

     if indic.noskill==0
        ceq(1*T+1:2*T) = zh*hh - (hhn + hhf+hhg); % high skill market clearing
        ceq(2*T+1:3*T) = ((1-zh))*hl - (hln + hlf+hlg);
        ceq(3*T+1:4*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
    else
        ceq(1*T+1:2*T) =h - (Ln + Lf+Lg); % high skill market clearing
     end

end
ceq = ceq';
end