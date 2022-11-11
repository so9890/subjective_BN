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

if indic.xgrowth==0
 x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T) = (y((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
 x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T) = (y((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
 x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T) = (y((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
end

if indic.target==1
 x((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T)   = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T)));
end

if indic.BN==1
    if indic.ineq==0
         x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T)   = B./(1+exp(y((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T)));
    else
         x((find(list.sp=='Ch')-1)*T+1:find(list.sp=='Ch')*T)   = Bh./(1+exp(y((find(list.sp=='Ch')-1)*T+1:find(list.sp=='Ch')*T)));
         x((find(list.sp=='Cl')-1)*T+1:find(list.sp=='Cl')*T)   = Bl./(1+exp(y((find(list.sp=='Cl')-1)*T+1:find(list.sp=='Cl')*T)));
    end
end
% variables
if indic.noskill==0
    if indic.xgrowth==0
 [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil, gammac]= SP_aux_vars_2S(x, list, params, T, init, indic);
    else
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_xgrowth(x, list, params, T, init, indic);
    end
        
else
    if indic.noneutral==0
            [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill(x, list, params, T, init, indic);
    elseif indic.noneutral==1
                [ xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag, sff, sg,  ...
            F, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
             pg, pf, pee, w, wsf, wsg, tauf, taul, lambdaa,...
             wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill_noneutral(x, list, params, T, init, indic);

    end
end
% inequality constraints
c=[];
if indic.xgrowth==0
    if indic.sep==0
        c(1:T)    = S-upbarH;
    else
        
        c(1:T)=sg-upbarH;
        c(T+1:2*T)=sff-upbarH;
       
        if indic.noneutral==0
            c(2*T+1:3*T)=sn-upbarH;
        end
    end
end


% equality constraints
ceq =[];
if indic.xgrowth==0
    if indic.noneutral==0
    if indic.ineq==0
        ceq(1:1*T)     = C - (Y-xn-xg-xf);
    else
        ceq(1:1*T)     = zh.*Ch+(1-zh).*Cl - (Y-xn-xg-xf);
    end
    if indic.noneutral==0
        ceq(1*T+1:2*T) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii); 
    else
        ceq(1*T+1:2*T) = zeros(T,1);
    end
    ceq(2*T+1:3*T) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
    ceq(3*T+1:4*T) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);

    if indic.noskill==0
        ceq(4*T+1:5*T) = zh*hh - (hhn + hhf+hhg); % high skill market clearing
        ceq(5*T+1:6*T) = (1-zh)*hl - (hln + hlf+hlg);
        ceq(6*T+1:7*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
    else
        ceq(4*T+1:5*T)       = h - (Ln + Lf+Lg); % high skill market clearing
    end
    elseif indic.noneutral==1
        if indic.ineq==0
            ceq(1:1*T)     = C - (Y-xg-xf);
        else
            ceq(1:1*T)     = zh.*Ch+(1-zh).*Cl - (Y-xg-xf);
        end
       
        
        ceq(1*T+1:2*T) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
        ceq(2*T+1:3*T) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);

        if indic.noskill==0
            ceq(3*T+1:4*T) = zh*hh - ( hhf+hhg); % high skill market clearing
            ceq(4*T+1:5*T) = (1-zh)*hl - ( hlf+hlg);
            ceq(5*T+1:6*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
        else
            ceq(3*T+1:4*T)       = h - (Lf+Lg); % high skill market clearing
        end

    end
else
    if indic.noneutral==0
    if indic.ineq==0
        ceq(1:1*T)     = C - (Y-xn-xg-xf);
    else
        ceq(1:1*T)     = zh.*Ch+(1-zh).*Cl - (Y-xn-xg-xf);
    end
     if indic.noskill==0
        ceq(1*T+1:2*T) = zh*hh - (hhn + hhf+hhg); % high skill market clearing
        ceq(2*T+1:3*T) = (1-zh)*hl - (hln + hlf+hlg);
        ceq(3*T+1:4*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
    else
        ceq(1*T+1:2*T)       = h - (Ln + Lf+Lg); % high skill market clearing
     end
    elseif indic.noneutral==1
        if indic.ineq==0
            ceq(1:1*T)     = C - (Y-xg-xf);
        else
            ceq(1:1*T)     = zh.*Ch+(1-zh).*Cl - (Y-xg-xf);
        end
         if indic.noskill==0
            ceq(1*T+1:2*T) = zh*hh - ( hhf+hhg); % high skill market clearing
            ceq(2*T+1:3*T) = (1-zh)*hl - (hlf+hlg);
            ceq(3*T+1:4*T) = F-xf.^alphaf.*(Af.*Lf).^(1-alphaf);
        else
            ceq(1*T+1:2*T)       = h - (Lf+Lg); % high skill market clearing
         end
    end
end
ceq = ceq';
end