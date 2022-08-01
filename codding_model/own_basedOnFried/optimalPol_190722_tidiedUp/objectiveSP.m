function f = objectiveSP(y, T, params, list, Ftarget, indic, init)

% read in stuff

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);
% % except for hours
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
    x((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T) = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T)));
end

if indic.noskill==0
    if indic.xgrowth==0
        % version below accounts for model without neutral good
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil, gammac]= SP_aux_vars_2S(x, list, params, T, init, indic);
    elseif indic.xgrowth==1
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_xgrowth(x, list, params, T, init, indic);
    end
else
    [ xn,xf,xg,Ag, An, Af,...
        Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
        F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
        pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
        wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill(x, list, params, T, init, indic);

end


%% social welfare

%- create discount vector
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;

%- vector of utilities

if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end


if indic.noskill==0
     Utillab = chii*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii*(h.^(1+sigmaa))./(1+sigmaa);
end

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
end

% continuation value
%- last period growth rate as proxy for future growth rates
gammay = Y(T)/Y(T-1)-1;
PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(T);
PVwork = 1/(1-betaa)*(Utillab(T)+Utilsci(T));
PV= betaa^T*(PVconsump-PVwork);

%Objective function value:
%!! Dot product!!! so no dot.*
% f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);
f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci-indic.extern*weightext*(omegaa.*F).^extexpp)+indic.PV*PV);


end





