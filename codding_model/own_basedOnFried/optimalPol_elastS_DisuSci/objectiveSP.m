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
    x((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T) = Ftarget'./(1+exp(y((find(list.sp=='F')-1)*T+1:find(list.sp=='F')*T)));
end

if indic.noskill==0
    [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsn, wsf,  taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S(x, list, params, T, init);
else
   [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf,  taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF]= SP_aux_vars_2S_noskill(x, list, params, T, init);
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
 Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
%!! Dot product!!! so no dot.*
f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);
end





