function f = objective(y, T, params, list, Ftarget, indic, init201519)

% pars
read_in_params;

%-- transform variables

x=exp(y);
% except for taus
x(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus')))=y(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus'))) ;
% hours
if indic.noskill==0
    x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = upbarH./(1+exp(y((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
    x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = upbarH./(1+exp(y((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
else
    x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T) = upbarH./(1+exp(y((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
end
% hours scientists

if indic.target == 0
    x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
    x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
    x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));

else
    if etaa<1
        x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
        x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
        x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH*indic.minn./(1+exp(y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
    else
         x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
         x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
         x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
    end

end
%- kuhn tucker on scientists
if indic.noneutral==1
    x(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))=(y(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))).^2;    
    x(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))=(y(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))).^2;
end
% F if bounded above
if indic.target==1
    x((find(list.opt=='F')-1)*T+1+2:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1+2:find(list.opt=='F')*T)));
end
if indic.BN==1
    if indic.ineq==0
        x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T)   = B./(1+exp(y((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T)));
    else
        x((find(list.opt=='Ch')-1)*T+1:find(list.opt=='Ch')*T)   = Bh./(1+exp(y((find(list.opt=='Ch')-1)*T+1:find(list.opt=='Ch')*T)));
        x((find(list.opt=='Cl')-1)*T+1:find(list.opt=='Cl')*T)   = Bl./(1+exp(y((find(list.opt=='Cl')-1)*T+1:find(list.opt=='Cl')*T)));
    end
end

 if indic.noskill==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S, gammac]= OPT_aux_vars_notaus_flex(x, list, params, T, init201519, indic);
 else
     if indic.noneutral==0
            [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init201519, indic);
     else
          [xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag,sff, sg,  ...
            F, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pg, pf, pee,  ws, wsf, wsg,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom_nn(x, list, params, T, init201519, indic);
            sn=zeros(size(sg));
     end
 end
%% social welfare

%- create discount vector
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;

%- vector of utilities
if indic.ineq==0
    if indic.BN==0
        if thetaa~=1
            Utilcon = (C.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = log(C);
        end
    else
        Utilcon=-(C-B).^(zetaa)./(zetaa); 
    end
else
    if indic.BN==0
        if thetaa~=1
            Utilcon = zh.*(Ch.^(1-thetaa))./(1-thetaa)+(1-zh).*(Cl.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = zh.*log(Ch)+(1-zh).*log(Cl);
        end
    else
        Utilcon=zh.*(-(Ch-Bh).^(zetaa)./(zetaa))+(1-zh).*(-(Cl-Bl).^(zetaa)./(zetaa)); 
    end

end


if indic.noskill==0
    Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
    Utillab = chii.*h.^(1+sigmaa)./(1+sigmaa);
end
if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
end

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
%!! Dot product!!! so no dot.*
f = (-1)*vec_discount*(Utilcon-Utillab- Utilsci-indic.extern*weightext*(omegaa.*F).^extexpp);
end





