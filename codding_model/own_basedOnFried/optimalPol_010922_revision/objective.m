function f = objective(y, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed)

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

%     if indic.sep==1
        x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
        x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
        x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
    if indic.sep==0
%         x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
%          x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
%          x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
         x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
    end
%- kuhn tucker on scientists

% F if bounded above
if indic.target==1
    x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end


 if indic.noskill==0
          [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, taus, lambdaa,...
            wln, wlg, wlf, SWF, S, GovCon, Tls, Tlsall, PV,PVSWF, objF]= OPT_aux_vars_notaus_flex_newTauf(x, list, params, T, init201519, indic, MOM, taulFixed);
 else
           [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, GovCon, Tls, PV,PVSWF, objF]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init201519, indic, MOM);
    
 end
% %% social welfare
% 
% %- create discount vector
%  disc=repmat(betaa, 1,T);
%  expp=0:T-1;
%  vec_discount= disc.^expp;
% 
% %- vector of utilities
% 
% if thetaa~=1
%     Utilcon = (C.^(1-thetaa))./(1-thetaa);
% elseif thetaa==1
%     Utilcon = log(C);
% end
% 
% 
% if indic.noskill==0
%     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
% else
%     Utillab = chii.*h.^(1+sigmaa)./(1+sigmaa);
% end
% if indic.sep==0
%       Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
%  else
%       Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
% end
% % % continuation value
% % %- last period growth rate as proxy for future growth rates
% gammay = Y(T)/Y(T-1)-1;
% PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(T);
% PVwork = 1/(1-betaa)*(Utillab(T)+Utilsci(T));
% PV= betaa^T*(PVconsump-PVwork);
% % 
% %Objective function value:
% %!! Dot product!!! so no dot.*
%  f(1) = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci-indic.extern*weightext*(omegaa.*F).^extexpp)+indic.PV*PV);
f=-1*objF;
end





