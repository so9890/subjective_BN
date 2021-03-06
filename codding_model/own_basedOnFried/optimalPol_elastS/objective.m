function f = objective(y, T, params, list, Ftarget, indic, upbS)

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
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbS./(1+exp(y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
else
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = (y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)).^2;
end
if indic.taus==1
    x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
end
% F if bounded above
if indic.target==1
    x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

% variables
if indic.noskill==0
    hl     = x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T);
    hh     = x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T);
else
    h     = x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T);
end

C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
S      = x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T);

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
    Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
    Utillab = chii.*h.^(1+sigmaa)./(1+sigmaa);
end
 Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
%!! Dot product!!! so no dot.*
f = (-1)*vec_discount*(Utilcon-Utillab- Utilsci);
end





