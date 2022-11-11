function f = objectiveSP_Ffixed(y, T, params, list, Ftarget, indic)

% read in stuff

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);
% % except for hours
 x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T) = upbarH./(1+exp(y((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
 x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T) = upbarH./(1+exp(y((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));

% variables
hl     = x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T);
hh     = x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T);
C      = x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T);


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
 Utillab = chii*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
%!! Dot product!!! so no dot.*
f = (-1)*vec_discount*(Utilcon-Utillab);
end





