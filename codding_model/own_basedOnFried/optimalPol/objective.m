function f = objective(x, T, params, list)
%% read in stuff
% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% C, F, G : intermediate outputs
% Af, Ag, An : technology
% hl, hh

C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
hl     = x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T);
hh     = x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T);
 
% parameters
thetaa = params(list.params=='thetaa');
sigmaa = params(list.params=='sigmaa');
betaa = params(list.params=='betaa');
zh  = params(list.params=='zh');

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
 Utillab = (zh.*hh.^(1-sigmaa)+(1-zh).*hl.^(1-sigmaa))./(1-sigmaa);

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
%!! Dot product!!! so no dot.*
f = (-1)*vec_discount*(Utilcon-Utillab);
end





