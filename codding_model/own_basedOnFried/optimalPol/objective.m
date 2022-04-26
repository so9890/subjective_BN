   function f = objective(x, T, params, list, indic)
%% read in stuff
% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% C, F, G : intermediate outputs
% Af, Ag, An : technology
% hl, hh

%  hhf    = x(1:T);
%  hhg    = x(  T+1:2*T);
%  hlf    = x(2*T+1:3*T);
%  hlg    = x(3*T+1:4*T);
 C      = x(4*T+1:5*T);
%  F      = x(5*T+1:6*T);
%  G      = x(6*T+1:7*T);
%  Af     = x(7*T+1:8*T);
%  Ag     = x(8*T+1:9*T);
%  An     = x(9*T+1:10*T);
 hl     = x(10*T+1:11*T);
 hh     = x(11*T+1:12*T);
 
 
% parameters
thetaa = params(list.params=='thetaa');
sigmaa = params(list.params=='sigmaa');
betaa = params(list.params=='betaa');
zl  = params(list.params=='zl');
zh  = params(list.params=='zh');
%% social welfare

%- create discount vector
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;

%- vector of utilities
if indic.util~=0
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
else
 Utilcon = log(C);
end
 Utillab = (zh*hh.^(1-sigmaa)+zl.*hl.^(1-sigmaa))./(1-sigmaa);

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
f = (-1)*vec_discount*(Utilcon-Utillab);
   end





