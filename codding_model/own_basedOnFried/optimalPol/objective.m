   function f = objective(x,T,periods,params, list)

% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% N, F, G : intermediate outputs
% Af, Ag, An : technology
% Sf, Sg    => Sn follows from scientist market clearing

 hhf    = x(1:T);
 hlf    = x(  T+1:2*T);
 hlf    = x(2*T+1:3*T);
 hlg    = x(3*T+1:4*T);
 C      = x(4*T+1:5*T);
 F      = x(5*T+1:6*T);
 G      = x(6*T+1:7*T);
 Af     = x(7*T+1:8*T);
 Ag     = x(8*T+1:9*T);
 An     = x(9*T+1:10*T);
 Sf     = x(10*T+1:11*T);
 Sg     = x(11*T+1:12*T);
 hl     = x(12*T+1:13*T);
 hh     = x(13*T+1:14*T);

 
% parameters
thetaa = params(list.params=='thetaa');
sigmaa = params(list.params=='sigmaa');
betaa = params(list.params=='betaa');

%% Preliminaries


   %% could do the following if I see there is a convergence of 
   %  sector shares; otherwise periods=0 and T=12
%%% Step 2: Compute Continuation Values of Allocations : AS FUNCTION OF LAST DIRECT PERIOD ALLOCATIONS%%%
% => SHE ASSUMES CONSTANT SHARES OF LABOUR AND CAPITAL ALLOCATED TO FINAL GOOD PRODUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Computation assumes T is sufficiently large that clean energy will be cost-competitive by then (year 2255 in DICE -> T=25 ~ 2265 ok)        

% YT = ((((1+theta1*(TC(T))^2)^(-1))*(Z(T))*(((x(4*T)*x(2*T)*N(T))^(1-alpha-v))*((E(T))^(v))*((N(T)*10000*x(5*T))^alpha))));
% pi1_k = x(5*T)/(x(5*T)+x(7*T+1));       %Period T share of capital in final goods production?> assumed to be constant?
% Kfut = ones(periods,1);                 %Continuation aggregate capital stock, bil. int. 2005 PPP dollars
% Kfut(1) = sT*(YT-Gct(T)+(1-Delta)*(N(T)*10000*(x(5*T)+x(7*T+1)))); % sonja: capital stock in first non-optimising 
%                                         period determined by savings
% L(T:1:T+periods) = L(T);                % labour supply is constant
% Yfut = zeros(periods,1);                % Continuation output, bil. int. 2005 PPP dollars
% for i = 1:1:(periods-1);
%   Yfut(i) = ((((1+theta1*(TC(T+i))^2)^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^i))^(1-alpha-v))...
%       *((E(T)*((1+gXt(T))^(i)))^(v))*((pi1_k*Kfut(i))^alpha))));
%   Kfut(i+1) = sT*(Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i));  
%   C(T+i) = ((Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i))*(1-sT))/(N(T)*10000);
% end
% Yfut(periods) =  (((1+theta1*(TC(T+periods)^2))^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^periods))^(1-alpha-v))*(((E(T))*((1+gXt(T))^periods))^(v))*((pi1_k*Kfut(periods))^alpha)));
% C(T+periods) = (Yfut(periods)-Gct(T+periods)-(Delta+gXt(T))*Kfut(periods))/(N(T)*10000);

%%% Step 3: Multiplier For Welfare Calculations%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: C(T) multiplier accounts for continuation infinite horizon utility (via line 89)
%     for i = 1:1:T+periods
%         C(i) = C(i)*(1+multip);
%     end

%%% Step 4: Compute PDV of Utility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Utilcon   = zeros(T+periods,1);
Utillab   = zeros(T+periods,1);
for i = 1:1:T+periods
   
    Utilcon(i) = (betaa^(i-1))*(C(i)^(1-thetaa))/(1-thetaa);
    Utillab(i) = (betaa^(i-1))*(zh*hh(i)^(1-sigmaa)+zl*hl(i)^(1-sigmaa))/(1-sigmaa);
end

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
% UtilTC_cont = (betaa^(T+periods-1))*N(T+periods)*(1/(1-betaa))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
% Util1_cont = (betaa^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-betaa*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
f = (-1)*(sum(Utilcon)-sum(Utillab));





