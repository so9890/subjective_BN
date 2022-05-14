   function f = COMET_Objective(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip)

%choice variables over direct periods
 C = x(1:T);
 L = x(T+1:2*T);
 E = x(2*T+1:3*T);
 pi1_l = x(3*T+1:4*T);
 K1t = x(4*T+1:5*T);
 ECleanPct = x(5*T+1:6*T);
 K2t = x(6*T+1+1:7*T+1);
 sT = x(6*T+1);
 
%%% Step 1: Compute Temperature Change %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempE = E;
for i = 1:1:T;
    TempE(i) = (1-ECleanPct(i))*TempE(i);
    TempE(i) = TempE(i)+ELand(i);
end
for i = (T+1):1:(T+periods);                        
   TempE(i) = (E(T)*(1-ECleanPct(T)))*((1+gXt(T))^(i));       %Assume balanced growth path after period T
   TempE(i) = TempE(i)+ELand(i);
end
    Zt = ones(T+periods,1); 
    Zt(1) = phi23*X0+phi33*Zt0;
    Xt = ones(T+periods,1); 
    Xt(1) = phi12*S_t0+phi22*X0+phi32*Zt0;
    St = ones(T+periods,1); 
    St(1) = ((E0*10)+ELand0)+phi11*S_t0+phi21*X0;
    for i = 2:1:T+periods;
        Zt(i) = phi23*Xt(i-1)+phi33*Zt(i-1);
        Xt(i) = phi12*St(i-1)+phi22*Xt(i-1)+phi32*Zt(i-1);
        St(i) = TempE(i-1)+phi11*St(i-1)+phi21*Xt(i-1);
    end
    Qt = ones(T+periods,1); 
    Qt(1) = Qt0*(1-ksi4)+ksi4*TC0;
    Ft = ones(T+periods,1); 
    Ft(1) = (eta*((log((((St(1)+St(2))/2)+0.000001)/Sbar))/log(2)))+Fx(1); 
    TC = ones(T+periods,1);
    TC(1) = TC0+ksi1*(Ft(1)-ksi2*TC0-ksi3*(TC0-Qt0));
    for i = 2:1:T-1+periods;
        Qt(i) = Qt(i-1)*(1-ksi4)+ksi4*TC(i-1);
        Ft(i) = (eta*((log((((St(i)+St(i+1))/2)+0.000001)/Sbar))/log(2)))+Fx(i);
        TC(i) = TC(i-1)+ksi1*(Ft(i)-ksi2*TC(i-1)-ksi3*(TC(i-1)-Qt(i-1)));
    end
    m = (T-1)+periods;
        Ft(m+1) = (eta*(log(((St(m+1)+0.000001)/Sbar))/log(2)))+Fx((m+1));
        TC(m+1) = TC(m)+ksi1*(Ft(m+1)-ksi2*TC(m)-ksi3*(TC(m)-Qt(m)));
        Qt(m+1) = Qt(m)*(1-ksi4)+ksi4*TC(m);       
        
        

%%% Step 2: Compute Continuation Values of Allocations : AS FUNCTION OF LAST DIRECT PERIOD ALLOCATIONS%%%
% => SHE ASSUMES CONSTANT SHARES OF LABOUR AND CAPITAL ALLOCATED TO FINAL GOOD PRODUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Computation assumes T is sufficiently large that clean energy will be cost-competitive by then (year 2255 in DICE -> T=25 ~ 2265 ok)        

YT = ((((1+theta1*(TC(T))^2)^(-1))*(Z(T))*(((x(4*T)*x(2*T)*N(T))^(1-alpha-v))*((E(T))^(v))*((N(T)*10000*x(5*T))^alpha))));
pi1_k = x(5*T)/(x(5*T)+x(7*T+1));       %Period T share of capital in final goods production?> assumed to be constant?
Kfut = ones(periods,1);                 %Continuation aggregate capital stock, bil. int. 2005 PPP dollars
Kfut(1) = sT*(YT-Gct(T)+(1-Delta)*(N(T)*10000*(x(5*T)+x(7*T+1)))); % sonja: capital stock in first non-optimising 
                                        % period determined by savings
L(T:1:T+periods) = L(T); % labour is assumed to be constant
Yfut = zeros(periods,1);                %Continuation output, bil. int. 2005 PPP dollars
for i = 1:1:(periods-1);
  Yfut(i) = ((((1+theta1*(TC(T+i))^2)^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^i))^(1-alpha-v))*((E(T)*((1+gXt(T))^(i)))^(v))*((pi1_k*Kfut(i))^alpha))));
  Kfut(i+1) = sT*(Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i)); % tomorrow's capital stock derived from 
                                                       % resources
                                                       % constraint with
                                                       % abatement costs =0; 
                                                       % savings rate
                                                       % defined as 
                                                       % K_t+1=sT*(K_t+1+C_t)
                                                   
  C(T+i) = ((Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i))*(1-sT))/(N(T)*10000);
end
Yfut(periods) =  (((1+theta1*(TC(T+periods)^2))^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^periods))^(1-alpha-v))*(((E(T))*((1+gXt(T))^periods))^(v))*((pi1_k*Kfut(periods))^alpha)));
C(T+periods) = (Yfut(periods)-Gct(T+periods)-(Delta+gXt(T))*Kfut(periods))/(N(T)*10000);

%%% Step 3: Multiplier For Welfare Calculations%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: C(T) multiplier accounts for continuation infinite horizon utility (via line 89)
    for i = 1:1:T+periods
        C(i) = C(i)*(1+multip);
    end

%%% Step 4: Compute PDV of Utility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Composite = zeros(T+periods,1);
Util1 = zeros(T+periods,1);  % consumption and labour part of utility
UtilTC = zeros(T+periods,1); % temperature change part of utility

for i = 1:1:T+periods; % this loop runs over all periods: direct and non direct ones: 
    Composite(i) = (C(i)*((1-phi_labor*L(i))^(gamma_labor)));
    Util1(i) = (beta^(i-1))*N(i)*(((Composite(i)^(1-sigma))/(1-sigma)));
    UtilTC(i) = (beta^(i-1))*N(i)*(((1+alpha0*(TC(i)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
end

%Infinite horizon PDV of utility after (T+periods) on balanced growth path (with no population growth)
UtilTC_cont = (beta^(T+periods-1))*N(T+periods)*(1/(1-beta))*(((1+alpha0*(TC(T+periods)^alpha1))^((-1)*(1-sigma)))/(1-sigma));
Util1_cont = (beta^(T+periods-1))*N(T+periods)*(((Composite(T+periods)^(1-sigma))/(1-sigma))*(1/(1-beta*(1+gXt(T))^(1-sigma)))); 

%Objective function value:
f = (-1)*(sum(Util1)+sum(UtilTC)+Util1_cont+UtilTC_cont);





