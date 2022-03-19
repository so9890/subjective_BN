function [c, ceq] = COMET_Constraints(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,b3,B0,energy_wedge_fix)

 C = x(1:T);
 L = x(T+1:2*T);
 E = x(2*T+1:3*T);
 pi1_l = x(3*T+1:4*T);
 K1t = x(4*T+1:5*T);
 ECleanPct = x(5*T+1:6*T);
 K2t = x(6*T+1+1:7*T+1);
 sT = x(6*T+1);

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 c = zeros((2*T)+1+1,1);

%%% 1. Initial capital constraint%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c((2*T)+1+1,1) = ((x(4*T+1)+x(6*T+1+1))-(K1/(N(1)*10000)));


%%% 2. Energy production constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:T;
   c(i) = (-1)*((A_E(i)*(((1-x(3*T+i))*x(T+i)*N(i))^alphaE)*((N(i)*10000*x(6*T+1+i))^(1-alphaE)))-(x(2*T+i)));
end

%%% 3. Aggregate resource constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Compute Temperature Change:
TempE = E;
for i = 1:1:T;
    TempE(i) = (1-x(5*T+i))*TempE(i);
    TempE(i) = TempE(i)+ELand(i);
end
    Zt = ones(T,1); %M_Lo
    Zt(1) = phi23*X0+phi33*Zt0;
    Xt = ones(T,1); %M_up
    Xt(1) = phi12*S_t0+phi22*X0+phi32*Zt0;
    St = ones(T,1); 
    St(1) = ((E0*10)+ELand0)+phi11*S_t0+phi21*X0;
    for i = 2:1:T;
        Zt(i) = phi23*Xt(i-1)+phi33*Zt(i-1);
        Xt(i) = phi12*St(i-1)+phi22*Xt(i-1)+phi32*Zt(i-1);
        St(i) = TempE(i-1)+phi11*St(i-1)+phi21*Xt(i-1);
    end
    Qt = ones(T,1); %T_Lo
    Qt(1) = Qt0*(1-ksi4)+ksi4*TC0;
    Ft = ones(T,1); 
    Ft(1) = (eta*((log((((St(1)+St(2))/2)+0.000001)/Sbar))/log(2)))+Fx(1); 
    TC = ones(T,1);
    TC(1) = TC0+ksi1*(Ft(1)-ksi2*TC0-ksi3*(TC0-Qt0));
    for i = 2:1:T-1;
        Qt(i) = Qt(i-1)*(1-ksi4)+ksi4*TC(i-1);
        Ft(i) = (eta*((log((((St(i)+St(i+1))/2)+0.000001)/Sbar))/log(2)))+Fx(i);
        TC(i) = TC(i-1)+ksi1*(Ft(i)-ksi2*TC(i-1)-ksi3*(TC(i-1)-Qt(i-1)));
    end
    m = (T-1);
        Ft(m+1) = (eta*(log(((St(m+1)+0.000001)/Sbar))/log(2)))+Fx((m+1));
        TC(m+1) = TC(m)+ksi1*(Ft(m+1)-ksi2*TC(m)-ksi3*(TC(m)-Qt(m)));
        Qt(m+1) = Qt(m)*(1-ksi4)+ksi4*TC(m);       

%Compute Abatement Costs & Resource Constraints:
abt_cost = zeros(T,1);
Eclean = zeros(T,1);
for j=0:1:T-2;
     Eclean(j+1) = x(5*T+1+j)*E(1+j);
     abt_cost(j+1) = ((gamma*Pc(1+j)*1000)/((1+(a1+a4*log(1+j))*exp(a2+a3*log(1+j)-(b1+b2*log(1+j))*(Eclean(1+j)^b3)))))*Eclean(1+j);
    c(T+j+1) = (-1)*((((1+theta1*(TC(1+j))^2)^(-1))*(Z(1+j))*(((x(T+1+j)*x(3*T+1+j)*N(j+1))^(1-alpha-v))*(((x(2*T+1+j)))^(v))*((N(1+j)*10000*x(4*T+1+j))^alpha)))-((x(4*T+1+1+j)+x(6*T+1+1+1+j))*10000*N(1+j+1))+(((1-delta)^10)*(x(4*T+1+j)+x(6*T+1+1+j))*10000*N(1+j))-((x(1+j)*N(1+j))*10000)-((1-trans_share)*G(j+1)*10)-abt_cost(1+j));
end
%Final Direct Optimization Period:
     Eclean(T) = x(5*T+T)*E(T);
     abt_cost(T) = ((gamma*Pc(T)*1000)/((1+(a1+a4*log(T))*exp(a2+a3*log(T)-(b1+b2*log(T))*(Eclean(T)^b3)))))*Eclean(T);  
     c(2*T) = (-1)*(((1-x(6*T+1))*(((1+theta1*(TC(T))^2)^(-1))*(Z(T))*(((x(4*T)*x(2*T)*N(T))^(1-alpha-v))*((x(3*T))^(v))*((N(T)*10000*x(5*T))^alpha))-abt_cost(T)-((1-trans_share)*G(T)*10)+(((1-delta)^10)*(x(5*T)+x(7*T+1))*10000*N(T))))-((x(T)*N(T))*10000));

     
%%% 4. Implementability constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute IMP elements:
Uct = ones(T,1);
Ult = ones(T,1);
betat = ones(T,1);
P1 = ones(T,1);
transfers = zeros(T,1);
 for i = 1:1:T;
     Uct(i) = (x(i)^(-sigma))*((1-phi_labor*x(T+i))^(gamma_labor*(1-sigma)));
     Ult(i) = (x(i)^(1-sigma))*(gamma_labor)*(-1)*(phi_labor)*((1-phi_labor*x(T+i))^(gamma_labor*(1-sigma)-1));
     betat(i) = N(i)*beta^(i-1);
     transfers(i) = (trans_share)*G(i)*10;  %billions of dollars per decade in aggregate
     transfers(i) = transfers(i)/N(i);  %dollars per decade per person (billions over billions)
     transfers(i) = transfers(i)/10000;  %tens of thousands of dollars per person per decade
         P1(i) = betat(i)*((Uct(i)*x(i))+(Ult(i)*x(T+i))-Uct(i)*transfers(i));
 end

Uc0 = Uct(1);
k0 = (K1/N(1))/10000;
r0 = (((1/(N(1)*10000*x(4*T+1)))*(alpha)*(((1+theta1*(TC(1))^2))^(-1))*(Z(1))*(((x(3*T+1)*x(T+1)*N(1))^(1-alpha-v))*((x(2*T+1))^(v))*((x(4*T+1)*10000*N(1))^alpha))));
b0 = (B0/N(1))/10000;

%IMP Constraint:
if distortionary==1
  c((2*T)+1) =(-1)*(sum(P1)-N(1)*Uc0*((k0*(1+(r0-(1-(1-delta)^10))*(1-tao_k_0)))+b0));
else
    c((2*T)+1)=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%Compute Preparatory Variables%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Yt = zeros(T,1);
  for j = 0:1:T-1
      Yt(1+j) = (((1+theta1*(TC(1+j))^2)^(-1))*(Z(1+j))*(((x(T+1+j)*x(3*T+1+j)*N(j+1))^(1-alpha-v))*(((x(2*T+1+j)))^(v))*((N(1+j)*10000*x(4*T+1+j))^alpha)));
  end 
 MPL = zeros(T,1);
 LaborTax = zeros(T,1);
 MPK = zeros(T,1);
 MRStime = ones(T,1);
 CapitalTax = ones(T,1);
 CapitalTax(1) = tao_k_0;
 Energy_Wedge = zeros(T,1);
 MPE = zeros(T,1);
 for i = 1:1:T;
     MPL(i) = ((1-alpha-v)*Yt(i))/(N(i)*pi1_l(i)*L(i));     %MPL in dollars
     LaborTax(i) = 1+((Ult(i)/Uct(i))/((MPL(i)/10000)));      
     MPK(i) = (alpha*Yt(i))/(K1t(i)*10000*N(i));
     MPE(i) = (v*Yt(i))/(E(i));    
     Energy_Wedge(i) = (MPE(i)-((MPL(i)*(1-pi1_l(i))*N(i)*L(i))/(alphaE*E(i))));
 end
 for i = 2:1:T;
   MRStime(i) = Uct(i-1)/(beta*Uct(i));
   CapitalTax(i) = 1 - ((MRStime(i)-1)/((MPK(i))-(1-(1-delta)^10)));
 end
  MAC = zeros(T,1);
  at = zeros(T,1);
  b0t = zeros(T,1);
  b1t = zeros(T,1);
  denom = zeros(T,1);
  for i = 1:1:T;
      at(i) = a1+a4*log(i);
      b0t(i) = a2+a3*log(i);
      b1t(i) = b1+b2*log(i);
      denom(i) = 1+at(i)*exp(b0t(i)-b1t(i)*(Eclean(i)^b3));
      MAC(i) = ((gamma*Pc(i)*1000)*((denom(i)^(-2))*b1t(i)*b3*(Eclean(i)^b3)*(denom(i)-1))+(gamma*Pc(i)*1000)*denom(i)^(-1));
  end
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   ceq = [];
 
%%% 1. Labor tax contraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isempty(tao_l_fix)==0
       for i = 1:1:T;   
          ceq(i) = (-1)*(((MPL(i)/10000)*(1-tao_l_fix))+(Ult(i)/Uct(i)));
      end
  end  
laenge0 = length(ceq);
  if tao_l_const==1
      for i = 1:1:T-1
          ceq(laenge0+i) = (LaborTax(1+i)-LaborTax(i));
      end
  end
   
%%% 2. Capital Tax Constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laenge1 = length(ceq);
  if isempty(tao_k_fix)==0
      for i = 1:1:T-1;
        ceq(laenge1+i) = (-1)*((1+(1-tao_k_fix)*((MPK(1+i))-(1-(1-delta)^10)))-MRStime(1+i));
     end
  end
  laenge2 = length(ceq);
  if tao_k_const==1
      for i = 1:1:T-2
          ceq(laenge2+i) = (CapitalTax(2+i) - CapitalTax(2+i-1));
      end
  end
  
%%% 3. Carbon Price Constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laenge3 = length(ceq);
  if isempty(tao_E_fix)==0
      for i = 1:1:T_tao_E_fix;
         ceq(laenge3+i) = (Energy_Wedge(i)-energy_wedge_fix(i)); 
         ceq(laenge3+T_tao_E_fix+i) = (MAC(i)-tao_E_fix(i));
      end
  end 
  
%%% 4. Intermediate Energy Input Tax Constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laenge4 = length(ceq);
 if no_interm_tax==1
      for i = 1:1:T-1
         ceq(laenge4+i) = (MAC(i)-Energy_Wedge(i));
      end
  end

ceq = ceq';