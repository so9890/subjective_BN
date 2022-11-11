% This function computes the utility for the case with DTC and without clean research subsidy tax depending on the vector x of input tax, the initial values for the productivity of the clean and dirty sector and the inital environmental quality.
function U = mysimopttaxnew2onlytau(x, Ac0, Ad0, S0)
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim S_bar phi
%%% Setting vectors' sizes
A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = 0.0001*ones(numsim,1);
s_c = zeros(numsim,1);
tau = x;
if min(x)<0
    tau=max(x,0);
end
%%% Initial values
S(1) = S0;
if fint(0,tau(1),Ac0,Ad0)>=1
    s_c(1)=1;
else
    if fint(1,tau(1),Ac0,Ad0)<1
        s_c(1)=0;
    else
        s_c(1)=((eta_c/eta_d)^(1/(phi+1))*(Ac0/Ad0)^(-phi/(phi+1))*(1+tau(1))^(epsilon/(phi+1))*(1+gamma*eta_d)-1)/(gamma*(eta_c+eta_d*(eta_c/eta_d)^(1/(phi+1))*(Ac0/Ad0)^(-phi/(phi+1))*(1+tau(1))^(epsilon/(phi+1))));
    end
end
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(1)*A_d(1)/((1+tau(1))^(1-epsilon)*A_c(1)^((1-alpha)*(1-epsilon))+A_d(1)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(1)*A_c(1)^((1-alpha)*(1-epsilon))/(A_c(1)^((1-alpha)*(1-epsilon))+(1+tau(1))^epsilon*A_d(1)^((1-alpha)*(1-epsilon))));
%%% Simulations
for n = 2:numsim
    if fint(0,tau(n),A_c(n-1),A_d(n-1))>=1
        s_c(n)=1;
    else
        if fint(1,tau(n),A_c(n-1),A_d(n-1))<1
            s_c(n)=0;
        else
            s_c(n)=((eta_c/eta_d)^(1/(phi+1))*(A_c(n-1)/A_d(n-1))^(-phi/(phi+1))*(1+tau(n))^(epsilon/(phi+1))*(1+gamma*eta_d)-1)/(gamma*(eta_c+eta_d*(eta_c/eta_d)^(1/(phi+1))*(A_c(n-1)/A_d(n-1))^(-phi/(phi+1))*(1+tau(n))^(epsilon/(phi+1))));
        end
    end
    A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
    A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
    C(n)=(alpha/psi)^(alpha/(1-alpha))*A_c(n)*A_d(n)/((1+tau(n))^(1-epsilon)*A_c(n)^((1-alpha)*(1-epsilon))+A_d(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c(n)^((1-alpha)*(1-epsilon))/(A_c(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d(n)^((1-alpha)*(1-epsilon))));

      S(n) = min(max(0.00000001,-qsi*(alpha/psi)^(alpha/(1-alpha))*(A_c(n-1)^(((1-alpha)*(1-epsilon))+alpha)*A_d(n-1))/(((1+tau(n-1))^(1-epsilon)*A_c(n-1)^((1-alpha)*(1-epsilon))+A_d(n-1)^((1-alpha)*(1-epsilon)))^(alpha/((1-alpha)*(1-epsilon)))*(A_c(n-1)^((1-alpha)*(1-epsilon))+(1+tau(n-1))^epsilon*A_d(n-1)^((1-alpha)*(1-epsilon)))) + (1+delta)*S(n-1)),S_bar); %notice that S should not be very close to zero
end
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;
 
Util=zeros(numsim,1);
for j=1:numsim
   Util(j)=-(1/(1-sigma))*Teste3(j)*(phiS(S(j))*C(j))^(1-sigma);
end
U=sum(Util(1:end));



