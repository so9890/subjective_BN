% This function computes the utility for the case without DTC depending on the vector x of input taxes, the initial values for the productivity of the clean and dirty sector and the inital environmental quality.
function U = mysimopttaxnew2noDTC(x, Ac0, Ad0, S0)
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim
%%% Setting vectors' sizes
A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
s_c = eta_d/(eta_c+eta_d)*ones(numsim,1);
tau = x;
q = zeros(numsim,1);

%%% Initial values
S(1) = S0;          
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(1)*A_d(1)/((1+tau(1))^(1-epsilon)*A_c(1)^((1-alpha)*(1-epsilon))+A_d(1)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(1)*A_c(1)^((1-alpha)*(1-epsilon))/(A_c(1)^((1-alpha)*(1-epsilon))+(1+tau(1))^epsilon*A_d(1)^((1-alpha)*(1-epsilon))));
%%% Simulations
for n = 2:numsim
A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
C(n)=(alpha/psi)^(alpha/(1-alpha))*A_c(n)*A_d(n)/((1+tau(n))^(1-epsilon)*A_c(n)^((1-alpha)*(1-epsilon))+A_d(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c(n)^((1-alpha)*(1-epsilon))/(A_c(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d(n)^((1-alpha)*(1-epsilon))));

  S(n) = max(0.00000000000000001,-qsi*(alpha/psi)^(alpha/(1-alpha))*(A_c(n-1)^(((1-alpha)*(1-epsilon))+alpha)*A_d(n-1))/(((1+tau(n-1))^(1-epsilon)*A_c(n-1)^((1-alpha)*(1-epsilon))+A_d(n-1)^((1-alpha)*(1-epsilon)))^(alpha/((1-alpha)*(1-epsilon)))*(A_c(n-1)^((1-alpha)*(1-epsilon))+(1+tau(n-1))^epsilon*A_d(n-1)^((1-alpha)*(1-epsilon)))) + (1+delta)*S(n-1)); %notice that S should not be very close to zero
end
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;

Util=zeros(numsim,1);
for j=1:numsim
   Util(j)=-(1/(1-sigma))*Teste3(j)*((2*(5*min(S(j),5))^.5-min(S(j),5))*C(j))^(1-sigma);
end
U=sum(Util(1:end));




