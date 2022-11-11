% This function computes the relevant parameters of the economy given a vector stacking the share of scientists in clean research and the carbon tax, the initial productivity parameters and the initial quality of the environment.
function Resp = mysimenvtaxnew2(x, Ac0, Ad0, S0)
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim phi S_bar

s_c=x(1:1:numsim);
tau=x((numsim+1):1:(2*numsim));

A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
q = zeros(numsim,1);
S(1) = S0;          
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
C(1)=(alpha/psi)^(alpha/(1-alpha))*A_c(1)*A_d(1)/((1+tau(1))^(1-epsilon)*A_c(1)^((1-alpha)*(1-epsilon))+A_d(1)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(1)*A_c(1)^((1-alpha)*(1-epsilon))/(A_c(1)^((1-alpha)*(1-epsilon))+(1+tau(1))^epsilon*A_d(1)^((1-alpha)*(1-epsilon))));

for n = 2:numsim
  A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
C(n)=(alpha/psi)^(alpha/(1-alpha))*A_c(n)*A_d(n)/((1+tau(n))^(1-epsilon)*A_c(n)^((1-alpha)*(1-epsilon))+A_d(n)^((1-alpha)*(1-epsilon)))^(1/((1-alpha)*(1-epsilon)))*(1-alpha+tau(n)*A_c(n)^((1-alpha)*(1-epsilon))/(A_c(n)^((1-alpha)*(1-epsilon))+(1+tau(n))^epsilon*A_d(n)^((1-alpha)*(1-epsilon))));

  S(n) = min(max(0.00000000000000001,-qsi*(alpha/psi)^(alpha/(1-alpha))*(A_c(n-1)^(((1-alpha)*(1-epsilon))+alpha)*A_d(n-1))/(((1+tau(n-1))^(1-epsilon)*A_c(n-1)^((1-alpha)*(1-epsilon))+A_d(n-1)^((1-alpha)*(1-epsilon)))^(alpha/((1-alpha)*(1-epsilon)))*(A_c(n-1)^((1-alpha)*(1-epsilon))+(1+tau(n-1))^epsilon*A_d(n-1)^((1-alpha)*(1-epsilon)))) + (1+delta)*S(n-1)),S_bar); %notice that S should not be very close to zero
end
% computing utility
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;
 
Util = zeros(1,numsim);
for j=1:numsim
    Util(j)=-(1/(1-sigma))*Teste3(j)*(phiS(S(j))*C(j))^(1-sigma);
end

% computing the clean research subsidy, see equation A17 for details
if  (eta_c*(1+tau(1))^epsilon*(1+gamma*eta_c*s_c(1))^(-phi-1)*Ac0^(-phi)/(eta_d*(1+gamma*eta_d*(1-s_c(1)))^(-phi-1)*Ad0^(-phi))>1 && s_c(1)<1) || (eta_c*(1+tau(1))^epsilon*(1+gamma*eta_c*s_c(1))^(-phi-1)*Ac0^(-phi)/(eta_d*(1+gamma*eta_d*(1-s_c(1)))^(-phi-1)*Ad0^(-phi))<1 && s_c(1)>0)
        q(1)=(eta_d*(1+gamma*eta_d*(1-s_c(1)))^(-phi-1)*Ad0^(-phi))/(eta_c*(1+tau(1))^epsilon*(1+gamma*eta_c*s_c(1))^(-phi-1)*Ac0^(-phi))-1;
    else
        q(1)=0;
end
for k=2:numsim
    if  (eta_c*(1+tau(k))^epsilon*(1+gamma*eta_c*s_c(k))^(-phi-1)*A_c(k-1)^(-phi)/(eta_d*(1+gamma*eta_d*(1-s_c(k)))^(-phi-1)*A_d(k-1)^(-phi))>1 && s_c(k)<1) || (eta_c*(1+tau(k))^epsilon*(1+gamma*eta_c*s_c(k))^(-phi-1)*A_c(k-1)^(-phi)/(eta_d*(1+gamma*eta_d*(1-s_c(k)))^(-phi-1)*A_d(k-1)^(-phi))<1 && s_c(k)>0)
        q(k)=(eta_d*(1+gamma*eta_d*(1-s_c(k)))^(-phi-1)*A_d(k-1)^(-phi))/(eta_c*(1+tau(k))^epsilon*(1+gamma*eta_c*s_c(k))^(-phi-1)*A_c(k-1)^(-phi))-1;
    else
        q(k)=0;
    end
end
Resp.Util = Util; % utility flow
Resp.C = C; % consumption
Resp.S = S; % environmental quality
Resp.Ac = A_c; % quality of clean machines
Resp.Ad = A_d; % quality of dirty machines
Resp.tau = tau; % input tax
Resp.S_c = s_c; % share of scientists in clean research
Resp.Q = q; % subsidy to clean research