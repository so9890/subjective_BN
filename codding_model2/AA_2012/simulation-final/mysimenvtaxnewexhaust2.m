% This function computes the relevant parameters of the economy given a vector stacking the share of scientists in clean research, the carbon tax and the resource tax, the initial productivity parameters, the initial quality of the environment and the initial stock of resource.
function Resp = mysimenvtaxnewexhaust2(x, Ac0, Ad0, S0, Q0)
global kappa rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim alpha2 alpha1 k1 k2 phi1 phi S_bar

s_c = x(1:1:numsim);
tau = x(numsim+1:1:2*numsim);
theta = x(2*numsim+1:1:3*numsim);

%%% Setting vectors' sizes
A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
Q = zeros(numsim,1);
Y_d = zeros(numsim,1);
R = zeros(numsim,1);
cleansubs = zeros(numsim,1);
profitratio=zeros(numsim,1);
S(1) = S0; 
Q(1) = Q0;
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
D1 = (k1*A_d(1)^phi1+((alpha^alpha)*(1+tau(1))*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))^(1/phi);
D2 = k1*(1+tau(1))^epsilon*A_d(1)^phi1+(alpha^alpha*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon);
C(1) = k2*A_c(1)*A_d(1)^((1-alpha1)/(1-alpha))*((1-alpha)*(1+tau(1))^epsilon*k1*A_d(1)^phi1+(1+tau(1)-(alpha1+alpha2/(1+theta(1))))*(alpha^alpha*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))/(D1*D2);
Y_d(1) = k2*alpha^(-alpha*epsilon)*(cost(Q(1))*(1+theta(1)))^(-alpha2*epsilon)*A_c(1)^(1-epsilon*(1-alpha))*A_d(1)^((1-alpha1)/(1-alpha))/((D1^alpha)*D2);
R(1) = k2*alpha^(alpha*(1-epsilon))*alpha2*(cost(Q(1))*(1+theta(1)))^(alpha2*(1-epsilon)-1)*(A_c(1)^(1+phi))*A_d(1)^((1-alpha1)/(1-alpha))/(D1*D2);
profitratio(1)=eta_d*A_c(1)^phi/(eta_c*kappa*(1+tau(1))^epsilon*(cost(Q(1))*(1+theta(1)))^(alpha2*(epsilon-1))*A_d(1)^phi1)*(1+gamma*eta_c*s_c(1))/(1+gamma*eta_d*(1-s_c(1)));
%%% Simulations
for n = 2:numsim
    A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
    A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
    Q(n) = max (0,Q(n-1)-R(n-1));
    if Q(n)==0
        R(n)=0;
        Y_d(n)=0;
        C(n) = (1-alpha)/alpha^(alpha/(1-alpha))*A_c(n);
    else
        D1 = (k1*A_d(n)^phi1+((alpha^alpha)*(1+tau(n))*(cost(Q(n))*(1+theta(n)))^alpha2*A_c(n)^(1-alpha))^(1-epsilon))^(1/phi);
        D2 = k1*(1+tau(n))^epsilon*A_d(n)^phi1+(alpha^alpha*(cost(Q(n))*(1+theta(n)))^alpha2*A_c(n)^(1-alpha))^(1-epsilon);
        R(n) = k2*alpha^(alpha*(1-epsilon))*alpha2*(cost(Q(n))*(1+theta(n)))^(alpha2*(1-epsilon)-1)*(A_c(n)^(1+phi))*A_d(n)^((1-alpha1)/(1-alpha))/(D1*D2);
        if R(n)>= Q(n)
            R(n)=Q(n);
            pd0 = alpha^alpha*cost(Q(n))^alpha2*A_c(n)^(1-alpha)/(k1*A_d(n)^phi1+ (alpha^alpha*(1+tau(n))*cost(Q(n))^alpha2)^(1-epsilon)*A_c(n)^phi)^(1/(1-epsilon));
            pd = fsolve (@ (pd) psi^(alpha1/(1-alpha1)-alpha/(1-alpha))*alpha^(alpha/(1-alpha))*(1 - pd^(1-epsilon)*(1+tau(n))^(1-epsilon))^(1/phi)*A_c(n)/(alpha1^(alpha1/(1-alpha1))*R(n)^(alpha2/(1-alpha1))*(1+tau(n))^(epsilon*alpha2/(1-alpha1))*(1-tau(n)*(1+tau(n))^(-epsilon)*pd^(1-epsilon))^(alpha2/(1-alpha1))*pd^((1-alpha2*(1-epsilon))/(1-alpha1))*A_d(n)) - 1, pd0);
            pc = (1-((1+tau(n))*pd)^(1-epsilon))^(1/(1-epsilon));
            C(n) = 1/((1+tau(n))^epsilon*pc^(1-epsilon)+pd^(1-epsilon))*((((alpha1/psi)^(alpha1/(1-alpha1))*R(n)^(alpha2/(1-alpha1))*pd^((phi+alpha1)/(1-alpha1))*A_d(n)*((1+tau(n))^epsilon*pc^(1-epsilon)+pd^(1-epsilon))^(alpha2/(1-alpha1)))^((epsilon-1)/epsilon)+((alpha/psi)^(alpha/(1-alpha))*(1+tau(n))^epsilon*pc^(1/(1-alpha)-epsilon)*A_c(n))^((epsilon-1)/epsilon))^(epsilon/(epsilon-1))-(alpha1^(1/(1-alpha1))*psi^(-1/(1-alpha1))*R(n)^(alpha2/(1-alpha1))*pd^((1+phi)/(1-alpha1))*((1+tau(n))^epsilon*pc^(1-epsilon)+pd^(1-epsilon))^(alpha2/(1-alpha1))*A_d(n)+alpha^(1/(1-alpha1))*psi^(-1/(1-alpha))*(1+tau(n))^epsilon*pc^((1+phi)/(1-alpha))*A_c(n)))-cost(Q(n))*R(n);
            Y_d(n) = 1/((1+tau(n))^epsilon*pc^(1-epsilon)+pd^(1-epsilon))*(alpha1/psi)^(alpha1/(1-alpha1))*R(n)^(alpha2/(1-alpha1))*pd^((phi+alpha1)/(1-alpha1))*A_d(n)*((1+tau(n))^epsilon*pc^(1-epsilon)+pd^(1-epsilon))^(alpha2/(1-alpha1));
            profitratio(n) = eta_d*(1-alpha1)/(eta_c*(1-alpha)*(1+tau(n))^epsilon)*(pd/pc)^(1-epsilon)*(1+gamma*eta_c*s_c(n))/(1+gamma*eta_d*(1-s_c(n)));
        else
            C(n) = k2*A_c(n)*A_d(n)^((1-alpha1)/(1-alpha))*((1-alpha)*(1+tau(n))^epsilon*k1*A_d(n)^phi1+(1+tau(n)-(alpha1+alpha2/(1+theta(n))))*(alpha^alpha*(cost(Q(n))*(1+theta(n)))^alpha2*A_c(n)^(1-alpha))^(1-epsilon))/(D1*D2);
            Y_d(n) = k2*alpha^(-alpha*epsilon)*(cost(Q(n))*(1+theta(n)))^(-alpha2*epsilon)*A_c(n)^(1-epsilon*(1-alpha))*A_d(n)^((1-alpha1)/(1-alpha))/((D1^alpha)*D2);
            profitratio(n)=eta_d*A_c(n)^phi/(eta_c*kappa*(1+tau(n))^epsilon*(cost(Q(n))*(1+theta(n)))^(alpha2*(epsilon-1))*A_d(n)^phi1)*(1+gamma*eta_c*s_c(n))/(1+gamma*eta_d*(1-s_c(n)));
        end
    end
    S(n) =  min(max(0.00000000000000001,-qsi*Y_d(n-1) + (1+delta)*S(n-1)),S_bar);
end

Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;
Util = zeros(1,numsim);
for j=1:numsim
    Util(j)=-(1/(1-sigma))*Teste3(j)*(phiS(S(j))*C(j))^(1-sigma);
end
for k=1:numsim
    if s_c(k)==0 && profitratio(k)>1
        cleansubs(k)=0;
    else
        if s_c(k)==1 && profitratio(k)<1
        cleansubs(k)=0;
        else
        cleansubs(k)= profitratio(k)-1;
        end
    end
end
Resp.Util = Util; % utility flow
Resp.C = C; % consumption
Resp.S = S; % quality of the environment
Resp.Ac = A_c; % productivity of the clean sector
Resp.Ad = A_d; % productivity of the dirty sector
Resp.tau = tau; % input tax
Resp.S_c = s_c; % share of scientists in clean research
Resp.Q = Q; % stock of resource
Resp.theta = theta; % resource tax
Resp.R = R; % resource extraction
Resp.Y_d =Y_d; % production of dirty input
Resp.cleansubs=cleansubs;% subsidy to clean research