% This function computes the utility for the exhaustible case depending on x a vector stacking paths for the allocation of scientists, the input tax and the resource tax, the initial values for the productivity of the clean and dirty sector, the inital environmental quality and the initial stock of resource.
function U = mysimopttaxnewexhaust2(x, Ac0, Ad0, S0, Q0)
global rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim alpha2 alpha1 k1 k2 phi phi1 S_bar

A_c = zeros(numsim,1);
A_d = zeros(numsim,1);
C = zeros(numsim,1);
S = zeros(numsim,1);
Q = zeros(numsim,1);
Y_d = zeros(numsim,1);
R = zeros(numsim,1);
s_c = x(1:1:numsim);
tau = x(numsim+1:1:2*numsim);
theta = x(2*numsim+1:1:3*numsim);

%%% Initial values
S(1) = S0;
Q(1)=Q0;
A_c(1)=(1+gamma*eta_c*s_c(1))*Ac0;
A_d(1)=(1+gamma*eta_d*(1-s_c(1)))*Ad0;
D1 = (k1*A_d(1)^phi1+((alpha^alpha)*(1+tau(1))*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))^(1/phi);
D2 = k1*(1+tau(1))^epsilon*A_d(1)^phi1+(alpha^alpha*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon);
C(1) = k2*A_c(1)*A_d(1)^((1-alpha1)/(1-alpha))*((1-alpha)*(1+tau(1))^epsilon*k1*A_d(1)^phi1+(1+tau(1)-(alpha1+alpha2/(1+theta(1))))*(alpha^alpha*(cost(Q(1))*(1+theta(1)))^alpha2*A_c(1)^(1-alpha))^(1-epsilon))/(D1*D2);
Y_d(1) = k2*alpha^(-alpha*epsilon)*(cost(Q(1))*(1+theta(1)))^(-alpha2*epsilon)*A_c(1)^(1-epsilon*(1-alpha))*A_d(1)^((1-alpha1)/(1-alpha))/((D1^alpha)*D2);
R(1) = k2*alpha^(alpha*(1-epsilon))*alpha2*(cost(Q(1))*(1+theta(1)))^(alpha2*(1-epsilon)-1)*(A_c(1)^(1+phi))*A_d(1)^((1-alpha1)/(1-alpha))/(D1*D2);
%%% Simulations
for n = 2:numsim
    A_c(n)=(1+gamma*eta_c*s_c(n))*A_c(n-1);
    A_d(n)=(1+gamma*eta_d*(1-s_c(n)))*A_d(n-1);
    Q(n) = max (0,Q(n-1)-R(n-1));
    if Q(n)<=0.0001
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
        else
            C(n) = k2*A_c(n)*A_d(n)^((1-alpha1)/(1-alpha))*((1-alpha)*(1+tau(n))^epsilon*k1*A_d(n)^phi1+(1+tau(n)-(alpha1+alpha2/(1+theta(n))))*(alpha^alpha*(cost(Q(n))*(1+theta(n)))^alpha2*A_c(n)^(1-alpha))^(1-epsilon))/(D1*D2);
            Y_d(n) = k2*alpha^(-alpha*epsilon)*(cost(Q(n))*(1+theta(n)))^(-alpha2*epsilon)*A_c(n)^(1-epsilon*(1-alpha))*A_d(n)^((1-alpha1)/(1-alpha))/((D1^alpha)*D2);
        end
    end
    S(n) =  min(max(0.00000000000000001,-qsi*Y_d(n-1) + (1+delta)*S(n-1)),S_bar);
end
Teste1 = (1+repmat(rho,numsim,1));
Teste2 = Teste1.^((0:numsim-1)');
Teste3 = 1./Teste2;

Util=zeros(numsim,1);
for j=1:numsim
   Util(j)=-(1/(1-sigma))*Teste3(j)*(phiS(S(j))*C(j))^(1-sigma);
end
U=sum(Util(1:end));




