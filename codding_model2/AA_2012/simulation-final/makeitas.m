% This function takes as arguments the paths for the share of scientists in the clean sector and the input tax, and delivers a tax which is the maximum between the tax necessary to implement the allocation of scientists - without any subsidy - and the tax given in the argument.
function tax = makeitas(s,t)
global numsim phi eta_c eta_d epsilon gamma Ac0 Ad0
taxx=zeros(numsim,1);
Acdec=Ac0*ones(numsim,1);
Addec=Ad0*ones(numsim,1);
if s(1)==1
       taxx(1)=(eta_d/eta_c*(1+gamma*eta_d)^(-phi-1)* (Acdec(1)/Addec(1))^phi)^(1/epsilon)-1;
elseif s(1)>0
       taxx(1)=(eta_d/eta_c*(1+gamma*eta_d*(1-s(1)))^(-phi-1)*(1+gamma*eta_c*s(1))^(phi+1)*(Acdec(1)/Addec(1))^phi)^(1/epsilon)-1;
end
for i=2:numsim
    Acdec(i)=(1+gamma*eta_c*s(i-1))*Acdec(i-1);
    Addec(i)=(1+gamma*eta_d*(1-s(i-1)))*Addec(i-1);
    if s(i)==1
       taxx(i)=(eta_d/eta_c*(1+gamma*eta_d)^(-phi-1)* (Acdec(i)/Addec(i))^phi)^(1/epsilon)-1;
    elseif s(i)>0
        taxx(i)=(eta_d/eta_c*(1+gamma*eta_d*(1-s(i)))^(-phi-1)*(1+gamma*eta_c*s(i))^(phi+1)*(Acdec(i)/Addec(i))^phi)^(1/epsilon)-1;
    end
end
tax=max(1.0*taxx,t);
end