function f = objfun_fitlogit_Abatement_Cost_Curve(x,TC_DICE,EVec,T,Pc,r)

%%% Minimize sum of squared differences between COMET and DICE abatement costs %%%
 
gamma = x(1);
b3 = x(2);
a1 = x(3);
a2 = x(4);
a3 = x(5);
a4 = x(6);
b1 = x(7);
b2 = x(8);

%%% Predicted Total Costs:
TC_hat = zeros(T,length(EVec));
diff = zeros(T,length(EVec));
diff_sq = zeros(T,length(EVec));
pv_diff_sq = zeros(T,length(EVec));
for i = 1:1:T;
    for j = 1:1:length(EVec);
      TC_hat(i,j) = ((gamma*Pc(i)*1000)/((1+(a1+a4*log(i))*exp(a2+a3*log(i)-(b1+b2*log(i))*(EVec(j)^b3)))))*EVec(j);
     diff(i,j) = TC_hat(i,j)-TC_DICE(i,j);
        diff_sq(i,j) = diff(i,j)^2;
        pv_diff_sq(i,j) = diff_sq(i,j)/((1+r)^(i-1));
    end
end

f = sum(sum(pv_diff_sq));

 TC_hat_1 = TC_hat(1,:);
 TC_1 = TC_DICE(1,:);
 
 plot(EVec,TC_hat_1,'-*b',EVec,TC_1)
