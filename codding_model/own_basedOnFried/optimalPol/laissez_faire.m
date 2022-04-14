function f=laissez_faire(x, params, list, pol, laggs, targets, Ems)
% Model
% equilibrium for one period!
% takes policy as given

%- read in policy and parameters
read_in_params;
read_in_pol;


% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% C, F, G : intermediate outputs
% Af, Ag, An : technology
% hl, hh

%- transform variables directly instead of in code
 hhf    = exp(x(list.choice=='hhf'));
 hhg    = exp(x(list.choice=='hhg'));
 hlf    = exp(x(list.choice=='hlf'));
 hlg    = exp(x(list.choice=='hlg'));
 C      = exp(x(list.choice=='C'));
 F      = exp(x(list.choice=='F'));
 G      = exp(x(list.choice=='G'));
 Af     = exp(x(list.choice=='Af'));
 Ag     = exp(x(list.choice=='Ag'));
 An     = exp(x(list.choice=='An'));
 hl     = barH-exp(x(list.choice=='hl'));
 hh     = barH-exp(x(list.choice=='hh'));
 


% auxiliary stuff
auxiliary_stuff;
%% model equations
q=0;

%1- household optimality (mu auxiliary variable determined above)
q=q+1;
f(q) = hh^(sigmaa+taul)-(mu*lambdaa*(1-taul)*wh^(1-taul)); %=> determines hh

q=q+1;
f(q) = hl^(sigmaa+taul)-(mu*lambdaa*(1-taul)*wl^(1-taul)); %=> determines hl

%3- budget
q=q+1;
f(q) = zh*lambdaa*(hh)^(1-taul)+zl*lambdaa*(hl)^(1-taul)+SGov-C; %=> determines C

%4- output fossil
q=q+1;
f(q) = ((1-tauf)*(alphaf*pf)).^(alphaf/(1-alphaf))*(Af.*Lf) -(F); 

%5- output neutral
q=q+1;
f(q) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); 

%6- demand green scientists
q=q+1;
f(q)= (1-taus)*(Ag./Ag_lag.*ws)-((gammaa*etaa*(A_lag./Ag_lag).^phii*rhog^etaa.*sg^(etaa-1).*pg.*G.*(1-alphag)));

%7- wage scientists neutral
q=q+1;
f(q)=((1-alphan)*etaa*gammaa*An_lag.^(1-phii).*A_lag.^phii.*sn.^etaa.*pn.*N)-ws.*sn.*An*rhon^etaa; 

%8- LOM neutral technology
q=q+1;
f(q) = An_lag.*(1+(sn./rhon).^etaa.*(A_lag./An_lag).^phii)-An;

%9- optimality labour input producers
q=q+1;
f(q) = thetan*Ln.*wln-wh.*hhn;
q=q+1;
f(q)= thetag*Lg.*wlg-wh.*hhg;
q=q+1;
f(q)=(1-thetan)*Ln.*wln-wl.*hln;
q=q+1;
f(q)=(1-thetag)*Lg.*wlg-wl.*hlg;

%12- market clearing consumption good=> numeraire
%q=q+1;
%f(q) = C+xf+xn+xg-Y;

%13- Kuhn Tucker Labour supply
q=q+1;
f(q)=
q=q+1;
f(q)=
%fprintf('number equations: %d; number variables %d', q, length(list.choice));
end
