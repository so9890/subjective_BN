% social planner problem analytical solution with target
function f= SP_target(x)


F=(Ems+deltaa)/omegaa;
read_in_params;

%- read in params
hhf
hlf
hhg
hlg
hln
hhn
xg
xn
% xf is fixed
sn
sff
sg

% auxiliary variables
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
S       = sff+sg +sn;

% loop over technology
Af=zeros(T,1);
Af_lag=[init201519(list.init=='Af0'); Af(1:T-1)];
Ag=zeros(T,1);
Ag_lag=[init201519(list.init=='Ag0'); Ag(1:T-1)];
An=zeros(T,1);
An_lag=[init201519(list.init=='An0'); An(1:T-1)];
A_lag=zeros(T,1);


for i=1:T
A_lag(i)   = (rhof*Af_lag(i)+rhon*An_lag(i)+rhog*Ag_lag(i))./(rhof+rhon+rhog);

Af(i)=Af_lag(i).*(1+gammaa*(sff(i)/rhof).^etaa.*(A_lag(i)/Af_lag(i))^phii);
Ag(i)=Ag_lag(i).*(1+gammaa*(sg(i)/rhog).^etaa.*(A_lag(i)/Ag_lag(i))^phii);
An(i)=An_lag(i).*(1+gammaa*(sn(i)/rhon).^etaa.*(A_lag(i)/An_lag(i))^phii);

%-update lags
Af_lag(i+1)=Af(i);
Ag_lag(i+1)=Ag(i);
An_lag(i+1)=An(i);

end


xf  =(F/(Af.*Lf).^(1-alphaf))^(1/alphaf);
N  = xn.^alphan.*(An.*Ln)^(1-alphan);
N       = xn.^alphan.*(An.*Ln).^(1-alphan); 
G       = xg.^alphag.*(Ag.*Lg).^(1-alphag); 

E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Y       = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% derivatives
dYdN  = Y.^(1/eppsy).*(1-deltay)^(1/eppsy).*N.^(-1/eppsy); 
dYdE  = Y.^(1/eppsy).*(deltay)^(1/eppsy).*E.^(-1/eppsy); 
dEdG  = E.^(1/eppse).*G.^(-1/eppse);
dEdF  = E.^(1/eppse).*F.^(-1/eppse);
dNdxn = (An.*Ln).^(1-alphan).*alphan.*xn.^(alphan-1);
dGdxg = (Ag.*Lg).^(1-alphag).*alphag.*xg.^(alphag-1);

dYdxg = dYdE.*dEdG.*dGdxg; 

dYdxn = dYdN.*dNdxn;

% equations
q=0;
%- xn
q=q+1;
f((q-1)*T+1:T*q) = C.^(-thetaa).*(-1+dYdxn);

%-xg
q=q+1;
f((q-1)*T+1:T*q) = C.^(-thetaa).*(-1+dYdxg);
end



