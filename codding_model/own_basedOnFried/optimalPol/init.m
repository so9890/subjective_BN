function [x0, ni, checkk]=init(Af, An, Ag, hhg, hhf, hlg, hlf, hh, hl, F, G,  gammalh, gammall, params, list,  pol, laggs, varrs, targets)

% function to get a good initial guess
% iterate over initial input variables

read_in_params;
read_in_pol;

%    read_in_pars_calib;


 checkk=10; 
 ni=0;
% while abs(checkk)>1e-5
    
%- definitions
%A       = max([Af, Ag, An]')'; 
Af_lag  = laggs(list.laggs=='Af_lag'); % shift Af backwards
Ag_lag  = laggs(list.laggs=='Ag_lag');
An_lag  = laggs(list.laggs=='An_lag');
A_lag   = max([Af_lag,Ag_lag,An_lag]);


hhn     = zh*hh-(hhf+hhg); % high skill market clearing
hln     = zl*hl-(hlf+hlg); % low skill market clearing

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

%muu      = C.^(-thetaa); % same equation in case thetaa == 1
% Muhh    = -zh*hh.^(-sigmaa);
% Muhl    = -zl*hl.^(-sigmaa);
% dIdhh   = lambdaa
% dIdhl 
 sf      = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
 sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
 sn      =  S-(sf+sg); 

E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

% prices and policy elements
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers

wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
ws      = (gammaa*etaa*(A_lag./Af_lag).^phii.*sf.^(etaa-1).*pf.*(1-tauf).*F*(1-alphaf))./(rhof^etaa.*Af./Af_lag); 


pe      = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
pn      = ((1-deltay^eppsy.*pe.^(1-eppsy))./(1-deltay)^(eppsy)).^(1/(1-eppsy)); % definition prices and numeraire

%
%taul    = (exp(wh./wl)-sigmaa*(hhhl))./(exp(hhhl)+exp(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
%lambdaa = hl.^(sigmaa+taul)./(muu.*(1-taul).*wl.^(1-taul));      % from FOC on low skill supply


% auxiliary stuff depending on prices
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +zl*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*omegaa*F;
            % subsidies and profits and wages scientists cancel
N       =  ((1-deltay)/deltay.*pe./pn)^(eppsy).*E; % demand N final good producers 
Y       =  (deltay.*E.^((eppsy-1)/eppsy)+(1-deltay).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

% wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan).*An); % price labour input neutral sector
% wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag).*Ag);


 xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
 xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
 xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;
 
 C   = Y-xn-xf-xg;
%  inc = zh*lambdaa*(hh*wh)^(1-taul)+zl*lambdaa*(hl*wl)^(1-taul)+SGov;
%  
%  %test
%  ni=ni+1;
%  if ni < 100
%  
%      checkk=C-inc;
%  
%  if checkk>1e-2
%      % then consumption too high => increase Y => Productivity; 
%      %                           => lower investment => lower prices =>
%      %                           lowr G
%      G=0.8*G; 
%      F=1.2*F;
%      An = 1.2*An; 
%  elseif checkk<1e-2
%       G=1.2*G; 
%       F=0.9*F;
%       An = 0.8*An; 
%  end
%  else
%      fprintf('did not converge')
%      break
%      
%  end

x0=eval(varrs);

x0(x0<0)=1; 

end