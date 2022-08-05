function  [ceq]= calibRem_upbar_gammaa(x, MOM, list, paramss, poll, Af, An, Ag)
% this function backs out missing functions:
% 

read_in_pars_calib

% calibration to 2019 (lag = 2010-2014)
Auf_lag  = exp(x(1));
Aug_lag  = exp(x(2));
Aun_lag  = exp(x(3));
gammaaup  = exp(x(4));


% auxiliary

%- full capacity growth 
suff=1; sug=1; sun=1;
Au =  (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Au_lag  = (rhof*Auf_lag+rhon*Aun_lag+rhog*Aug_lag)/(rhof+rhon+rhog);

%% target equations


q=0;
%-- targets
q=q+1;
ceq(q) = Af- Auf_lag*(1+gammaaup*(suff/rhof)^etaa*(Au_lag/Auf_lag)^phii);
q=q+1;
ceq(q) = Ag- Aug_lag*(1+gammaaup*(sug/rhog)^etaa*(Au_lag/Aug_lag)^phii);
q=q+1;
ceq(q) = An- Aun_lag*(1+gammaaup*(sun/rhon)^etaa*(Au_lag/Aun_lag)^phii);

q=q+1;
ceq(q) = (Au/Au_lag-1)-MOM.grup;


% fprintf('number equations %d, number variables %d', q, length(x));
end