function f=aux_calib_skill(x, MOM, paramss, list, poll)
% matches equilibrium equations to 
% 1) share of high skill in green sector
% 2) share of high skill employment to total
% 3) 

read_in_pars_calib;

% variables
thetag  = 1/(1+exp(x(1))); % bounded by 0 and 1;
eleh    = exp(x(2));
thetaf   = 1/(1+exp(x(3)));
thetan   = 1/(1+exp(x(4)));

% auxiliary
hhhl   = (MOM.whwl/eleh)^((1-taul)/(taul+sigmaa)); % skill supply as fcn of wage ratio
%hlghhg = (1-thetag)/thetag*MOM.whwl; % skill demand as fcn of wage ratio
%hlfhhf = (1-thetaf)/thetaf*MOM.whwl; % skill demand as fcn of wage ratio
%hlnhhn = (1-thetan)/thetan*MOM.whwl; % skill demand as fcn of wage ratio
% HHHL   = hhhl/eleh*zh/(1-zh);

% model
% share high skill in green sector
f(1) = MOM.hhg_hhghlg-(1-(1-thetag)/(thetag/MOM.whwl+(1-thetag))); %=> determines thetag
% total skill employment
f(2) = MOM.hhehzh_total-1/(1+zh/(1-zh)*hhhl); % => determines eleh
f(3) = MOM.sharehighnongreen- 1/(1+(1-thetan)/thetan/MOM.whwl); % assuming equal shares
f(4) = thetaf-thetan;
% f(3)= hhgHH - MOM.hg_total*(1+1/HHHL)/(1+hlghhg); %=> determines share of high skill in green
% f(4)= hlgHL - MOM.hg_total*(1+HHHL)/(1+1/hlghhg); %=> determines share of high skill in green

end