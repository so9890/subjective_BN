function EMSALL= targetSOL(list, symms, T, indexx, init201519)
% function to solve policy under target
syms taul taus pf sff sg sn real
symms.targetALL = [ taul taus pf sff sg sn];
list.targetALL = string(symms.targetALL);

helper=load('SP_target.mat');
sp_t=helper.sp_all;

% initial guess
x0= zeros(T*length(list.targetALL),1);

x0((find(list.targetALL=='pf')-1)*T+1:find(list.targetALL=='pf')*T)= log(sp_t(1:T,list.allvars=='pf'));
x0((find(list.targetALL=='taus')-1)*T+1:find(list.targetALL=='taus')*T)= sp_t(1:T,list.allvars=='taus');
x0((find(list.targetALL=='taul')-1)*T+1:find(list.targetALL=='taul')*T)= log(1-0.9); %log(1-sp_t(1:T,list.allvars=='taul'));
x0((find(list.targetALL=='sff')-1)*T+1:find(list.targetALL=='sff')*T)=log(sp_t(1:T,list.allvars=='taul'));
x0((find(list.targetALL=='sg')-1)*T+1:find(list.targetALL=='sg')*T)= log(sp_t(1:T,list.allvars=='sg'));
x0((find(list.targetALL=='sn')-1)*T+1:find(list.targetALL=='sn')*T)= log(sp_t(1:T,list.allvars=='sn'));

% laggs
laggs=init201519;
f= rd_target(x0,  list,  params, laggs, Ems, T, indexx);
rdf = @(x)rd_target(x,  list,  params, laggs, Ems, T, indexx);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sollab, fval] = fsolve(rdf, x0, options);
%
pf=
end