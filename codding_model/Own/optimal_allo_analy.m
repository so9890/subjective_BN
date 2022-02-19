function [ybgp, xbgp, solution]= optimal_allo_analy(params, list, resultsopt, x_init)


Ac=x_init(list.x=='Ac');
Ad=x_init(list.x=='Ad');

% read in required parameters from input

thetac=params(list.params=='thetac');
thetad=params(list.params=='thetad');
sigmaa=params(list.params=='sigmaa');
zetaa=params(list.params=='zetaa');
eppsilon=params(list.params=='eppsilon');
alphaa=params(list.params=='alphaa');
psii=params(list.params=='psii');
G=params(list.params=='G');
betaa=params(list.params=='betaa');

%gammaa=params(list.params=='gammaa');
%etaa= params(list.params=='etaa');

% targets
deltaa =symms.targets(list.targets=='deltaa');
kappaa =symms.targets(list.targets=='kappaa');

Muc =symms.marginals(list.marginals=='Muc');
Muhh =symms.marginals(list.marginals=='Muhh');
%Muhl =symms.marginals(list.marginals=='Muhl');

solution_