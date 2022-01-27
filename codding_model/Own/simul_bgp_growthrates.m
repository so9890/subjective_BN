function []= simul_bgp_growthrates(params, pols_num, list)

% calculates evolution of variables
% based on analytically derived 


vd=params(list.params=='vd');
vc=params(list.params=='vc');
eppsilon=params(list.params=='eppsilon');
alphaa=params(list.params=='alphaa');
tauul=pols_num(list.pol=='tauul');

1+pic= (1+v_c)