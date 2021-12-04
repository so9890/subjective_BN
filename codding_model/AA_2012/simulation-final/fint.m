% This is the f function of lemma 1 with an input tax. s is the share of scientists in the clean sector, t is the input tax, ac and ad are the productivity levels in the previous period
% See equation A17 in the paper.
function f = fint(s,t,ac,ad)
global eta_c eta_d phi epsilon gamma
f = eta_c/eta_d * (1+ t)^epsilon * ((1+gamma*eta_c*s)/(1+gamma*eta_d*(1-s)))^(-phi-1)*(ac/ad)^(-phi);