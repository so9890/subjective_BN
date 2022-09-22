function [c, ceq]=laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic, Emlim, t)

c=[];
if indic.sep<=1
    ceq=laissez_faire_nows_sep(x, params, list, pol, laggs, indic, Emlim, t);
elseif indic.sep==2
    ceq=laissez_faire_nows_partialS(x, params, list, pol, laggs, indic, Emlim, t);
elseif indic.sep==3
    ceq=laissez_faire_nows_sepSe(x, params, list, pol, laggs, indic, Emlim, t);
end    
end
