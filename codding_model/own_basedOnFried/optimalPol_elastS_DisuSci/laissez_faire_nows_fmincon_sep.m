function [c, ceq]=laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic)

c=[];
if indic.noneutral==0
    ceq=laissez_faire_nows_sep(x, params, list, pol, laggs, indic);
else
    ceq=laissez_faire_nows_sep_non(x, params, list, pol, laggs, indic);
end
end
