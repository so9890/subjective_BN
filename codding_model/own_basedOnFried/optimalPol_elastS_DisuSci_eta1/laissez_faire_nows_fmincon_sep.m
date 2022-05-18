function [c, ceq]=laissez_faire_nows_fmincon_sep(x, params, list, pol, laggs, indic)

c=[];
ceq=laissez_faire_nows_sep(x, params, list, pol, laggs, indic);
end
