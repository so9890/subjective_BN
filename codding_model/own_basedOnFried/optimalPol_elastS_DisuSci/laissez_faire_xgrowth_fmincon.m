function [c, ceq]=laissez_faire_xgrowth_fmincon(x, params, list, pol, laggs, indic)

c=[];
if indic.noneutral==0
    ceq=laissez_faire_xgrowth(x, params, list, pol, laggs, indic);
else
    f=laissez_faire_xgrowth_nn(x, params, list, pol, laggs, indic);
end
end
