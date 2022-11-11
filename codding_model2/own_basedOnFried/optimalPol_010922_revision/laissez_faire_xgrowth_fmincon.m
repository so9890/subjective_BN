function [c, ceq]=laissez_faire_xgrowth_fmincon(x, params, list, pol, laggs, indic, MOM,t, Emlim)

c=[];

    ceq=laissez_faire_xgrowth(x, params, list, pol, laggs, indic, MOM,t, Emlim);

end
