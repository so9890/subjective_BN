function [c,ceq]=laissez_faireVECT_sep_fmincon(x, params, list, varrs, laggs,T, indic)
c=[];

ceq=laissez_faireVECT_sep_NoRed(x, params, list, varrs, laggs,T, indic);
end