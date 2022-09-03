function [c,ceq]=laissez_faireVECT_sep_fmincon(x, params, list, varrs, laggs,T, indic, Ems)
c=[];

ceq=laissez_faireVECT_sep_NoRed(x, params, list, varrs, laggs,T, indic, Ems);
end