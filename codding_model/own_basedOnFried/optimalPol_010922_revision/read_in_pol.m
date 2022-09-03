
if indic.notaul~=6 % bcs then taul is an equbm object
    taul    = pol(:,list.pol=='taul');
else
    lambdaa = pol(:, list.pol =='lambdaa');
end

taus    = pol(:,list.pol=='taus');
if indic.limit_LF==0 % otherwise, tauf is choice variable
    tauf    = pol(:,list.pol=='tauf');
end