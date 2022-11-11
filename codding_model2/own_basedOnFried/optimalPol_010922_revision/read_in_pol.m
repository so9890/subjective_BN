
if indic.notaul~=6 % bcs then taul is an equbm object
    if indic.taul0==0
        taul    = pol(:,list.pol=='taul');
    else
        taul    = zeros(size(pol(:,list.pol=='taul')));
    end
else
    lambdaa = pol(:, list.pol =='lambdaa'); % in scenario 6 lambdaa is fixed
end

if indic.notaul~=7
    taus    = pol(:,list.pol=='taus');
end
if indic.limit_LF==0 % otherwise, tauf is choice variable
    tauf    = pol(:,list.pol=='tauf');
end