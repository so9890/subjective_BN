if length(pol)>1
    %vector version for dynamic policy
    taul    = pol(t,list.pol=='taul');
    taus    = pol(t,list.pol=='taus');
    tauf    = pol(t,list.pol=='tauf');
else
    taul    = pol(list.pol=='taul');
    taus    = pol(list.pol=='taus');
    tauf    = pol(list.pol=='tauf');
end