read_in_params;

% auxiliary stuff
Af=(1+vf)*init201519(list.init=='Af0');
Ag=(1+vg)*init201519(list.init=='Ag0');

if indic.extern==0
    % no externality
    s=eppsy;
    weff=((Af*s)^eppsy*(Ag*(1-s))^(1-eppsy));
    h= ((weff)^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
else
    error('version with externality not yet coded')
    % with externality
    s0=log(eppsy);
    Uf = -weightext*extexpp*(omegaa)^extexpp*(Af*s*h).^(extexpp-1);
    ff=@x[(Af^eppsy*Ag^(1-eppsy))^(1-thetaa)*chii^(thetaa/sigmaa)/()];
end
Lf=s*h;
Lg=(1-s)*h;
F=Af*Lf;
G=Ag*Lg;
Y=(F)^(eppsy)*(G)^(1-eppsy);
C=Y;

