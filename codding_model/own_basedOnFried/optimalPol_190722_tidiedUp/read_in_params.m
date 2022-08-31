if indic.util==0
    thetaa = params(list.params=='thetaa');
else
    if indic.Bop==1 % income effect dominates!
        thetaa=2;
    else
        thetaa=0.4;
    end
end
sigmaa = params(list.params=='sigmaa');
sigmaas = params(list.params=='sigmaas');

upbarH  = params(list.params=='upbarH');
chii  = params(list.params=='chii');
zh  = params(list.params=='zh');

if sum(list.params=='betaa')==1
    betaa = params(list.params=='betaa');
end
chiis  = params(list.params=='chiis');
% upbS   = params(list.params=='upbS');
thetaf = params(list.params=='thetaf');
thetan = params(list.params=='thetan');
thetag = params(list.params=='thetag');

alphag = params(list.params=='alphag');
alphaf = params(list.params=='alphaf');
alphan = params(list.params=='alphan');

if indic.subs==0
    eppsy = params(list.params=='eppsy');
else
    eppsy=1.3;
end
eppse = params(list.params=='eppse');
deltay = params(list.params=='deltay');
gammaa = params(list.params=='gammaa');
etaa = params(list.params=='etaa');
rhof = params(list.params=='rhof');
rhog = params(list.params=='rhog');
rhon = params(list.params=='rhon');
if indic.noknow_spill==1
    phii=0;
else
    phii = params(list.params=='phii');
end

omegaa = params(list.params=='omegaa'); % carbon content of fossil energy
deltaa = params(list.params=='deltaa'); % natural sink

% utility externality from emissions
extexpp=1.02; 
weightext=0.01; % high weight: 

% growth rates in case of exogenous growth
vn=0.1;
vg=0.1;
vf=0.1;
