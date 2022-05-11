thetaa = params(list.params=='thetaa');
sigmaa = params(list.params=='sigmaa');
sigmaas = params(list.params=='sigmaas');

upbarH  = params(list.params=='upbarH');
chii  = params(list.params=='chii');
zh  = params(list.params=='zh');

if sum(list.params=='betaa')==1
    betaa = params(list.params=='betaa');
end
chiis  = params(list.params=='chiis');
upbS   = params(list.params=='upbS');
thetaf = params(list.params=='thetaf');
thetan = params(list.params=='thetan');
thetag = params(list.params=='thetag');

alphag = params(list.params=='alphag');
alphaf = params(list.params=='alphaf');
alphan = params(list.params=='alphan');

eppsy = params(list.params=='eppsy');
eppse = params(list.params=='eppse');
deltay = params(list.params=='deltay');
gammaa = params(list.params=='gammaa');
etaa = params(list.params=='etaa');
rhof = params(list.params=='rhof');
rhog = params(list.params=='rhog');
rhon = params(list.params=='rhon');
phii = params(list.params=='phii');

if exist('init')~=7
    Af0 = init(list.init=='Af0');
    Ag0 = init(list.init=='Ag0');
    An0 = init(list.init=='An0');
end

omegaa = params(list.params=='omegaa'); % carbon content of fossil energy
deltaa = params(list.params=='deltaa'); % natural sink