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

if indic.BN_red==1
    B  = 0.75*params(list.params=='B');
    Bl = params(list.params=='Bl');
    Bh = 0.75* params(list.params=='Bh');
else
    B= params(list.params=='B');
    Bl = params(list.params=='Bl');
    Bh = params(list.params=='Bh');
end


zetaa = params(list.params=='zetaa');
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
if indic.noneutral==0
    deltay = params(list.params=='deltay');
else
    deltay =1; % weight on energy good
end
gammaa = params(list.params=='gammaa');
% if indic.dim==1
%     etaa =0.79; % 
% else
etaa = params(list.params=='etaa');
% end
rhof = params(list.params=='rhof');
rhog = params(list.params=='rhog');
rhon = params(list.params=='rhon');
phii = params(list.params=='phii');

% if exist('init')~=7
%     Af0 = init(list.init=='Af0');
%     Ag0 = init(list.init=='Ag0');
%     An0 = init(list.init=='An0');
% end

omegaa = params(list.params=='omegaa'); % carbon content of fossil energy
deltaa = params(list.params=='deltaa'); % natural sink

% utility externality from emissions
extexpp=1.02; 
weightext=0.01; % high weight: 

% growth rates
if indic.zero==0
    vn=0.1;
    vg=0.1;
    vf=0.1;
else
    vn=0; vg=0; vf=0; 
end