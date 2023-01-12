if indic.GOV==1
    GovRev=params(list.params=='GovRev');
else
    GovRev=0;
end


if indic.sigmaWorker==0
    sigmaa = params(list.params=='sigmaa');
elseif indic.sigmaWorker==1 % more elastic labor supply
    sigmaa = 1/1.5;
elseif indic.sigmaWorker==2
    sigmaa = 1/0.5;
elseif indic.sigmaWorker==3 % inelasitc labor supply
    sigmaa = 1/0.0001;
end

sigmaas = params(list.params=='sigmaas');

if indic.util==0
    thetaa = params(list.params=='thetaa');
else
    if indic.Bop==1 % income effect dominates! => labor supply less responsive
        thetaa=(0.2+1/sigmaa)./((1-0.2)/sigmaa);
    else
        thetaa=0.4;
    end
end

if indic.sep==2
    % relative to sep==0 calibration
    wspar       = 1; %0.0683; %1; %exp(x(list.choice=='wsf'));
%     wsnpar     = 0.1; %exp(x(list.choice=='wsn'));
%     wsgpar     = 0.01; %exp(x(list.choice=='wsg'));
end
upbarH  = params(list.params=='upbarH');
upbarS  = params(list.params=='upbarS'); 
chii  = params(list.params=='chii');
zh  = params(list.params=='zh');
zs  = params(list.params=='zs');
if sum(list.params=='betaa')==1
    betaa = params(list.params=='betaa');
end
chiis  = params(list.params=='chiis');
thetaf = params(list.params=='thetaf');
thetan = params(list.params=='thetan');
thetag = params(list.params=='thetag');

if indic.labshareequ==0
    alphag = params(list.params=='alphag');
    alphaf = params(list.params=='alphaf');
    alphan = params(list.params=='alphan');
else
    alphag = 1/3*(params(list.params=='alphag')+params(list.params=='alphaf')+params(list.params=='alphan'));
    alphaf = alphag;
    alphan = alphag; 
end
if indic.subs==0
    eppsy = params(list.params=='eppsy');
else
    eppsy=1.3;
end

if indic.elasE==0
    eppse = params(list.params=='eppse');

else 
    eppse =10;
end

deltay = params(list.params=='deltay');
gammaa = params(list.params=='gammaa');
etaa = params(list.params=='etaa');
rhon = params(list.params=='rhon');

if indic.sizeequ ==0
    rhof = params(list.params=='rhof');
    rhog = params(list.params=='rhog');
else
    rhog=rhon;
    rhof=rhon;
end
if indic.noknow_spill==1
    phii=0;
elseif indic.noknow_spill ==0
    phii = 0.5;
elseif indic.noknow_spill==2
    phii =0.25;
elseif indic.noknow_spill==3
    phii =0.75; % params(list.params=='phii');
end

omegaa = params(list.params=='omegaa'); % carbon content of fossil energy
deltaa = params(list.params=='deltaa'); % natural sink

% utility externality from emissions
extexpp=1.02; 
weightext=0.01; % high weight: 

% growth rates in case of exogenous growth
vn=0;
vg=0;
vf=0;

% for exogenous growth fix research inputs at initial ss value
if indic.xgrowth==1
% hhelper= load(sprintf('params_0209_sep%d', indic.sep));
if indic.sep==0 || indic.sep==2
    sg0 =1.0305e-06; %hhelper.x0LF(hhelper.list.choice=='sg');% 
    sff0=5.5660e-08; %hhelper.x0LF(hhelper.list.choice=='sff');
    sn0= 0.3364; %hhelper.x0LF(hhelper.list.choice=='sn'); %
elseif indic.sep==1
    sg0 =0.1220;
    sn0= 0.7922;
    sff0= 0.0949; 
else
    error('initial values for scientists not given')
end
end