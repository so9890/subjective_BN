function model_num = model_tosolve(in1)
%MODEL_TOSOLVE
%    MODEL_NUM = MODEL_TOSOLVE(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    08-Feb-2022 12:10:38

Acp = in1(:,1);
Adp = in1(:,2);
G = in1(:,3);
H = in1(:,4);
Lc = in1(:,5);
Ld = in1(:,6);
Y = in1(:,7);
c = in1(:,8);
hh = in1(:,9);
hl = in1(:,10);
lhc = in1(:,11);
lhd = in1(:,12);
llc = in1(:,13);
lld = in1(:,14);
mu = in1(:,15);
pc = in1(:,16);
pcL = in1(:,17);
pd = in1(:,18);
pdL = in1(:,19);
wh = in1(:,20);
wl = in1(:,21);
xc = in1(:,22);
xd = in1(:,23);
yc = in1(:,24);
yd = in1(:,25);
t2 = hh.*wh;
t3 = hl.*wl;
t4 = 1.0./wh;
t5 = 1.0./wl;
t6 = sqrt(3.0);
t7 = H.^(4.0./3.0);
t8 = pc.^(3.0./5.0);
t9 = pd.^(3.0./5.0);
t10 = t2+t3;
t11 = t8+t9;
t12 = t11.^(5.0./3.0);
t13 = t10.^6.442858343500681e-17;
t14 = t10.^1.0;
model_num = [-t14+c.*t12,-mu.*t12+1.0./c,t7.*(7.0./5.0)-mu.*t13.*wh.*1.0,t7-mu.*t13.*wl.*1.0,-H+hh.*(7.0./5.0)+hl,Y-1.0./(1.0./yc.^(3.0./2.0)+1.0./yd.^(3.0./2.0)).^(2.0./3.0),yd-yc.*(pc./pd).^(2.0./5.0),yc-Lc.*sqrt(pc).*t6.*3.470865302265647e+3,pcL-pc.^(3.0./2.0).*t6.*2.313910201510431e+3,yd-Ld.*sqrt(pd).*t6.*6.941730604531284e+3,pdL-pd.^(3.0./2.0).*t6.*4.627820403020856e+3,Lc-lhc.^(7.0./1.0e+1).*llc.^(3.0./1.0e+1),Ld-lhd.^(1.4e+1./2.5e+1).*lld.^(1.1e+1./2.5e+1),lhc-llc.*(pcL.*t4.*(7.0./1.0e+1)).^(1.0e+1./3.0),llc-lhc.*(pcL.*t5.*(3.0./1.0e+1)).^(1.0e+1./7.0),lhd-lld.*(pdL.*t4.*(1.4e+1./2.5e+1)).^(2.5e+1./1.1e+1),lld-lhd.*(pdL.*t5.*(1.1e+1./2.5e+1)).^(2.5e+1./1.4e+1),Acp-3.713825873424242e+3,Adp-7.427651746848474e+3,t12-1.0,-hh+lhc+lhd,-hl+llc+lld,xd-Ld.*(pd.*3.0).^(3.0./2.0).*6.941730604531284e+3,xc-Lc.*(pc.*3.0).^(3.0./2.0).*3.470865302265647e+3,G-t10+t14];
