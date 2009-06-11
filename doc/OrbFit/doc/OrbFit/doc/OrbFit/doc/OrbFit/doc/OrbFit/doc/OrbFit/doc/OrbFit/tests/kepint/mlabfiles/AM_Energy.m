clear;clf;
load q_coeff;
ind1=q_coeff(:,1)+1;
ind2=q_coeff(:,2)+1;
qcoe=q_coeff(:,3);
[sq1 sq2]=size(q_coeff);
qcoef=0.d0;
for m=1:sq1;
qcoef(ind1(m),ind2(m)) = qcoe(m);
end;
clear q_coeff;

load p_coeff;
ind1=p_coeff(:,1)+1;
ind2=p_coeff(:,2)+1;
pcoe=p_coeff(:,3);
[sp1 sp2]=size(p_coeff);
pcoef=0.d0;
for m=1:sp1;
pcoef(ind1(m),ind2(m)) = pcoe(m);
end;
clear p_coeff;

% grid
% Ida (243)
%rho1=1.8:0.003:2.8;
%rho2=1.8:0.003:2.8;
rho1=1.8:0.005:4.2;
rho2=1.8:0.005:4.2;

% 96217
%rho1=3:0.002:3.5;
%rho2=3:0.002:3.5;

% 1220T-2
%rho1=1.9 :0.001: 2;
%rho2=2.1 :0.001: 2.2;

% 96217
%rho1=0.15 :0.001: 0.25;
%rho2=0.25 :0.001: 0.35;

% asteroid (1) Ceres
%rho1=3 :0.005: 4;
%rho2=3 :0.005: 4;

% asteroid (100000) 
%rho1=0.8 :0.004: 2;
%rho2=0.8 :0.004: 1.8;

[vr1 vr2]=size(rho1);
[vs1 vs2]=size(rho2);

figure(1);
hold off;

for i=1:vr2;
for j=1:vs2;
AM(j,i) = 0.d0;
Energy(j,i) = 0.d0;

for h=1:3;
for k=1:3;
AM(j,i) = AM(j,i) + qcoef(h,k)*rho1(i)^(h-1)*rho2(j)^(k-1);
end;
end;
for h=1:21;
for k=1:21;
Energy(j,i) = Energy(j,i) + pcoef(h,k)*rho1(i)^(h-1)*rho2(j)^(k-1);
end;
end;

end;
end;

contour(rho1,rho2,AM,[0 0],'k');
hold on;
contour(rho1,rho2,Energy,[0 0],'r');

% Ida (243)
rho1T = 2.028773646921539
rho2T = 1.899886333422677


% 96217
%rho1T = 2.1438499950424d0;
%rho2T = 3.0340706822328d0;

%rho1T = 3.224267129284632;
%rho2T = 3.261285850669222;
%rho2T = 3.407260200297148
%rho4T = 3.424694273358647
%rho1T = 3.424694273358647
%rho2T = 2.977217781199927

%1220T-2
%rho1T = 1.9556774;
%rho2T = 2.1786839;

% (433) Eros
%rho1T = 2.0611562;
%rho2T = 1.8281632;

% (99942) Apophis
%rho1T = 0.1962117;
%rho2T = 0.2943724;

% (1) Ceres
%rho1T = 3.6433296;
%rho2T = 3.4398916;

% (100000) Astronautica
%rho1T =  1.333809493311238;
%rho1T =  1.308520836365218;
%rho1T =  1.252657040445881;
%rho2T =  1.156018454302229;
%rho1T2 = 1.580386324100532;
%rho2T2 = 1.458886141247294;
%rho1T3 = 0.9665141114006981;
%rho2T3 = 0.8919822336834913;
%rho1T4 = 0.9706698689262245;
%rho2T4 = 0.8958164670471922;
%rho1T5 = 1.911035511445894;
%rho2T5 = 1.764940703936710;

plot(rho1T,rho2T,'b*')
%plot(rho1T2,rho2T2,'b*')
%plot(rho1T3,rho2T3,'b*')
%plot(rho1T4,rho2T4,'b*')
%plot(rho1T5,rho2T5,'b*')
