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

% grid for (8825)
rho1=1:0.008:2.8;
rho2=1:0.008:2.8;

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

%r1s=16;
%r2s=16;
%contour(rho1(1:r1s:vr2),rho2(1:r2s:vs2),AM(1:r1s:vr2,1:r2s:vs2),[0 0],'k--');
contour(rho1,rho2,AM,[0 0],'k');
hold on;
contour(rho1,rho2,Energy,[0 0],'r:');

% (52229)
%rho1T = 3.89832787671254;
%rho2T = 4.04283097603517;

% (8825)
% geocentric
%rho1T = 1.31600266183168;
%rho2T= 1.55321412150356
%rho2T = 1.81055223607844;

% topocentric
rho1T = 1.31600113592494;
rho2T= 1.81052775335609;

plot(rho1T,rho2T,'b*')

xlabel('\rho_1');
ylabel('\rho_2');

%title('Linkage for (8825) 1988MF, geocentric obs., \Deltat = 30 days');
title('Linkage for (8825) 1988MF, topocentric obs., \Deltat = 60 days');
