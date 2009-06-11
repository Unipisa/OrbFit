clear;clf;
load q_coeff;
indq1=q_coeff(:,1)+1;
indq2=q_coeff(:,2)+1;
qcoe=q_coeff(:,3);
[sq1 sq2]=size(q_coeff);
qcoef=0.d0;
for m=1:sq1;
qcoef(indq1(m),indq2(m)) = qcoe(m);
end;
clear q_coeff;

load p_coeff;
indp1=p_coeff(:,1)+1;
indp2=p_coeff(:,2)+1;
pcoe=p_coeff(:,3);
[sp1 sp2]=size(p_coeff);
pcoef=0.d0;
for m=1:sp1;
pcoef(indp1(m),indp2(m)) = pcoe(m);
end;
clear p_coeff;

aa(1:21) = 0;
for i = 1:21;
for j = 1:21;
aa(i) = pcoef(i,j)+aa(i);
end;
end;

bb(1:3) = 0;
for i =1:3;
for j = 1:3;
bb(i) = qcoef(i,j)+bb(i);
end;
end;

for i = 1:21
SYLV(i,1) = aa(22-i);
SYLV(i+1,2) = aa(22-i);
end
for j = 3:22
SYLV(j-2,j) = bb(3);
SYLV(j-1,j) = bb(2);
SYLV(j,j) = bb(1);
end

det(SYLV);


% 96217
    rho1 = 2.1438499950424d0;
    rho2 = 3.0340706822328d0;

%1220T-2
%rho1 = 1.9556774; 
%rho2 = 2.1786839;

%plot(rho1,rho2,'b*');

