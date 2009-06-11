clear;
load coeff;
co00=coeff(1,1);
co01=coeff(1,2);
co02=coeff(1,3);
co10=coeff(1,4);
co11=coeff(1,5);
co20=coeff(1,6);
co21=coeff(1,7);
co30=coeff(1,8);
co40=coeff(1,9);
clear coeff;
load admis_reg;
c0=admis_reg(1);
c1=admis_reg(2);
c2=admis_reg(3);
c3=admis_reg(4);
c4=admis_reg(5);
c5=admis_reg(6);
k2=admis_reg(7);
clear admis_reg;

r1=0:0.005:10;
rd1=-0.5:0.002:0.3;
[s1 s2]=size(r1);
[t1 t2]=size(rd1);

f3=co30;
f4=co40;
for j=1:t2;
  f0(j)=co00+co01*rd1(j)+co02*rd1(j)^2;
  f1(j)=co10+co11*rd1(j);
  f2(j)=co20+co21*rd1(j);
  for i=1:s2;
    f(j,i)=f0(j)+f1(j)*r1(i)+f2(j)*r1(i)^2+f3*r1(i)^3+f4*r1(i)^4;
  end;
end;

for i=1:s2;
  Cr = c2*r1(i)^2 + c3*r1(i) + c4;
  Sr = r1(i)^2 + c5*r1(i) + c0;
  for j=1:t2;
%    dEEarth(j,i) = rd1(j)^2 + r1(i)^2*c2 - 2*k2*muE/r1(i)  ;
    dESun(j,i) = rd1(j)^2 + c1*rd1(j) + Cr - 2*k2/sqrt(Sr);
  end;
end;

r1true=input('r1=');
%r1true = 1.78331154669379;
rd1true=input('rd1=');
%rd1true = -2.023782577794471E-004;
r2true = 2.0379380;
rd2true = 0.0101617;
figure(1);
hold off;
contour(r1,rd1,f,[0 0]);
hold on;
plot(r1true,rd1true,'r*');
%plot(r2true,rd2true,'g*');
%contour(r1,rd1,dEEarth,[0 0],'b');
contour(r1,rd1,dESun,[0 0],'k');



