
rx=max(rI);rm=min(rI);nnr=size(rI);nr=nnr(2)
rIp=rI;dr=(rx-rm)/(nr-1);rI=rx:-dr:rm;
V4=[-.5 .5]
V2=[-2 0 2]
V1=[-1 0 1]
V5=[-5 -4 -3 -2 -1 0 1 2 3 4 5]
d10=gq+sq-g6-s6;
figure(1);
hold off
contour(ra,rI,d10,V4)
%clear d10
%pause
hold on
d13=gq-2*g6+g5;
contour(ra,rI,d13,V4)
%clear d13
'd13'
%pause
d46=sq-s6-g5+g6;
contour(ra,rI,d46,V4)
'd46'
%clear d46
pause
d2=gq-g6;
contour(ra,rI,d2,V2)
'd2'
clear d2
pause
d1=gq-g5;
contour(ra,rI,d1,V2)
'd1'
clear d1
pause
d4=sq-s6;
contour(ra,rI,d4,V2)
'd4'
clear d4
pause
z2=(gq-g6)*2+(sq-s6);
contour(ra,rI,z2,V4)
'z2'
%clear z2
pause
z3=(gq-g6)*3+(sq-s6);
contour(ra,rI,z3,V4)
'z3'
%clear z3
pause
d17=(gq-sq-g5+s6);
contour(ra,rI,d17,V4)
'd17'
%clear d17

title('Secular resonances and proper elements vers. 9.0, e=0')
xlabel('Proper a (AU)')
ylabel('Proper inclination (DEG)')
pause
load numb.fla;
N=numb(:,1);a=numb(:,3);e=numb(:,4);sI=numb(:,5);
g=numb(:,7);s=numb(:,8);
LCE=numb(:,9);H=numb(:,2);I=asin(sI)*180/pi;
clear numball;
load numbsig.fla;
sa=numbsig(:,2);se=numbsig(:,3);ssI=numbsig(:,4);
clear numbsig;
h=H<15&e>0.05&e<0.1;
f=LCE>50&h;
fs=(se>0.005|ssI>0.005)&h;
plot(a(h),I(h),'.');

pause

plot(a(fs),I(fs),'.g');
plot(a(f),I(f),'.r');


figure(2);
hold off
plot(g,s,'.')
hold on;


plot(g(fs),s(fs),'.g')

plot(g(f),s(f),'.r')

