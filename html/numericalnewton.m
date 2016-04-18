%acceleration
clear all

%
dt=0.001;
t=0:dt:2
GM = 1;  %F=m*a; F=-GmM/r^2;  a=-GM/r^2 = -GM/r^3 *x -GM/r^3 *y; 
x(1)=0.5;
vx(1)=0;
y(1)=0;
vy(1)=1.63;
r(1)=sqrt(x(1)^2 + y(1)^2);
ax(1)=-GM*x(1)/r^3;
ay(1)=-GM*y(1)/r^3;
vx(1)=vx(1)+ ax(1)*dt/2;
vy(1)= vy(1)+ay(1)*dt/2;
for i = 2:size(t,2)
    x(i)=x(i-1)+vx(i-1)*dt; y(i)=y(i-1)+vy(i-1)*dt;
    r(i)=sqrt(x(i)^2 + y(i)^2);
    ax(i)=-GM*x(i)/r(i)^3;
    ay(i)=-GM*y(i)/r(i)^3;
    vx(i)=vx(i-1)+ax(i)*dt; vy(i)=vy(i-1) + ay(i)*dt;
end
figure(1);plot(t,x)
figure(2);plot(x,y)
figure(3);plot(t,.5*(vx.^2+vy.^2)-1./r)