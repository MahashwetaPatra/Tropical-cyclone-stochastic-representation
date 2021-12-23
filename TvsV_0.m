 clc; close all; clear all;
 ne=1;p1=200;r=0.25;s=0.1;dt=0.01;fn=[];v_0s=[];Tc=[];
 for trial=0.004:0.005:0.05,% trial is the initial value of v_0
t(1)=0;  %initializing x,y,z,t
x(1)=-0.01; z(1)=0.0001;
y(1)=trial;
v_0s=[v_0s;trial];
 for k = 1:ne

a11=-y(1); a12=2*p1*y(1)-x(1); a13=-(p1+1); a21=-y(1); a22=-x(1)-2*y(1);  a23=0; a31=z(1)+s;  a32=1; a33=x(1)-r;
     A(:,:,1)=[a11, a12, a13; a21, a22, a23; a31, a32, a33];

h=0.001;   %step size
t=0:h:30;
f=@(t,x,y,z) p1*y*y-(p1+1)*z-x*y;  %ode 
g=@(t,x,y,z) -x*y-y*y; 
p=@(t,x,y,z) z*x+s*x+y-r*z;
for i=1:(length(t)-1) %loop
    k1=f(t(i),x(i),y(i),z(i));
    l1=g(t(i),x(i),y(i),z(i));
    m1=p(t(i),x(i),y(i),z(i));
      k2=f(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));     
      l2=g(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
      m2=p(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
      k3=f(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
      l3=g(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
      m3=p(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
      k4=f(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
      l4=g(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
      m4=p(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
       x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6; %final equations
      y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6;
      z(i+1) = z(i) + h*(m1+2*m2 +2*m3  +m4)/6;

      check=y(i+1);
      if check<0.1
         T1= i+1;
      end
a11=-y(i+1); a12=p1*2*y(i+1)-x(i+1); a13=-(p1+1); a21=-y(i+1); a22=-x(i+1)-2*y(i+1); a23=0; a31=z(i+1)+s; a32=1; a33=x(i+1)-r;
A(:,:,i+1)=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
end
   set(gca, 'GridLineStyle', ':') %dotted grid lines
        set(gca,'FontSize',14,'LineWidth',2.75)

figure(1);
plot3(x,y,z, '.','markersize',10)
 xlabel('normalized tangential wind u'); 
        ylabel('normalized vertical wind v'); 
        zlabel('normalized temperature anomaly b')
        grid on
        view(58,9)
        set(gca, 'GridLineStyle', ':') %dotted grid lines
        set(gca,'FontSize',14,'LineWidth',2.75)
hold on;
figure(2);
plot(t,y,'markersize',15)
hold on;

  
 end

 t(T1+1); %% time where the v component first time touches 0.1
 Tc=[Tc;t(T1+1)];

%  figure(3);
%  plot(trial, t(T1+1), '.b')
%  hold on;
%    xlabel('v_0'); 
%         ylabel('T'); 
%     set(gca, 'GridLineStyle', ':') %dotted grid lines
%         set(gca,'FontSize',14,'LineWidth',2.75)

 end
  figure(3);
%  plot(trial, t(T1+1), '.b')
%  hold on;
   
 plot(v_0s,Tc, '.-b', 'markersize', 10)
  xlabel('v_0'); 
         ylabel('T'); 
     set(gca, 'GridLineStyle', ':') %dotted grid lines
         set(gca,'FontSize',14,'LineWidth',2.75)