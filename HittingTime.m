%
% NOTE: This is for calculating variance from equation generated mathematically 
%       and plots the variance with different initial v_0 and T_critical
%
% HIST:  - 10 Oct, 2020: Created by Patra
%        - 11 Nov, 2020: multiplied by time step h 
%        - 12 Jan, 2021: modified, that plots H(t) vs s
%=========================================================================
tic 
clc; close all; clear all;
ne=1;p1=200;r=0.25;s=0.1;dt=0.01;fn=[];Vs=[];
e=0.001;% epsilon
Tmax=[]; %array of the critical time T's for different initial v_0
for vini=0.005:0.001:0.05 % initial value of v_0 is changing
    t(1)=0;  %initializing x,y,z,t
    x(1)=-0.01; z(1)=0.0001;
    y(1)=vini;
    Vs=[Vs;vini];
    for k = 1:ne
        a11=-y(1);
        a12=2*p1*y(1)-x(1);
        a13=-(p1+1);
        a21=-y(1);
        a22=-x(1)-2*y(1);
        a23=0;
        a31=z(1)+s;
        a32=1;
        a33=x(1)-r;
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
            a11=-y(i+1);
            a12=p1*2*y(i+1)-x(i+1); 
            a13=-(p1+1); a21=-y(i+1);
            a22=-x(i+1)-2*y(i+1);
            a23=0;
            a31=z(i+1)+s;
            a32=1;
            a33=x(i+1)-r;
            A(:,:,i+1)=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
        end
    end
    Tcritical=t(T1+1); %% time where the v component first time touches 0.1
    Tmax=[Tmax;Tcritical];
    U_T=x(T1+1);
    sigma(:,:,1)=[0,0,0;0,0,0;0,0,0]; %% Initial sigma matrix
    I=eye(3); %% Identitty matrix
    for j=1:T1 %loop %% Eular scheme to caculate the sigma matrix
        F(:,:,j)=A(:,:,j)*sigma(:,:,j)+I+sigma(:,:,j)*(A(:,:,j))';
        sigma(:,:,j+1)=sigma(:,:,j)+F(:,:,j)*0.001;
    end
    variance= sigma(:,:,T1); %%Final sigma matrix at time t(1260) where the v component reaches the 0.1
    yvariance=variance(2,2);
    fn=[fn;(h*yvariance)/(0.01*(U_T+0.1)*(U_T+0.1))];% epsilon e is taken as 1e-3
end
figure(1);
plot(Vs, fn, '.-b', 'markersize', 10)
xlabel('v_0'); 
ylabel('H(T)'); 
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
figure(2);
plot(Tmax, fn, '.-b', 'markersize', 10)
xlabel('T (0.005<v_0<0.05)'); 
ylabel('H(T)'); 
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
toc