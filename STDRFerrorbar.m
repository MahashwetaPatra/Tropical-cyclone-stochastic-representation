%
% NOTE: This is for calculating standard deviation in the presence of 
%       Random Forcing that plots the standard deviation with the noise intensity a.       
%
% HIST:  - Sep, 2020: Created by Patra
%        - 04 Dec, 2020: + cleaned up, indentation and added 
%                          more notes by Patra
%=========================================================================
clc; close all; clear all;
a_array=[];
for a=0.02:0.01:0.1
    a_array=[a_array;a];
    Fraction=[];param=[];
    for N=1:10,
    % a=0.01;
    ne=100;
    temp=zeros(1,ne);
    var_tave=zeros(1,ne);
    ie = 1;
    for k = 1:ne
        u_0s = [-0.01, -1.0, -1.0, -0.1];
        v_0s = [0.01,  1.0,  1.4, -0.2];
        b_0s = [0.01,   0.5,  1.0,  0.1];
        trial=1;
        t(1)=0;  %initializing x,y,z,t
        x(1)=u_0s(trial); y(1)=v_0s(trial); z(1)=b_0s(trial);
        p=200;r=0.25;s=0.1;cd=1.0;Ts=1.0;
        dt = 0.001;
        h=0.001;   %step size
        t=0:h:30;
        f1=@(t,x,y,z) p*y*y-(p+1)*z-cd*x*y;  %ode 
        g1=@(t,x,y,z) -x*y-cd*y*y; 
        p1=@(t,x,y,z) z*x+s*x+Ts*y-r*z;
        n=(length(t)-1);
        for i=1:n %loop
            y0=y(i);
                k1=f1(t(i),x(i),y(i),z(i));
                l1=g1(t(i),x(i),y(i),z(i));
                m1=p1(t(i),x(i),y(i),z(i));
                k2=f1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
                l2=g1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
                m2=p1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
                k3=f1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                l3=g1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                m3=p1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                k4=f1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
                l4=g1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
                m4=p1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
                 
            x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6+a*sqrt(dt)*randn; %final equations
            y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6+a*sqrt(dt)*randn;
            z(i+1) = z(i) + h*(m1+2*m2 +2*m3  +m4)/6+a*sqrt(dt)*randn;
            
        end
        % set(gca, 'GridLineStyle', ':') %dotted grid lines
        % set(gca,'FontSize',14,'LineWidth',2.75)
        % figure(1);
        % plot3(x,y,z, '.','markersize',10)
        % xlabel('normalized tangential wind u');
        % ylabel('normalized vertical wind v');
        % zlabel('normalized temperature anomaly b')
        % grid on
        % view(58,9)
        % set(gca, 'GridLineStyle', ':') %dotted grid lines
        %  set(gca,'FontSize',14,'LineWidth',2.75)
        % hold on;
        % figure(2);
        % plot(t,y,'markersize',15)
        % hold on;
        var_tave(k) = std(y(n/2:n));
        if (y(n)>0)
            temp(1,ie) = y(n);
            ie = ie + 1;
        end
         
    end
    ie = ie - 1;
    v_save(1:ie)=temp(1:ie);
    var_eave = std(v_save);
    Fraction=[Fraction;var_eave];
    %fprintf('Ensemble variance: %11.5f \n',mean(var_tave));
    %fprintf('Temporal variance: %11.5f \n',var_eave);
%     figure(3);
%     plot(a,var_eave,'.b','markersize',15)
%     xlabel('Noise intensity');
%     ylabel('Standard deviation of v-wind');

%    hold on;
end
hold on;
P=mean(Fraction);
SF=std(Fraction);
ZF=1.96;
nF=10;
P1=P+(ZF*SF)/sqrt(10);
P2=P-(ZF*SF)/sqrt(10);
PB=[a;a];
EB=[P1;P2];
figure(4);
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
plot(PB,EB, '.-black', 'markersize', 10)
hold on;
plot(a,P, '.-red', 'markersize', 10)
xlabel('T (150<p<250)'); 
ylabel('H(T)'); 
end