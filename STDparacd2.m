%
% NOTE: This is for calculating standard deviation in the presence of 
%       Random parametric noise when added to parameter cd throughout the 
%       time evolution and STD vs noise amplitude is plotted. 95% Confidence
%       Interval is added with the mean STD
%
% HIST:  - March, 2020: Created by Patra
%=========================================================================
tic
clc; close all; clear all;
P_array=zeros(9,1);
a_array=zeros(9,1);
exp=1;
for a=0.01:0.01:0.09
    Fraction=zeros(10,1);
    Vmax_array=zeros(10,1);
    for N=1:10
        ne=100;
        temp=zeros(1,ne);
        var_tave=zeros(1,ne);
        ie = 1;
        for k = 1:ne
            u_0s = [-0.01, -1.0, -1.0, -0.1];
            v_0s = [0.01,  1.0,  1.4, -0.2];
            b_0s = [0.01,   0.5,  1.0,  0.1];
            trial=1;
            dt = 0.001;
            h=0.001;   %step size
            t=0:h:30;
            n=(length(t)-1);
            x=zeros(1,length(t));
            y=zeros(1,length(t));
            z=zeros(1,length(t));
            x(1)=u_0s(trial); y(1)=v_0s(trial); z(1)=b_0s(trial);
            p=200;r=0.25;s=0.1;cd=1.0; Ts=1.0;
            %cd=1.0+a*sqrt(dt)*randn;
            cd = min(max(1+(a*randn),0.5),2.5);
            %Ts = min(max(1+(a*randn),0.5),2.0);
            %s = min(max(0.1+(a*randn),0.0),1.0);
            f1=@(t,x,y,z,cd,Ts,s) p*y*y-(p+1)*z-cd*x*y;  %ode 
            g1=@(t,x,y,z,cd,Ts,s) -x*y-cd*y*y; 
            p1=@(t,x,y,z,cd,Ts,s) z*x+s*x+Ts*y-r*z;
            for i=1:n %loop
                %cd = min(max(1+(a*randn),0.5),2.5);
                %Ts = min(max(1+(a*randn),0.5),2.0);
                %s = min(max(0.1+(a*randn),0.0),1.0);
                k1=f1(t(i),x(i),y(i),z(i),cd,Ts,s);
                l1=g1(t(i),x(i),y(i),z(i),cd,Ts,s);
                m1=p1(t(i),x(i),y(i),z(i),cd,Ts,s);
                k2=f1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),cd,Ts,s);
                l2=g1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),cd,Ts,s);
                m2=p1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)),cd,Ts,s);
                k3=f1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),cd,Ts,s);
                l3=g1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),cd,Ts,s);
                m3=p1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)),cd,Ts,s);
                k4=f1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),cd,Ts,s);
                l4=g1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),cd,Ts,s);
                m4=p1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h),cd,Ts,s);
                x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6; %final equations
                y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6;
                z(i+1) = z(i) + h*(m1+2*m2 +2*m3  +m4)/6;       
            end
            var_tave(k) = std(y(n/2:n));
            if (y(n)>0)
                temp(1,ie) = y(n);
                ie = ie + 1;
            end
        end
        ie = ie - 1;
        v_save(1:ie)=temp(1:ie);
        v_avg=mean(v_save);
        var_eave = std(v_save);
        Vmax_array(N)=v_avg;
        Fraction(N)=var_eave;
    end
    P=mean(Fraction);
    P_array(exp)=P;
    a_array(exp)=a;
    exp=exp+1;
    SF=std(Fraction);
    ZF=1.96;
    nF=10;
    P1=P+(ZF*SF)/sqrt(10);
    P2=P-(ZF*SF)/sqrt(10);
    PB=[a;a];
    EB=[P1;P2];
    figure(1);
    set(gca, 'GridLineStyle', ':') %dotted grid lines
    set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
    plot(PB,EB, '.-black', 'markersize', 10)
    hold on;
    xlabel('T (150<p<250)'); 
    ylabel('H(T)'); 
end
hold on;
plot(a_array,P_array, '.red', 'markersize', 20)    
toc