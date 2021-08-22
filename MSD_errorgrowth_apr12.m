%
% NOTE: This is for calculating error growth in the presence of 
%       Random Forcing that plots the error growth with tau.       
%
% HIST:  - 05 Nov, 2020: Created by Patra
%        - 11 Nov, 2020: + introduced nef member, cleaned up, and added 
%                          more notes by CK
%                        + set a control orbit to be a deterministic one
%                          by removing random forcing, and set nef=1
%        - 12 Apr, 2021: + Update Cd, Ts parameters by CK
%                        + Optimize the code to speed up workflow
%                        + Add options for other random methods   
%=========================================================================
clc; close all; clear all;
%
% Random forcing options
%
eps_RF=0.0;                % random var added to model forcing
eps_IC=0.000;              % random var added to ini. condition at time tau
eps_Cd=0.1;                 % random var added to Cd parameter
eps_Ts=0.0;                 % random var added to Ts parameter
eps_Ss=0.0;                % random var added to static stability s parameter
%
% Model parameters
%
p=200;r=0.25;s=0.1;         % (p,r,s) parameters
h=0.001;                    % time step dt
dt=h;                       % time step again??? need to clean up
t=0:h:32;                   % Integration interval
n=(length(t)-1);            % number of steps in RK4
ne=1000;                     % Monte-Carlo members
nef=1;                      % number of members for a reference orbit
delta=1.0;                  % error growth interval (lead time) 
Cd=1;                       % drag coefficient
Ts=1;                       % SST effects
%
% Construct first an reference orbit, along which perturbations 
% will be added. Note that tau value must be less than the time t
%
tau_array=[];error_array=[];
U_mean=[];V_mean=[];B_mean=[];
for k = 1:nef
    t(1)=0; 
    x(1)=-0.01;  % initial condition for u
    y(1)=0.01;  % initial condition for v
    z(1)=0.01;  % initial condition for b
    %
    % define the deterministic ODE
    %
    f1=@(t,x,y,z) p*y*y-(p+1)*z-Cd*x*y;  % MSD+ ode 
    g1=@(t,x,y,z) -x*y-Cd*y*y; 
    p1=@(t,x,y,z) z*x+s*x+Ts*y-r*z;
    
    f2=@(t,x,y,z) p*y*y-(p+1)*z+Cd*x*y;  % MSD- ode 
    g2=@(t,x,y,z) -x*y+Cd*y*y; 
    p2=@(t,x,y,z) z*x+s*x-Ts*y-r*z;
    %
    % integrating the ODE using RK4
    %
    for i=1:n
            y0=y(i);
            if (y0 >=0)
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
            else
                k1=f2(t(i),x(i),y(i),z(i));
                l1=g2(t(i),x(i),y(i),z(i));
                m1=p2(t(i),x(i),y(i),z(i));
                k2=f2(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));     
                l2=g2(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
                m2=p2(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
                k3=f2(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                l3=g2(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                m3=p2(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
                k4=f2(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
                l4=g2(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
                m4=p2(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));                    
            end
            x(i+1) = x(i) + h*(k1 + 2*k2 + 2*k3 + k4)/6; 
            y(i+1) = y(i) + h*(l1 + 2*l2 + 2*l3 + l4)/6;
            z(i+1) = z(i) + h*(m1 + 2*m2 + 2*m3 + m4)/6;
    end
    x1(:,k)=x;
    y1(:,k)=y;
    z1(:,k)=z;                 
end
%
% compute the reference trajectory either by ensemble mean (nef > 1) 
% or picking one specific member (nef = 1)
%
U_mean = mean(x1,2);
V_mean = mean(y1,2);
B_mean = mean(z1,2);
%
% Loop thru all instances at which random error are added
%
for tau=0.5:1.0:30             
    tau_array=[tau_array;tau];
    %
    % define some parameters for each ensemble of error growth
    %
    ta1=round(tau/h);     % convert to time step
    u_tau=U_mean(ta1);    % initial condition for each tau
    v_tau=V_mean(ta1);    % initial condition for each tau
    b_tau=B_mean(ta1);    % initial condition for each tau
    tb=tau+delta;         % reference time tau+\Delta for checking error 
    tb1=round(tb/h);      % convert to time step
    vb=V_mean(tb1);       
    up=u_tau+eps_IC;      % perturbed initial condition for u 
    vp=v_tau+eps_IC;      % perturbed initial condition for v
    bp=b_tau+eps_IC;      % perturbed initial condition for b
    %
    % plot for each tau 
    %
%     figure(1);
%     xlabel('\tau (\Delta = 1.4)'); 
%     ylabel('Error Growth'); 
%     set(gca, 'GridLineStyle', ':') %dotted grid lines
%     set(gca,'FontSize',14,'LineWidth',2.75)
%     plot(tau, v_tau, '+red', 'markersize', 18)
%     hold on;
%     plot(tb, vb, '+blue', 'markersize', 18)
%     hold on;
    %
    % running an ensemble intergration for each tau and plot
    %
    t2=tau:h:(tau+delta);
    for k = 1:ne
        n2=length(t2);
        t2(1)=tau; 
        x2(1)=up; 
        y2(1)=vp; 
        z2(1)=bp;     
        %Cdr=1.0+0.1*sqrt(dt)*randn;
        Cdr = min(max(Cd+(eps_Cd*randn),0.5),2.5);
        Tsr = min(max(Ts+(eps_Ts*randn),0.5),2.0);
        Sr = min(max(s+(eps_Ss*randn),0.0),1.0);
%         %
        % define the deterministic ODE
        %
        f1=@(t2,x2,y2,z2,Cdr,Tsr,Sr) p*y2*y2-(p+1)*z2-Cdr*x2*y2;  % MSD+ ode 
        g1=@(t2,x2,y2,z2,Cdr,Tsr,Sr) -x2*y2-Cdr*y2*y2; 
        p1=@(t2,x2,y2,z2,Cdr,Tsr,Sr) z2*x2+Sr*x2+Tsr*y2-r*z2;

        f2=@(t2,x2,y2,z2,Cdr,Tsr,Sr) p*y2*y2-(p+1)*z2+Cdr*x2*y2;  % MSD- ode 
        g2=@(t2,x2,y2,z2,Cdr,Tsr,Sr) -x2*y2+Cdr*y2*y2; 
        p2=@(t2,x2,y2,z2,Cdr,Tsr,Sr) z2*x2+Sr*x2-Tsr*y2-r*z2;
        %
        % integrating the ODE using RK4
        %
            for i=1:n2-1
                y0=y2(i);
%                 Cdr = min(max(Cd+(eps_Cd*randn),0.5),2.5);
%                 Tsr = min(max(Ts+(eps_Ts*randn),0.5),2.0);
%                 Sr = min(max(s+(eps_Ss*randn),0.0),1.0);
%         
                if (y0 >=0)
                    k21=f1(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    l21=g1(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    m21=p1(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    k22=f1(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);     
                    l22=g1(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);
                    m22=p1(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);
                    k23=f1(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr);
                    l23=g1(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr);
                    m23=p1(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr); 
                    k24=f1(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                    l24=g1(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                    m24=p1(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                else
                    k21=f2(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    l21=g2(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    m21=p2(t2(i),x2(i),y2(i),z2(i),Cdr,Tsr,Sr);
                    k22=f2(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);     
                    l22=g2(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);
                    m22=p2(t2(i)+h/2,(x2(i)+0.5*k21*h),(y2(i)+(0.5*l21*h)),(z2(i)+(0.5*m21*h)),Cdr,Tsr,Sr);
                    k23=f2(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr);
                    l23=g2(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr);
                    m23=p2(t2(i)+h/2,(x2(i)+0.5*k22*h),(y2(i)+(0.5*l22*h)),(z2(i)+(0.5*m22*h)),Cdr,Tsr,Sr);
                    k24=f2(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                    l24=g2(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                    m24=p2(t2(i)+h,(x2(i)+k23*h),(y2(i)+l23*h),(z2(i)+m23*h),Cdr,Tsr,Sr);
                end
                x2(i+1) = x2(i) + h*(k21 +2*k22  +2*k23   +k24)/6+eps_RF*sqrt(dt)*randn; %final equations with noise amplitude 0.01
                y2(i+1) = y2(i) + h*(l21  +2*l22   +2*l23    +l24)/6+eps_RF*sqrt(dt)*randn;
                z2(i+1) = z2(i) + h*(m21+2*m22 +2*m23  +m24)/6+eps_RF*sqrt(dt)*randn;                
            end
            x3(:,k)=x2;
            y3(:,k)=y2;
            z3(:,k)=z2;
            if (y2(i+1)>0)
                v(k)=y2(i+1);% stored v values at time tau+delta
            end
            %plot(t2(1:n2), y3(1:n2,k), '-g', 'markersize', 3)
    end
    % 
    % compute error growth rate here
    %
    sum=0.0;
    for k=1:ne
        sum1=((v(k)-vb)*(v(k)-vb))/(delta*delta);
        %sum1=((v(k)-vb)*(v(k)-vb));
        sum=sum+sum1;
    end
    error=sqrt(sum/ne);
    error_array=[error_array;error];
end
figure(2);
plot(tau_array, error_array, '-black', 'markersize', 10)
xlabel('\tau (\Delta = 1.4)'); 
ylabel('Error Growth'); 
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontSize',14,'LineWidth',2.75)
