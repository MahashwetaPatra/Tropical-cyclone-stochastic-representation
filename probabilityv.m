%
% [NOTE]: This program solves for the stochastic MSD system in (Kieu and Wang
%         2017), using RK4 for deterministic forcings, and Euler method 
%         for the stochastic forcing components. 
%
% [AUTHOR]: Chanh Kieu, Indiana University
%
% [HIST]: - May 26, 2019: created from MSD system by CK
%         - May 14, 2020: re-designed for adding Monte-Carlo integration
%                         directly without using Bash shell loop
%         - May 17, 2020: added few statistical analyses and historgram
%                         and boxplot figure options.
%         - May 22, 2020: corrected time label (nondimensional time)
%         - Jun 15, 2020: add plotting 3D option and constraints for ini
%         - February 26, 2021: For different v_0 and \epsilon, calculates
%           the fractions of realizations that reaches l before it reaches 0,
%           by Patra
%         -May 23, 2021 Error Bars with 95% confidence interval is plotted
%           by Patra 
% [REF]: Kieu, C. Q., and Q. Wang, 2017: JAS, doi/pdf/10.1175/JAS-D-17-0028.1
%
%==========================================================================
clc; close all; clear all;
% Initialize parameters, using full Eqs. (61)-(63) in KW17.
tic
p1 = 200;    % p parameter: ratio of PBL over depth of troposphere 
p2 = p1+1;   % aspect ratio R/H
p3 = 1.0;    % storm size scale
r = 0.25;    % radiative forcing per day
s = 0.1;     % s parameter: stratification
f2 = 0.00;   % Coriolis force
f1 = p1*f2;  % Coriolis force
n = 30000;   % number of integrations
dt = 0.01;  % time step
a = 0.0001;    % std of u stochastic forcing -> v-variance is ~ 0.1%
b = 0.0001;    % std of v stochastic forcing -> v-variance is ~ half
c = 0.0001;    % std of b stochastic forcing -> v-variance is ~ half
ne = 100;    % number of Monte-Carlo integrations
v_array=[];STD=[];P_array=[];
for trial_0s=0.001:0.001:0.01,% trial is the initial value of v_0
    t(1)=0;  %initializing x,y,z,t
    %v_array=[v_array;trial];
    v_array=[v_array;trial_0s];
    %
    % Set HSD initial conditions by creating 4 different initial points in the
    % phase space of (u,v,b)
    %
    u_0s = [-0.01, -1.0, -1.0, -0.1];
    %trial_0s = [0.003,  0.005, 0.01, 0.05];
    b_0s = [0.0001,   0.5,  1.0,  0.1];
    %v_0s(1)=trial_0s(trial);
    v_0s(1)=trial_0s;
    Fraction=[];param=[];
    for N=1:10,
        ne = 100;    % number of Monte-Carlo integrations
        onset_stat=zeros(ne,1); %Preallocating the dimensions of onset-vectors
        onset_stat_plus=zeros(ne,1);
        onset_stat_zero=zeros(ne,1);
%
% searching and sorting which and when trajectories hit M^+ or M^0.
% For tauplus we take into account only the trajectories that remain
% positive all the time.
%   
        for k = 1:ne
            [t1,u1,v1,b1] = tc_ri_onset_rk4(u_0s(1), v_0s(1), b_0s(1), p1, p2, p3, r, s, f1, f2, dt, n,a,b,c);
            onset_stat(k) = NaN;
            onset_stat_plus(k) = NaN;
            onset_stat_zero(k) = NaN;
            for i = 1:n
                if (v1(i)<0)
                    onset_stat_plus(k) = 0;
                    onset_stat_zero(k) = t1(i);
                    onset_stat(k)=t1(i);
                    %Here onset_stat just identifies the time it hit M^+ or M^0
                    break
                elseif (v1(i)>0.1)
                    onset_stat_plus(k) = t1(i);
                    onset_stat_zero(k) = 0;
                    onset_stat(k)=t1(i);
                    break 
                elseif (0<v1(i)<0.1)
                    onset_stat_plus(k) = 0;
                    onset_stat_zero(k) = 0;
                    onset_stat(k)=0;
                end
            end       
        end
        onset_stat=onset_stat_plus+onset_stat_zero; %Array of hitting times for the manifold M+0
        Hitting_times_Mplus = onset_stat_plus(onset_stat_plus~=0); %modifying the arrays cropping zeros
        Hitting_times_Mzero = onset_stat_zero(onset_stat_zero~=0);
        %Probabilities to hit M^+ or M^-
        [dplus , numcolplus] = size(Hitting_times_Mplus);
        [dzero , numcolzero] = size(Hitting_times_Mzero);
        ProbHitMplus = dplus/ne
        ProbHitMzero = dzero/ne ;
        Fraction=[Fraction;ProbHitMplus];
        mplus=mean (Hitting_times_Mplus);
        sigmaplus= std(Hitting_times_Mplus);
        Varplus= sigmaplus^2;
        STD=[STD;Varplus]; %% keeping the variance values in an aaray
    end
    % v_array
    % Fraction  %%array of probalities to reach l before 0 for different v_0
    % figure(4);
    % plot(v_array, Fraction, '.-b', 'markersize', 10)
    % xlabel('T (150<p<250)'); 
    % ylabel('H(T)'); 
    % set(gca, 'GridLineStyle', ':') %dotted grid lines
    % set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
    P=mean(Fraction);
    SF=std(Fraction);
    ZF=1.96;
    nF=10;
    P1=P+(ZF*SF)/sqrt(10);
    P2=P-(ZF*SF)/sqrt(10);
    PB=[trial_0s;trial_0s];
    EB=[P1;P2];
    P_array=[P_array;P];
    figure(4);
    set(gca, 'GridLineStyle', ':') %dotted grid lines
    set(gca,'FontName','Times','FontSize',24,'LineWidth',2.75)
    plot(PB,EB, '.-black', 'markersize', 20)
    hold on;
%    plot(trial_0s,P, '.-red', 'markersize', 10)
    xlabel('T (150<p<250)'); 
    ylabel('H(T)'); 
end
    plot(v_array,P_array, '.-red', 'markersize', 20)
 
