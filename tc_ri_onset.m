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
%
% [REF]: Kieu, C. Q., and Q. Wang, 2017: JAS, doi/pdf/10.1175/JAS-D-17-0028.1
%
%==========================================================================
close all; clear
%
% Initialize parameters, using full Eqs. (61)-(63) in KW17.
%
p1 = 200;    % p parameter: ratio of PBL over depth of troposphere 
p2 = p1+1;   % aspect ratio R/H
p3 = 1.0;    % storm size scale
r = 0.25;    % radiative forcing per day
s = 0.1;     % s parameter: stratification
f2 = 0.00;   % Coriolis force
f1 = p1*f2;  % Coriolis force
n = 30000;   % number of integrations
dt = 0.001;  % time step
a = 0.01;    % std of u stochastic forcing -> v-variance is ~ 0.1%
b = 0.01;    % std of v stochastic forcing -> v-variance is ~ half
c = 0.01;    % std of b stochastic forcing -> v-variance is ~ half
ne = 1;    % number of Monte-Carlo integrations
plot_1D = 0;
plot_3D = 1;
%
% Set HSD initial conditions by creating 4 different initial points in the
% phase space of (u,v,b)
% 
u_0s = [-0.01, -1.0, -1.0, -0.1];
v_0s = [0.05,  1.0,  1.4, -0.2];
b_0s = [0.01,   0.5,  1.0,  0.1];
%b_0s(1) = v_0s(1)^2;    % constrain the initial condition for b_0
%u_0s(1) = -v_0s(1);     % constrain the initial condition for u_0
%
% integrate the MSD system in time for the first initial condition, 
% using RK4 order and plot right away
%
onset_stat=zeros(ne,1);
figure('Position',[0 0 700 550]);
for trial=1:1
for k = 1:ne
    [t1,u1,v1,b1] = tc_ri_onset_rk4(u_0s(trial), v_0s(trial), b_0s(trial), p1, p2, p3, r, s, f1, f2, dt, n,a,b,c);
%
% plot time series for each member
%
    if (plot_1D == 0)
        plot(t1(1:n),v1(1:n), 'LineWidth', 1.8)
        hold on
        xlabel('Nondimentional time'); 
        ylabel('Nondimentional maximum tangential wind v');
        title(['Tangential wind (v) time series'])
    else
%
% plot 3D for each member
%
        plot3(u1(1:n),v1(1:n),b1(1:n), '.', 'MarkerSize', 2)
        xlabel('normalized tangential wind u'); 
        ylabel('normalized vertical wind v'); 
        zlabel('normalized temperature anomaly b')
        grid on
        view(58,9)
        set(gca, 'GridLineStyle', ':') %dotted grid lines
        set(gca,'FontSize',14,'LineWidth',2.75)
        hold on
        plot3(u1(1),v1(1),b1(1), '+', 'Color','red',  'MarkerSize', 25)
        hold on
    end    
%
% searching for the onset time that first hit the level 0.1. By default,
% the hitting time is set to -999 to signal that the onset time does not
% exist.
%
    onset_stat(k) = NaN;
    for i = 1:n
       if (v1(i)>0.1)
           onset_stat(k) = t1(i);
           break
       end 
    end       
end
end
print('fig_V_time_series_MonteCarlo', '-dpng', '-r450');
%
% plot histogram
%
figure('Position',[0 0 500 300]);
edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ...
         18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, inf];
histogram(onset_stat,edges);
xlabel('RI onset nondimensional time'); 
ylabel('Frequency');
fprintf('Mean RI onset time hit level 0.1 is %11.3f \n',mean(onset_stat));
fprintf('Standard deviation RI onset time hit level 0.1 is %11.3f \n',std(onset_stat));
title(['RI onset hit time histogram'])
print('fig_RI_onset historgram', '-dpng', '-r450');
%
% print box plot
%
figure3=figure;
axes1 = axes('Parent',figure3);
boxplot(onset_stat,'BoxStyle','outline','Widths',0.2);
box(axes1,'on');
set(axes1,'XTickLabel',{'b=0.01'},'LineWidth',2);
set(gca,'FontSize',11,'LineWidth',2.)
xlabel('Experiments') 
ylabel('RI onset nondimensional time')
title(['RI boxplot distribution'])
print('fig_ri_onset_boxplot_dist', '-dpng', '-r450');
