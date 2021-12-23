function [T,U,V,B] = tc_ri_onset_rk4(u_0, v_0, b_0, p1, p2, p3, r, s, f1, f2, dt, m,av,bv,cv)
%---------------------------------------------------------------------------
%MPI_RK4  Runge-Kotta order 4 solutions 
%         for all three of the three HSD equations
%            u_dot
%            v_dot
%            b_dot
%
% Sample call
%   [T,U,V,B] = rk4(...)
%     the capital letters show that they are vectors of many values 
%     and do not imply that they are not nondimensional
%
% Inputs
%   *_0     initial values for u, v, and b
%   s       squared ratio of height of tropopause to depth of boundary layer
%   k       nondimensional forcing
%   gamma   
%   dt      width of a step
%   m       number of steps
%   a       sqrt of the variance of u-forcing amplitude
%   b       sqrt of the variance of v-forcing amplitude
%   c       sqrt of the variance of b-forcing amplitude
%
% Return
%   T       solution: col vector of times (abscissas)
%   U,V,B   solutions (ordinates)
%---------------------------------------------------------------------------

%our set of differential equations for MSD+
up_dot = @(u,v,b) ( p1*v.^2 - p2*b - p3*u.*v + f1*v);
vp_dot = @(u,v,b) ( -u.*v - v.^2 - f2*u );
bp_dot = @(u,v,b) ( u.*b + s*u + v - r*b );
%our set of differential equations for MSD-
um_dot = @(u,v,b) ( p1*v.^2 - p2*b + p3*u.*v + f1*v);
vm_dot = @(u,v,b) ( -u.*v + v.^2 - f2*u );
bm_dot = @(u,v,b) ( u.*b + s*u - v - r*b );


%preallocate space to hold the times and solutions
T = zeros(m+1, 1);
U = zeros(m+1, 1);
V = zeros(m+1, 1);
B = zeros(m+1, 1);

%start at t=0 with the initial values
T(1) = 0;
U(1) = u_0;
V(1) = v_0;
B(1) = b_0;

%use the RK4 method in a loop to form solutions to the diff eqns
for j=1:m
    
    %step 0
    u0 = U(j);
    v0 = V(j);
    b0 = B(j);
    
    if (v0 >=0)
        %step 1
        u = u0;
        v = v0;
        b = b0;
        u1 = up_dot(u,v,b) * dt;
        v1 = vp_dot(u,v,b) * dt;
        b1 = bp_dot(u,v,b) * dt;
    
        %step 2
        u = u0 + u1 / 2;
        v = v0 + v1 / 2;
        b = b0 + b1 / 2;
        u2 = up_dot(u,v,b) * dt;
        v2 = vp_dot(u,v,b) * dt;
        b2 = bp_dot(u,v,b) * dt;
    
        %step 3
        u = u0 + u2 / 2;
        v = v0 + v2 / 2;
        b = b0 + b2 / 2;
        u3 = up_dot(u,v,b) * dt;
        v3 = vp_dot(u,v,b) * dt;
        b3 = bp_dot(u,v,b) * dt;
    
        %step 4
        u = u0 + u3;
        v = v0 + v3;
        b = b0 + b3;
        u4 = up_dot(u,v,b) * dt;
        v4 = vp_dot(u,v,b) * dt;
        b4 = bp_dot(u,v,b) * dt;
    else
        %step 1
        u = u0;
        v = v0;
        b = b0;
        u1 = um_dot(u,v,b) * dt;
        v1 = vm_dot(u,v,b) * dt;
        b1 = bm_dot(u,v,b) * dt;
    
        %step 2
        u = u0 + u1 / 2;
        v = v0 + v1 / 2;
        b = b0 + b1 / 2;
        u2 = um_dot(u,v,b) * dt;
        v2 = vm_dot(u,v,b) * dt;
        b2 = bm_dot(u,v,b) * dt;
    
        %step 3
        u = u0 + u2 / 2;
        v = v0 + v2 / 2;
        b = b0 + b2 / 2;
        u3 = um_dot(u,v,b) * dt;
        v3 = vm_dot(u,v,b) * dt;
        b3 = bm_dot(u,v,b) * dt;
    
        %step 4
        u = u0 + u3;
        v = v0 + v3;
        b = b0 + b3;
        u4 = um_dot(u,v,b) * dt;
        v4 = vm_dot(u,v,b) * dt;
        b4 = bm_dot(u,v,b) * dt;         
    end
    
    %combine the steps to give final answers
    T(j+1) = T(1) + dt*j;
    U(j+1) = u0 + (u1 + 2*u2 + 2*u3 + u4)/6 + av*sqrt(dt)*randn;
    V(j+1) = v0 + (v1 + 2*v2 + 2*v3 + v4)/6 + bv*sqrt(dt)*randn;
    B(j+1) = b0 + (b1 + 2*b2 + 2*b3 + b4)/6 + cv*sqrt(dt)*randn;
%NoNoise 
    %U(j+1) = u0 + (u1 + 2*u2 + 2*u3 + u4)/6 ;
    %V(j+1) = v0 + (v1 + 2*v2 + 2*v3 + v4)/6 ;
    %B(j+1) = b0 + (b1 + 2*b2 + 2*b3 + b4)/6 ; 
end

 %U(j+1) = u0 + (u1 + 2*u2 + 2*u3 + u4)/6 ;
    %V(j+1) = v0 + (v1 + 2*v2 + 2*v3 + v4)/6 ;
    %B(j+1) = b0 + (b1 + 2*b2 + 2*b3 + b4)/6 ; 

%U(j+1) = u0 + (u1 + 2*u2 + 2*u3 + u4)/6 + av*sqrt(dt)*randn;
   % V(j+1) = v0 + (v1 + 2*v2 + 2*v3 + v4)/6 + bv*sqrt(dt)*randn;
    %B(j+1) = b0 + (b1 + 2*b2 + 2*b3 + b4)/6 + cv*sqrt(dt)*randn; 

