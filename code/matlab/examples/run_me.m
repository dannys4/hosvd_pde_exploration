%  A script to compute a POD model for Burgers equation by
%
%   (1)  computing simulation data:  burgers_1d_periodic
%   (2)  using the simulation to find a good low-dimensional basis
%        generate_pod
%   (3)  running test_model that itself builds the low-dimensional dynamical
%        system, simu lates it, reconstructs an approximation to the solution
%        to Burgers equation, then computes the error in the approximation
%        (difference between FEM solution and POD solution).

%  addpath('/Volumes/borggaard/Software/ROM_Software/POD')

% Compute basis at p1
%-----------------------------------------------------------------------------
R  = 1000;
q1 = 0.25;
q2 = 0.37;
r_dim  = 12;   % the number of POD basis functions to use...
Nx = 80;
Nt = 201;

epsilon = 0.01;
p1 = [ epsilon; q1; q2 ];

% the 80 and 201 below are spatial and temporal FEM discretization parameters
[z,x,t,e_conn] = burgers_1d_periodic(epsilon,q1,q2,Nx,Nt);

%  Comment out the eye candy below for parametric experiments...
figure(1)
mesh(x,t,z')
view([.7 -.8 .6])
pause(0.001)
xlabel('x'); ylabel('t'); zlabel('z')
title('Finite Element Solution')

[M] = compute_mass_matrix(x,e_conn);

[POD1,lam1] = generate_pod(z,r_dim,M);


[ROM_error] = test_model(POD1,r_dim,epsilon,q1,q2,M,z);
w_norm = 0.5*z(:,1)'*M*z(:,1);
for i=2:Nt-1
    w_norm = w_norm + z(:,i)'*M*z(:,i);
end
w_norm = w_norm + 0.5*z(:,Nt)'*M*z(:,Nt);

w_norm = sqrt(10*w_norm/Nt);  % ( \int_0^T \| error(t,.) \|^2 dt )^{1/2}
fprintf('At P_1, p_1 = %g, p_2 = %g, and p_3 = %g,\n',epsilon,q1,q2);
fprintf('  the ROM_error is %g\n',ROM_error);
fprintf('  the relative error is %g\n', w_norm);


