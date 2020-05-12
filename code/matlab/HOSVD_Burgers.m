% This code creates the sample_soln.mat. Make sure to run this and save your
% environment before using *any* other main functionality examining the
% HOSVD.

Nx = 80;
Nt = 201;
Nq1 = 10;
Nq2 = 10;
Neps = 20;
epsilon = linspace(0.01,1,Neps);
q1 = linspace(0,1,Nq1);
q2 = linspace(0,1,Nq2);
addpath('./examples')

Y = zeros(2*Nx + 1, Nt, Neps, Nq1, Nq2);
for j=1:numel(epsilon)
    for k=1:numel(q1)
        for l=1:numel(q2)
            [z,x,t,e_conn] = burgers_1d_periodic(epsilon(j),q1(k),q2(l),Nx,Nt);
            Y(:,:,j,k,l) = z;
        end
    end
end