function [b,A,B,a0] = burgers_rom(U,POD,x,e_conn,p,epsilon,w0,M)
%-------------------------------------------------------------------------------
%  burgers_rom - Builds a reduced-order model for the 1d viscous Burgers eqn.
%                using a given set of M-orthogonal basis functions.
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [b,A,B] = burgers_rom(U,POD,x,e_conn,p,epsilon,w0)
%
%  Variables:  U
%                        a steady-state centering solution
%              POD
%                        reduced-order basis for Galerkin projection of weak
%                        form of Burgers equation
%              x,e_conn
%                        finite element nodes and connectivity information
%              p
%                        order of the reduced-order model
%              epsilon
%                        viscosity parameter (default 0.01)
%              w0
%                        the uncentered initial solution (optional)
%                        used to compute the (optional) a0 initial
%                        conditions
%              M
%                        mass matrix, or discrete weighted inner product
%                        (optional, will be computed as mass matrix if not
%                        present)
%%
%  Example usage:
%              [b,A,B,a0] = test_model(U,POD,6,0.01,M,w0);
%
%       The reduced-order model (ROM) is computed by Galerkin projection of
%       Burgers equation onto the POD basis.  The result is a differential
%       equation in the following form:
%       
%       \dot{a} = b + A a + a'*B*a 
%
%          a(0) = Phi'*M*(w(0)-U)   (where w(0) is the vector of FEM 
%                                    coefficients corresponding to the 
%                                    initial condition for Burgers equation)
%% -----------------------------------------------------------------------------

  addpath('/Volumes/borggaard/Software/FEM/fem_functions')

  error('in progress')

  if ( nargin<5 )
    epsilon = 0.01;
    p       = 6;
  end

    % number of points used in Gaussian integration
  n_gauss               = 3;


  [n_nodes   , n_dimensions] = size(x     );
  [n_elements, nel_dof     ] = size(e_conn);


  %% ---------------------------------------------------------------------------
  %  Compute the ROM Matrices
  %-----------------------------------------------------------------------------
  ide(1:n_nodes-1) = 1:n_nodes-1;  ide(n_nodes) = 1;   % periodic boundary condition

  % Generate T(i).tensor matrices...
  II     = zeros(n_elements*nel_dof^2,1 );
  JJ     = zeros(n_elements*nel_dof^2,1 );
  KK     = zeros(n_elements*nel_dof^2,1 );
  TB.mat = zeros(n_elements*nel_dof^2,1 );

  TT     = repmat(TB, p, 1);

  [r,w] = oned_gauss(n_gauss);
  ep    = epsilon*ones(size(w));

  entry_counter = 0;
  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points
    nodes_local       = e_conn(n_el,:);
    x_local           = x(nodes_local,:);
    [ ~, w_g, h, h_x] = oned_shape(x_local,r,w);

    pod = h*POD(nodes_local,:);

    K_loc = oned_bilinear( ep, h_x, h_x, w_g );

    %-------------------------------------------------------------------
    %  Assemble contributions into the global system matrix
    %-------------------------------------------------------------------
    for n_t=1:nel_dof
      for n_u=1:nel_dof
        entry_counter = entry_counter + 1;
        II(entry_counter) = ide(nodes_local(n_t));
        JJ(entry_counter) = ide(nodes_local(n_u));
        KK(entry_counter) = K_loc(n_t,n_u);
        for k=1:p
          T_loc = oned_bilinear(pod(:,k),h,h_x,w_g);
          TT(k).mat(entry_counter) = T_loc(n_t,n_u);
        end
      end
    end
  end

  K = sparse( II(1:entry_counter), ...
              JJ(1:entry_counter), ...
              KK(1:entry_counter), n_nodes-1, n_nodes-1 );                  %#ok

  for k=1:p
    T(k).mat = sparse( II       (1:entry_counter), ...
                       JJ       (1:entry_counter), ...
                       TT(k).mat(1:entry_counter), n_nodes-1, n_nodes-1 );  %#ok
  end

  %  Assemble System Matrices
%  b   = POD(1:end-1,1:p)'*A_c*w_mean(1:end-1);
  A   = POD(1:end-1,1:p)'*K*POD(1:end-1,1:p);
  % for i=1:p
  %   b(i)   = b(i)   + w_mean(1:end-1)'*T(i).mat *w_mean(1:end-1);
  %   A(i,:) = A(i,:) + w_mean(1:end-1)'*T(i).mat *POD(1:end-1,1:p) ...
  %                   + w_mean(1:end-1)'*T(i).mat'*POD(1:end-1,1:p);
  % end

  for i=1:p
    B(i).mat = POD(1:end-1,1:p)'*T(i).mat*POD(1:end-1,1:p);                 %#ok
  end


   
  %  Integrate the amplitude equations (explicit Euler method)
  w0 = zeros(n_nodes,1);
  for n_nd=1:n_nodes
    if (x(n_nd)<=0.5)
      w0(n_nd) = q1*sin(2*pi*x(n_nd)) + q2*sin(4*pi*x(n_nd));
    else
      w0(n_nd) = 0;
    end
  end

  a0  = POD(:,1:p)'*M*w0;

  t = linspace(0,10,(Nt+1));

  n_sub = 200;
  delta_t = (t(2)-t(1))/n_sub;
  n_timesteps = (length(t)-1)*n_sub;
  w_pod_store = zeros(n_nodes,Nt+1);
  kount       = 1;
  w_pod_store(:,kount) = POD(:,1:p)*a0;
  a_c = a0;
  
  rhs = zeros(p,1);
  
  for n=1:n_timesteps
    for i=1:p
      rhs(i) = -A(i,:)*a_c - a_c'*B(i).mat*a_c;
    end
    a_n = a_c + delta_t*rhs;   % solution at n*delta_t
    if (~mod(n,n_sub))
      kount = kount + 1;
      w_pod_store(:,kount) = POD(:,1:p)*a_n;
      %plot(x,POD(:,1:p)*a_n)
      %pause(.001)
    end
    a_c = a_n;
  end

%   figure
%   mesh(x,t,w_pod_store')
%   view([.7 -.8 .6])

  
  difference = w_save-w_pod_store;
%   figure(3)
%   mesh(x,t,difference')


  ROM_error = 0.5*difference(:,1)'*M*difference(:,1);
  for i=2:Nt
    ROM_error = ROM_error + difference(:,i)'*M*difference(:,i);
  end
  ROM_error = ROM_error + 0.5*difference(:,Nt+1)'*M*difference(:,Nt+1);
  
  ROM_error = sqrt(delta_t*ROM_error);  % ( \int_0^T \| error(t,.) \|^2 dt )^{1/2}
end