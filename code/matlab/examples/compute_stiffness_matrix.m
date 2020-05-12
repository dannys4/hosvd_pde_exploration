function [K] = compute_stiffness_matrix(x,e_conn)
%%----------------------------------------------------------------------
%  compute_stiffness_matrix - computes a stiffness matrix 
%                        \int(\nabla h_i(x) \cdot \nabla h_j(x)) dx
%                        from given nodal coordinates and element 
%                        connectivity
%                        
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [K] = compute_stiffness_matrix( x , e_conn );
%
%  Variables:  x      
%                      nodal coordinates
%              e_conn
%                      element connectivity
%
%              K
%                      (sparse) stiffness matrix
%% ---------------------------------------------------------------------


  [N,          dim    ] = size(x     );
  [n_elements, nel_dof] = size(e_conn);

  % Due to memory limitations, we will fill the mass matrix by integrating
  % over partitions of elements.  The constant (max_elem_per_partition)
  % should be adjusted based on available memory.
  max_elem_per_partition = 10000;
  
  n_part = floor( n_elements / (max_elem_per_partition+1) ) + 1;
  elem_segment = floor( linspace(0,n_elements,n_part+1) );
  max_part_size = max( diff( elem_segment ) );

  if (dim==1)
    [r,w]   = oned_quadrature(3);
    
    for n_pt = 1:n_part
      II = zeros( max_part_size*nel_dof^2,1 );  
      JJ = zeros( max_part_size*nel_dof^2,1 );
      XX = zeros( max_part_size*nel_dof^2,1 );
      entry_counter = 0;

      for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
        nodes_local      = e_conn(n_el,:);
        x_local          = x(nodes_local,:);
        [~, w_g, ~, p_x] = oned_shape(x_local,r,w);

        one = ones(size(w_g));
        K_loc = oned_bilinear( one, p_x, p_x, w_g );

        %  Assemble into the global system matrix
        for i=1:nel_dof
          for j=1:nel_dof
            entry_counter = entry_counter + 1;
            II(entry_counter) = nodes_local(i);
            JJ(entry_counter) = nodes_local(j);
            XX(entry_counter) = K_loc(i,j);
          end
        end
      end
    
      if ( n_pt==1 )
        K = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      else
        K = K + ...
            sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      end
      
    end % n_pt loop
    
    clear II JJ XX
  end % dim==1
  
  if (dim==2)
    [r,s,w] = twod_quadrature(7);
    one = ones(size(w));
    
    for n_pt = 1:n_part
      II = zeros( max_part_size*nel_dof^2,1 );
      JJ = zeros( max_part_size*nel_dof^2,1 );
      XX = zeros( max_part_size*nel_dof^2,1 );
      entry_counter = 0;

      for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
        nodes_local               = e_conn(n_el,:);
        x_local                   = x(nodes_local,:);
        [x_g, w_g, ~, p_x, p_y] = twod_shape(x_local,r,s,w);

        K_loc = twod_bilinear( one, p_x, p_x, w_g ) + ...
                twod_bilinear( one, p_y, p_y, w_g );

        %  Assemble into the global system matrix
        for i=1:nel_dof
          for j=1:nel_dof
            entry_counter = entry_counter + 1;
            II(entry_counter) = nodes_local(i);
            JJ(entry_counter) = nodes_local(j);
            XX(entry_counter) = K_loc(i,j);
          end
        end
      end
    
      if ( n_pt==1 )
        K = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      else
        K = K + ...
            sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      end
      
    end % n_pt loop
    
    clear II JJ XX
    
  end % dim==2

  if (dim==3)
    [r,s,t,w] = threed_quadrature(5);

    for n_pt = 1:n_part
      II = zeros( max_part_size*nel_dof^2,1 );  
      JJ = zeros( max_part_size*nel_dof^2,1 );
      XX = zeros( max_part_size*nel_dof^2,1 );
      entry_counter = 0;

      for n_el=elem_segment(n_pt)+1:elem_segment(n_pt+1)
      nodes_local                    = e_conn(n_el,:);
      x_local                        = x(nodes_local,:);
      [x_g, w_g, ~, p_x, p_y, p_z] = threed_shape(x_local,r,s,t,w);

      one = ones(size(w_g));
      K_loc = threed_bilinear( one, p_x, p_x, w_g ) + ...
              threed_bilinear( one, p_y, p_y, w_g ) + ...
              threed_bilinear( one, p_z, p_z, w_g );

        %  Assemble into the global system matrix
        for i=1:nel_dof
          for j=1:nel_dof
            entry_counter = entry_counter + 1;
            II(entry_counter) = nodes_local(i);
            JJ(entry_counter) = nodes_local(j);
            XX(entry_counter) = K_loc(i,j);
          end
        end
      end
    
      if ( n_pt==1 )
        K = sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      else
        K = K + ...
            sparse( II(1:entry_counter), JJ(1:entry_counter),...
                    XX(1:entry_counter),...
                    N, N );
      end
      
    end % n_pt loop
    
    clear II JJ XX
    
  end % dim==3
  
  
end % function compute_stiffness_matrix
