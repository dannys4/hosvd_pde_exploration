function [soln] = interpSVDBurgers(SVD,epsilon,q1,q2,epsvec,q1vec,q2vec)
    %% Interpolate the SVD using convex combinations
    % First, find closest values in tested values
    [~,epsind] = min(abs(epsvec-epsilon));
    [~,q1ind] = min(abs(q1vec-q1));
    [~,q2ind] = min(abs(q2vec-q2));
    
    % Take care of case where minimum distance is last number
    if(epsvec(epsind) > epsilon)
        epsind = epsind - 1;
    end
    if(q1vec(q1ind) > q1)
        q1ind = q1ind - 1;
    end
    if(q2vec(q2ind) > q2)
        q2ind = q2ind - 1;
    end
    
    if epsind < 1 || q1ind < 1 || q2ind < 1
        fprintf('problem found\n');
    end
    if epsind > numel(epsvec) || q1ind > numel(q1vec) || q2ind > numel(q2vec)
        fprintf('problem found\n');
    end
    
    % Then, find weights for the U matrices
    weps = [epsilon-epsvec(epsind), epsvec(epsind+1)-epsilon];
    wq1 = [q1-q1vec(q1ind), q1vec(q1ind+1)-q1];
    wq2 = [q2-q1vec(q2ind), q2vec(q2ind+1)-q2];
    weps = weps/norm(weps,1); wq1 = wq1/norm(wq1,1); wq2 = wq2/norm(wq2,1);
    
    % Recreate U with interpolation, so as to make a 2-D ttensor
    Ueps = weps*SVD.U{3}(epsind:(epsind+1),:);
    Uq1 = wq1*SVD.U{4}(q1ind:(q1ind+1),:);
    Uq2 = wq2*SVD.U{5}(q2ind:(q2ind+1),:);
    U = {SVD.U{1}, SVD.U{2}, Ueps, Uq1, Uq2};
    
    % Create a solution with this interpolation.
    soln = double(ttensor(SVD.core, U));
end