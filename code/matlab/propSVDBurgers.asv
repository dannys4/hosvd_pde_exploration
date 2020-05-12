function [soln] = propSVDBurgers(SVD,epsilon,q1,q2,epsvec,q1vec,q2vec)
    %% Estimate weights of different vectors via proximity
    % Then, find weights for the U matrices
    weps = 1./abs(epsilon-epsvec); weps=weps/norm(weps,1);
    wq1 = 1./abs(q1-q1vec); wq1=wq1/norm(wq1,1);
    wq2 = 1./abs(q2-q2vec); wq2=wq2/norm(wq2,1);
    
    % Create a "flattened" U, so as to make a 2-D ttensor
    Ueps = weps*SVD.U{3};
    Uq1 = wq1*SVD.U{4};
    Uq2 = wq2*SVD.U{5};
    U = {SVD.U{1}, SVD.U{2}, Ueps, Uq1, Uq2};
    
    % Create a solution with this interpolation.
    soln = double(ttensor(SVD.core, U));
end