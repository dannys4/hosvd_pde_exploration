function [interp_err, prop_err, rand_vals] = HOSVD_error_estimate(SVD,...
    epsvec, q1vec, q2vec, trials)
addpath('./examples');
% Come up with an average error for each of the solution estimate methods
sizeSVD = size(SVD);
Nx = floor(sizeSVD(1)/2);
Nt = sizeSVD(2);

% Create random values for epsilon, q1, and q2
eps_vals = 0.99*rand(trials,1) + 0.01;
q1_vals  = rand(trials,1);
q2_vals  = rand(trials,1);
rand_vals = [eps_vals, q1_vals, q2_vals];

% Initialize the error for each of the methods
interp_err = zeros(trials,1);
prop_err = zeros(trials,1);

% Try a number of trials to get a gauge on the error
for trial=1:trials
    epsilon = eps_vals(trial);
    q1 = q1_vals(trial);
    q2 = q2_vals(trial);
    
    % Get "true" solution
    [z,~,~,~] = burgers_1d_periodic(epsilon, q1, q2, Nx, Nt);
    
    % Get approximate solutions
    interp_soln = interpSVDBurgers(SVD,epsilon,q1,q2,epsvec,q1vec,q2vec);
    prop_soln = propSVDBurgers(SVD,epsilon,q1,q2,epsvec,q1vec,q2vec);
    
    % Calculate error
    interp_err(trial) = norm(z-interp_soln)/norm(z);
    prop_err(trial) = norm(z-prop_soln)/norm(z);
end

end