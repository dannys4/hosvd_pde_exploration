%% Danny Sharp Final Paper Associated Code
% Under Dr. Jeff Borggaard, Spring 2020
% Investigation into the structure of the HOSVD in PDE Solutions
addpath('./examples');
fig_folder = '../../finalPaper/figures/';
fig_num = 1;
%% Using the Heat Equation
% Solve the heat equation
[x,t,p,SVD,Y] = HOSVD_heat_equation('N');
% Set the analytical solutions at the points
X = cos(pi*x);
T = exp(-pi^2*t);
P = p;

% Normalize the vectors and change sign accordingly
core = double(SVD.core);
c = sign(core)*norm(X)*norm(T)*norm(P);
X = sign(SVD.U{1}(1))*X/norm(X);
T = sign(SVD.U{2}(1))*T/norm(T);
P = sign(SVD.U{2}(1))*P/norm(P);

% Analyze error
rel_err = abs(c-core)/abs(c);
fprintf('Relative difference in Analytical Core and Actual Core:\n')
fprintf('c=%f, core=%f, |c-core|/|c| = %f\n', c, core, rel_err)

%% Comparison plots
clf;

% Comparison in x
figure(fig_num); hold on;
plot(x,SVD.U{1},'r','LineWidth',2);
plot(x,X,'b--','LineWidth',2);
leg = legend('HOSVD vector in $x$','Analytical $X$');
leg.Location = 'northwest';
set(leg, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
title('Comparing Analytical value of $X$ to the HOSVD Approximation',...
    'Interpreter', 'latex');
hold off;
fig_heat_x = sprintf('%s%s', fig_folder, 'heat_eqn_comp_x.pdf');
printpdf(figure(fig_num), fig_heat_x); fig_num=fig_num+1;

% Comparison in t
figure(fig_num); hold on;
plot(t,SVD.U{2},'r','LineWidth',2);
plot(t,T,'b--','LineWidth',2);
leg = legend('HOSVD vector in $t$','Analytical $T$');
leg.Location = 'northwest';
set(leg, 'Interpreter', 'latex');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$f(t)$', 'Interpreter', 'latex');
title('Comparing Analytical value of $T$ to the HOSVD Approximation',...
    'Interpreter', 'latex');
hold off;
fig_heat_t = sprintf('%s%s', fig_folder, 'heat_eqn_comp_t.pdf');
printpdf(figure(fig_num), fig_heat_t); fig_num=fig_num+1;

% Comparison in p
figure(fig_num); hold on;
plot(p,SVD.U{3},'r','LineWidth',2);
plot(p,P,'b--','LineWidth',2);
leg = legend('HOSVD vector in $p$','Analytical $P$');
leg.Location = 'northwest';
set(leg, 'Interpreter', 'latex');
xlabel('$p$', 'Interpreter', 'latex');
ylabel('$f(p)$', 'Interpreter', 'latex');
title('Comparing Analytical value of $P$ to the HOSVD Approximation',...
    'Interpreter', 'latex');
hold off;
fig_heat_p = sprintf('%s%s', fig_folder, 'heat_eqn_comp_p.pdf');
printpdf(figure(fig_num), fig_heat_p); fig_num=fig_num+1;

%% Burger's Equation Section
load('sampleSoln.mat'); clf;

%% Examination in x
fig_num=1;
% u_{1_x} through u_{4_x}
figure(fig_num); hold on;
plot(x, SVD.U{1}(:,1:4));
leg = legend('$\mathbf{u}_{1_x}$','$\mathbf{u}_{2_x}$',...
    '$\mathbf{u}_{3_x}$', '$\mathbf{u}_{4_x}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_x}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_x$',...
    'Interpreter', 'latex');
fig_ux_1_4 = sprintf('%s%s', fig_folder, 'burgers_ux_1_4.pdf');
printpdf(figure(fig_num), fig_ux_1_4); fig_num=fig_num+1;
hold off;


% u_{7_x} through u_{10_x}
figure(fig_num); hold on;
plot(x, SVD.U{1}(:,7:10));
leg = legend('$\mathbf{u}_{7_x}$','$\mathbf{u}_{8_x}$',...
    '$\mathbf{u}_{9_x}$', '$\mathbf{u}_{10_x}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_x}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_x$',...
    'Interpreter', 'latex');
fig_ux_7_10 = sprintf('%s%s', fig_folder, 'burgers_ux_7_10.pdf');
printpdf(figure(fig_num), fig_ux_7_10); fig_num=fig_num+1;
hold off;

%% Examination in t
fig_num=1;
% u_{1_t} through u_{4_t}
figure(fig_num); hold on;
plot(t, SVD.U{2}(:,1:4));
leg = legend('$\mathbf{u}_{1_t}$','$\mathbf{u}_{2_t}$',...
    '$\mathbf{u}_{3_t}$', '$\mathbf{u}_{4_t}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_t}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_t$',...
    'Interpreter', 'latex');
fig_ut_1_4 = sprintf('%s%s', fig_folder, 'burgers_ut_1_4.pdf');
printpdf(figure(fig_num), fig_ut_1_4); fig_num=fig_num+1;
hold off;


% u_{7_t} through u_{10_x}
figure(fig_num); hold on;
plot(t, SVD.U{2}(:,7:10));
leg = legend('$\mathbf{u}_{7_t}$','$\mathbf{u}_{8_t}$',...
    '$\mathbf{u}_{9_t}$', '$\mathbf{u}_{10_t}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_t}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_t$',...
    'Interpreter', 'latex');
fig_ut_7_10 = sprintf('%s%s', fig_folder, 'burgers_ut_7_10.pdf');
printpdf(figure(fig_num), fig_ut_7_10); fig_num=fig_num+1;
hold off;

%% Examination in epsilon
fig_num=1;
% u_{1_eps} through u_{4_eps}
figure(fig_num); hold on;
plot(epsilon, SVD.U{3}(:,1:4));
leg = legend('$\mathbf{u}_{1_\varepsilon}$','$\mathbf{u}_{2_\varepsilon}$',...
    '$\mathbf{u}_{3_\varepsilon}$', '$\mathbf{u}_{4_\varepsilon}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$\varepsilon$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_\varepsilon}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_\varepsilon$',...
    'Interpreter', 'latex');
fig_ueps_1_4 = sprintf('%s%s', fig_folder, 'burgers_ueps_1_4.pdf');
printpdf(figure(fig_num), fig_ueps_1_4); fig_num=fig_num+1;
hold off;


% u_{5_eps} through u_{7_eps}
figure(fig_num); hold on;
plot(epsilon, SVD.U{3}(:,5:7));
leg = legend('$\mathbf{u}_{5_\varepsilon}$','$\mathbf{u}_{6_\varepsilon}$',...
    '$\mathbf{u}_{7_\varepsilon}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$\varepsilon$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_\varepsilon}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_\varepsilon$',...
    'Interpreter', 'latex');
fig_ueps_5_7 = sprintf('%s%s', fig_folder, 'burgers_ueps_5_7.pdf');
printpdf(figure(fig_num), fig_ueps_5_7); fig_num=fig_num+1;
hold off;

%% Examination in q1
fig_num=1; clf;
% u_{1_q1} through u_{3_q1}
figure(fig_num); hold on;
plot(q1, SVD.U{4}(:,1:3));
leg = legend('$\mathbf{u}_{1_{q_1}}$','$\mathbf{u}_{2_{q_1}}$',...
    '$\mathbf{u}_{3_{q_1}}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$q_1$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_{q_1}}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_{q_1}$',...
    'Interpreter', 'latex');
fig_uq1_1_3 = sprintf('%s%s', fig_folder, 'burgers_uq1_1_3.pdf');
printpdf(figure(fig_num), fig_uq1_1_3); fig_num=fig_num+1;
hold off;


% u_{4_q1} through u_{6_q1}
figure(fig_num); hold on;
plot(q1, SVD.U{4}(:,4:6));
leg = legend('$\mathbf{u}_{4_{q_1}}$','$\mathbf{u}_{5_{q_1}}$',...
    '$\mathbf{u}_{6_{q_1}}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$q_1$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_{q_1}}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_{q_1}$',...
    'Interpreter', 'latex');
fig_uq1_4_6 = sprintf('%s%s', fig_folder, 'burgers_uq1_4_6.pdf');
printpdf(figure(fig_num), fig_uq1_4_6); fig_num=fig_num+1;
hold off;
%% Examination in q2
fig_num=1; clf;
% u_{1_q2} through u_{4_q2}
figure(fig_num); hold on;
plot(q2, SVD.U{5}(:,1:4));
leg = legend('$\mathbf{u}_{1_{q_2}}$','$\mathbf{u}_{2_{q_2}}$',...
    '$\mathbf{u}_{3_{q_2}}$', '$\mathbf{u}_{4_{q_2}}$');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$q_2$', 'Interpreter', 'latex');
ylabel('$\mathbf{u}_{j_{q_2}}$', 'Interpreter', 'latex');
title('Examining Singular Vectors in $\mathbf{U}_{q_2}$',...
    'Interpreter', 'latex');
fig_uq2_1_4 = sprintf('%s%s', fig_folder, 'burgers_uq2_1_4.pdf');
printpdf(figure(fig_num), fig_uq2_1_4); fig_num = fig_num+1;
hold off;

%% PREDICTION PERFORMANCE
trials = 500;
[interp_err, prop_err, rand_vals] = HOSVD_error_estimate(SVD,...
    epsilon, q1, q2, trials);
interp_err = interp_err*100;
prop_err = prop_err*100;
mean_interp = mean(interp_err);
mean_prop = mean(prop_err);
med_interp = median(interp_err);
med_prop = median(prop_err);
std_interp = std(interp_err);
std_prop = std(prop_err);
%% Approximation Methods
% Interpolation Method
fig_num=1; clf;
figure(fig_num); hold on;
histogram(interp_err,50); 
title('Histogram of Error in Interpolation Method', 'Interpreter', 'latex');
xlabel('Error of Interpolative Model', 'Interpreter', 'latex');
ylabel('Trials', 'Interpreter', 'latex');
ylim([0,325]);
fig_interp_error = sprintf('%s%s', fig_folder, 'interp_error.pdf');
printpdf(figure(fig_num), fig_interp_error);
hold off; fig_num = fig_num+1;

% Proportional Method
figure(fig_num); hold on;
histogram(prop_err, 50); 
title('Histogram of Error in Proportional Method', 'Interpreter', 'latex');
xlabel('Error of Proportional Model (Percent)', 'Interpreter', 'latex');
ylabel('Trials', 'Interpreter', 'latex');
ylim([0,325]);
fig_prop_error = sprintf('%s%s', fig_folder, 'prop_error.pdf');
printpdf(figure(fig_num), fig_prop_error); fig_num = fig_num+1;
hold off;
%% Plotting parameters compared to error
fig_num=1; clf;
app_err = [interp_err,prop_err];
% Epsilon vs. approximation error
figure(fig_num); hold on;
plot(rand_vals(:,1), app_err, 'o');
leg = legend('Interpolation', 'Proportional');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$\varepsilon$', 'Interpreter', 'latex');
ylabel('Error($\%$)', 'Interpreter', 'latex');
title('Error vs. $\varepsilon$ Value by Approximation Method',...
    'Interpreter', 'latex');
fig_eps_err = sprintf('%s%s', fig_folder, 'eps_err.pdf');
printpdf(figure(fig_num), fig_eps_err); fig_num = fig_num+1;
hold off;

% q1 vs. approximation error
figure(fig_num); hold on;
plot(rand_vals(:,2), app_err, 'o');
leg = legend('Interpolation', 'Proportional');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$q_1$', 'Interpreter', 'latex');
ylabel('Error($\%$)', 'Interpreter', 'latex');
title('Error vs. $q_1$ Value by Approximation Method',...
    'Interpreter', 'latex');
fig_q1_err = sprintf('%s%s', fig_folder, 'q1_err.pdf');
printpdf(figure(fig_num), fig_q1_err); fig_num = fig_num+1;
hold off;

% q2 vs. approximation error
figure(fig_num); hold on;
plot(rand_vals(:,3), app_err, 'o');
leg = legend('Interpolation', 'Proportional');
leg.Location = 'northeast';
set(leg, 'Interpreter', 'latex');
xlabel('$q_2$', 'Interpreter', 'latex');
ylabel('Error($\%$)', 'Interpreter', 'latex');
title('Error vs. $q_2$ Value by Approximation Method',...
    'Interpreter', 'latex');
fig_q2_err = sprintf('%s%s', fig_folder, 'q2_err.pdf');
printpdf(figure(fig_num), fig_q2_err); fig_num = fig_num+1;
hold off;
%% Sample Solution Burgers
% Take Sample parameter values
hateps = 0.07;
hatq1 = 0.42;
hatq2 = 0.23;

% Values to have on-hand
sizeSVD = size(SVD);
Nx = floor(sizeSVD(1)/2);
Nt = sizeSVD(2);

% Take approximate and true solutions
[z,~,~,~] = burgers_1d_periodic(hateps, hatq1, hatq2, Nx, Nt);
interp_z = interpSVDBurgers(SVD, hateps, hatq1, hatq2, epsilon, q1, q2);
prop_z = propSVDBurgers(SVD, hateps, hatq1, hatq2, epsilon, q1, q2);

% Evaluate relative error
interp_err = norm(z-interp_z)/norm(z);
prop_err = norm(z-prop_z)/norm(z);

fprintf('Error in interpolated approximation:%f\n', interp_err);
fprintf('Error in proportional approximation:%f\n', prop_err);
%% Plotting True and Approximate Solutions
fig_num = 1;

% Plot "true" solution
figure(fig_num); hold on;
mesh(x,t,z');
view([.7 -.8 .6]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$\hat{w}$', 'Interpreter', 'latex');
title('Isometric View of $\hat{w}(x,t;\hat{\varepsilon},\hat{q}_1,\hat{q}_2)$',...
    'Interpreter', 'latex');
fig_true_soln = sprintf('%s%s', fig_folder, 'true_soln.pdf');
printpdf(figure(fig_num), fig_true_soln); fig_num = fig_num+1;
hold off;

% Plot interpolated approximation
figure(fig_num); hold on;
mesh(x,t,interp_soln');
view([.7 -.8 .6]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$\mathcal{T}_i$', 'Interpreter', 'latex');
title('Isometric View of Interpolated Approximation to $\hat{w}$',...
    'Interpreter', 'latex');
fig_interp_soln = sprintf('%s%s', fig_folder, 'interp_soln.pdf');
printpdf(figure(fig_num), fig_interp_soln); fig_num = fig_num+1;
hold off;

% Plot proportional approximation
figure(fig_num); hold on;
mesh(x,t,prop_soln');
view([.7 -.8 .6]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex');
zlabel('$\mathcal{T}_p$', 'Interpreter', 'latex');
title('Isometric View of Proportional Approximation to $\hat{w}$',...
    'Interpreter', 'latex');
fig_prop_soln = sprintf('%s%s', fig_folder, 'prop_soln.pdf');
printpdf(figure(fig_num), fig_prop_soln); fig_num = fig_num+1;
hold off;