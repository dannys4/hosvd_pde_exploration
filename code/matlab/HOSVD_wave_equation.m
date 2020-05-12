% Solve the following equation:
% u_tt = f(p) u_xx
% u(0,x;p) = v(x,p), u_t(0,x;p) = w(x,p)
% With associated boundary conditions
conds = input('What boundary conditions? Options: D(irichlet), N(eumann), R(obin):\n');
v0 = @(x,p) sin(pi*x)*(p');
w0 = @(x,p) cos(pi*x)*(p');
f = @(p) p;
finv = @(u) 1/u;


% Must be careful with initial conditions
minx = 0; maxx = 1; delx = 0.005;
mint = 0; maxt = 1; delt = 0.005;

% Finite difference dictates the maximum value of p
minp = 0; maxp = finv((delx^2)/(4*delt^2)); delp = (maxp-minp)/1000;

% Set x,t,p vectors
x = (minx:delx:maxx)';
t = (mint:delt:maxt)';
p = (minp:delp:maxp)';
fp = f(p);

% Create spatial finite difference matrix and initial condition
if conds == 'D'
    A = dirichletFD(numel(x),delx);
elseif conds == 'N'
    A = neumannFD(numel(x),delx);
else
    A = robinFD(numel(x),delx);
end

% Initialize multi-dim array for tensor
Y = zeros(numel(x), numel(t), numel(fp));
Y(:,1,:) = theta0(x,p);
I = speye(numel(x));

for k=1:numel(fp)
    K = zeros(numel(x),numel(t));
    K(:,1) = theta0(x,p(k));
    for j=2:numel(t)
          K(:,j) = (I-delt^2*fp(k)*A)\K(:,j-1);
    end
    Y(:,:,k) = K;
end

% Create tensor and SVD of it
T = tensor(Y);
SVD = hosvd(T,0.1);
