function [x,t,p,SVD,Y] = HOSVD_heat_equation(conds)

    % Solve following equation
    % u_t(t,x;p) = f(p)*u_xx(t,x;p)
    % f(p) = 1, u(0,x;p) = u_0(x,p);
    % Solution: u(t,x;p) = p*exp(-pi^2*t)*cos(pi*x)
    % Create theta_0, f(p)

    theta0 = @(x,p) cos(pi*x)*(p');
    f = @(p) ones(numel(p),1);
    finv = @(u) ones(numel(u),1);


    % Must be careful with initial conditions
    minx = 0; maxx = 1; delx = 0.005;
    mint = 0; maxt = 1; delt = 0.005;

    % Finite difference dictates the maximum value of p
    minp = 0; maxp = finv((delx^2)/(4*delt)); delp = (maxp-minp)/1000;

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
              K(:,j) = (I-delt*fp(k)*A)\K(:,j-1);
        end
        Y(:,:,k) = K;
    end

    % Create tensor and SVD of it
    T = tensor(Y);
    SVD = hosvd(T,0.1);
end

