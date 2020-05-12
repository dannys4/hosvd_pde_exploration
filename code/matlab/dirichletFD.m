function A = dirichletFD(numx,delx)
    % via MathWorks spdiags page
    e = ones(numx,1)/(delx^2);
    A = spdiags([e -2*e e], -1:1, numx, numx);
    %A(:,1) = 0; A(:,numx) = 0;
end