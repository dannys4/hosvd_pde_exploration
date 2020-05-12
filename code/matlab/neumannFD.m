function A = neumannFD(numx,delx)
    A = -2*eye(numx) +...
    diag(diag(eye(numx-1)),-1) +...
    diag(diag(eye(numx-1)),1);

    A(1,2) = 2; A(numx, numx-1) = 2; A = A/(delx^2);
end