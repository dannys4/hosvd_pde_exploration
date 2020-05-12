function A = robinFD(numx,delx)
    A = -2*eye(numx) +...
    diag(diag(eye(numx-1)),-1) +...
    diag(diag(eye(numx-1)),1);

    A(1,2) = 2; A(numx, numx-1) = 2;
    A(1,1) = -2+2*delx; A(numx,numx) = -2-2*delx;
    A = A/(delx^2);
end