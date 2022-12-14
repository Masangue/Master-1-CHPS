function [x,iter,err] = poisson1D(A,b,x0,alpha,maxiter,tol)
    iter = 0;
    res = b - A*x0;
    err(iter+1) = norm(res)/norm(b);
    x=x0;
    while (err(iter+1) > tol) && (iter < maxiter)
        iter = iter + 1;
        x = x + alpha*(res);
        res = b - A*x;
        err(iter+1) = norm(res)/norm(b);
    end
endfunction

n = 10;
A = rand(n,n);
xex = rand(n,1);
b = A*xex;
x0 = ones(n,1);
maxiter = 100000;
tol = 1e-5;

vp = spec(A);
A = 2*diag(ones(n,1)) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
disp(A);
disp(vp);
lambdamax = max(vp);
lambdamin = min(vp);
alpha = 2/(lambdamax+lambdamin);

[x,iter,err] = poisson1D(A,b,x0,alpha,maxiter,tol);

disp(x,xex,iter,err);

quit;