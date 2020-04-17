function [x,cost,gradnorm] = bfgs_method(f,grad,x0,tol,iterations,minstepsize,verbosity)

cost = zeros(iterations+1,1);
gradnorm = zeros(iterations+1,1);

[d,m] = size(x0);
cost(1) = f(x0);
gradnorm(1) = norm(grad(x0),'fro');

i = 1;
x = x0;
c = 0.0001;
apprhess = eye(d*m);
while gradnorm(i) > tol && i < iterations+1
  
    g = grad(x);
    p = reshape(-apprhess*g(:),[d,m]);
    
    super = 0;
    alpha = 1;
    while f(x + alpha*p) - f(x) > c*alpha*dot(g(:),p(:)) && alpha > minstepsize
        alpha = alpha/2;
    end
    
    if alpha < minstepsize
        beta = 1;
        while f(x - beta*g) - f(x) > -c*beta*dot(g(:),g(:)) && beta > minstepsize
            beta = beta/2;
        end
        if beta < minstepsize
            disp('Stepsize smaller than allowed!');
            break
        else
            x = x - beta*g;
            apprhess = eye(d*m);
            s = -beta*g(:);
        end
    else
        super = 1;
        x = x + alpha*p;
        s = alpha*p(:);
    end
    
    cost(i+1) = f(x);
    gradnorm(i+1) = norm(grad(x0),'fro');
    i = i+1;
    if verbosity
        fprintf('Step: %i, f: %d, gradnorm: %d , super: %i \n',i,f(x),norm(grad(x0),'fro'),super);
    end
        
    y = grad(x) - g;
    y = y(:);
    if (s'*y)^2 < 1e-16
        apprhess = eye(d*m);
    else
        apprhess = apprhess + (s'*y + y'*apprhess*y)*(s*s')/((s'*y)^2) - (apprhess*y*s' + s*y'*apprhess)/(s'*y);
    end
        
end

end


