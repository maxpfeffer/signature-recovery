function [X,cost] = recover_linear_transformation(S,m,options)
% Recovers the linear transformation of the axis signature with m pieces
% that results in the signature S.
%
% This method uses a trust-region solver from the manopt toolbox (cite as
% below) with a preliminary run of a BFGS-solver. If the number of linear
% pieces m is lower than or equal to the dimension d, the path will be
% recovered uniquely. The algorithm will be performed several times with
% different starting points to avoid local minima. If there are more pieces
% than the dimension, the algorithm will attempt to obtain the shortest
% such path.
%
% The initial iterate is X0 if it is provided. Otherwise, a random point is
% picked. To specify options whilst not specifying an initial iterate, give
% X0 as [] (the empty matrix).
%
% The two outputs 'x' and 'cost' are the last reached point on the manifold
% and its cost. 
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tries (5)
%       Number of tries with different starting points.
%   bfgs (true)
%       Set to false if no preliminary BFGS is required.
%   display (true)
%       Set to false to disable all outputs.
%   decreases (20)
%       If the shortest path is to be calculated, this is the number of
%       times lambda will be decreased.
%   startlambda (0.001)
%       If the shortest path is to be calculated, the is the initial value 
%       of lambda.
%   plot (false)
%       Set to true if the path should be plotted after recovery (only
%       possible if d = 2,3).
%
%   FURTHER options OF THE MANOPT TOOLBOX (THAT CAN REASONABLY BE ALTERED):
%
%   tolgradnorm (1e-10)
%       The algorithm terminates if the norm of the gradient drops below
%       this.
%   tolcost (1e-20)
%       The algorithm terminates if the value of the cost function drops
%       below this. For exact recovery, it is desired that this is very
%       small.
%   maxiter (1000)
%       The algorithm terminates if maxiter (outer) iterations were 
%       executed. 
%   minstepsize (1e-20)
%       The algorithm terminates if the stepsize is smaller than this
%       value. It is desired that even very small steps are possible.
%	maxinner (10*problem.M.dim() : the manifold's dimension)
%       Maximum number of inner iterations (for tCG).
%	Delta_bar (0.1*sqrt(problem.M.dim()))
%       Maximum trust-region radius. If you specify this parameter but not
%       Delta0, then Delta0 will be set to 1/8 times this parameter.
%   Delta0 (Delta_bar/8)
%       Initial trust-region radius. If you observe a long plateau at the
%       beginning of the convergence plot (gradient norm VS iteration), it
%       may pay off to try to tune this parameter to shorten the plateau.
%       You should not set this parameter without setting Delta_bar too (at
%       a larger value).
%   verbosity (0)
%       Integer number used to tune the amount of output the algorithm
%       generates during execution (mostly as text in the command window).
%       The higher, the more output. 0 means silent. 3 and above includes a
%       display of the options structure at the beginning of the execution.
%
% Please cite the Manopt paper:
%       @article{manopt,
%         author  = {Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
%         title   = {{M}anopt, a {M}atlab Toolbox for Optimization on Manifolds},
%         journal = {Journal of Machine Learning Research},
%         year    = {2014},
%         volume  = {15},
%         pages   = {1455--1459},
%         url     = {http://www.manopt.org},
%       }
% Please also cite the research paper:
%       @article{PfefferSeigalSturmfels2019,
%         author  = {Max Pfeffer and Anna Seigal and Bernd Sturmfels},
%         title   = {{Learning paths from signature tensors}},
%         doi     = {10.1137/18M1212331},
%         journal = {SIAM journal on matrix analysis and applications},
%         pages   = {394--416},
%         year    = {2019},
%         volume  = {40},
%         number  = {2},
%         issn    = {0895-4798},
%       }

C = piecewise_lin_sign(m);
d = size(S,1);
% If d >= m, the recovered path will be unique. Otherwise, we seek the
% shortest path
if d < m
    short_path = true;
    fprintf('The number of pieces is greater than the dimension. Looking for shortest path now.');
else
    short_path = false;
end

% Create the problem structure.
manifold = euclideanfactory(d,m);
problem.M = manifold;

% Set default options here
defaults.tries = 5;
defaults.bfgs = true;
defaults.display = true;
defaults.decreases = 20;
defaults.startlambda = 0.001;
defaults.plot = false;
defaults.tolgradnorm = 1e-10;
defaults.tolcost = 1e-20;
defaults.maxiter = 1000;
defaults.minstepsize = 1e-20;
defaults.maxinner = 10*problem.M.dim();
defaults.Delta_bar = 0.1*sqrt(problem.M.dim());
defaults.verbosity = 0;
    
% Merge defaults with user options, if any.
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(defaults, options);

function f = costfun(A,lambda)
    cos = trip_prod(A,C) - S; 
    if short_path
        f = lambda*path_len(A) + 0.5*dot(cos(:),cos(:));
    else
        f = 0.5*dot(cos(:),cos(:));
    end
end

function p = path_len(A)
    p = 0;
    for i = 1:m
        p = p + norm(A(:,i)); 
    end
end

function diff = diff_len(A)
    diff = zeros(d,m);
    for i = 1:m
        diff(:,i) = A(:,i)./norm(A(:,i),'fro'); 
    end
end
    
function g = eucgrad(A,lambda)
    X12 = nmodeproduct(nmodeproduct(C,A,1),A,2);
    X23 = nmodeproduct(nmodeproduct(C,A,2),A,3);
    X13 = nmodeproduct(nmodeproduct(C,A,1),A,3);
    res = trip_prod(A,C) - S;
    g12 = reshape(res,[d*d,d])'*reshape(X12,[d*d,m]);
    g23 = reshape(permute(res,[2 3 1]),[d*d,d])'*reshape(permute(X23,[2 3 1]),[d*d,m]);
    g13 = reshape(permute(res,[1 3 2]),[d*d,d])'*reshape(permute(X13,[1 3 2]),[d*d,m]);
    if short_path
        g = lambda*diff_len(A) + g12 + g23 + g13;
    else
        g = g12 + g23 + g13;
    end
end

function diff_form = diff2_len(A,U)
    diff2 = zeros(m*d,m*d);
    for j = 1:m
        k = (j-1)*d + 1;
        for i = 1:d
            diff2(k+i-1,k+i-1) = 1/norm(A(:,j),'fro') - (A(i,j)^2)/(norm(A(:,j),'fro')^3);
            for l = (i+1):d
                diff2(k+i-1,k+l-1) = - A(l,j)*A(i,j)/(norm(A(:,j),'fro')^3);
                diff2(k+l-1,k+i-1) = diff2(k+i-1,k+l-1);
            end
        end
    end
    diff_form = reshape(diff2*U(:),[d,m]);
end 

function h = euchess(A,U,lambda)
    X12 = nmodeproduct(nmodeproduct(C,A,1),A,2);
    X23 = nmodeproduct(nmodeproduct(C,A,2),A,3);
    X13 = nmodeproduct(nmodeproduct(C,A,1),A,3);
    X12U = nmodeproduct(X12,U,3);
    X23U = nmodeproduct(X23,U,1);
    X13U = nmodeproduct(X13,U,2);
    right = X12U + X23U + X13U;
    h12 = reshape(X12,[d*d,m])'*reshape(right,[d*d,d]);
    h23 = reshape(permute(X23,[2 3 1]),[d*d,m])'*reshape(permute(right,[2 3 1]),[d*d,d]);
    h13 = reshape(permute(X13,[1 3 2]),[d*d,m])'*reshape(permute(right,[1 3 2]),[d*d,d]);

    res = trip_prod(A,C) - S;
    X1U2 = nmodeproduct(nmodeproduct(C,A,1),U,2);
    X1U3 = nmodeproduct(nmodeproduct(C,A,1),U,3);
    X2U1 = nmodeproduct(nmodeproduct(C,A,2),U,1);
    X2U3 = nmodeproduct(nmodeproduct(C,A,2),U,3);
    X3U1 = nmodeproduct(nmodeproduct(C,A,3),U,1);
    X3U2 = nmodeproduct(nmodeproduct(C,A,3),U,2);
    h3 = reshape(res,[d*d,d])'*reshape(X1U2 + X2U1,[d*d,m]);
    h1 = reshape(permute(res,[2 3 1]),[d*d,d])'*reshape(permute(X2U3 + X3U2,[2 3 1]),[d*d,m]);
    h2 = reshape(permute(res,[1 3 2]),[d*d,d])'*reshape(permute(X1U3 + X3U1,[1 3 2]),[d*d,m]);
    
    h = h12' + h23' + h13' + h3 + h1 + h2;
    
    if short_path
        h = h + lambda*diff2_len(A,U);
    end
end 

if short_path
    fprintf('\n----------------------------------------------------- \n');
    Y = randn(d,m);
    for tr = 1:options.tries
        X0 = randn(d,m);
        lambda = options.startlambda;
        fprintf('Try  %i,  initial lambda  %d \n----------------------------------------------------- \n',tr,lambda);
        for it = 1:options.decreases

            % Define the problem cost function and its Euclidean gradient.
            problem.cost  = @(A) costfun(A,lambda);
            problem.egrad = @(A) eucgrad(A,lambda);
            problem.ehess = @(A,U) euchess(A,U,lambda);

            % Solve.
            if it == 1 && options.bfgs
                X = bfgs_method(problem.cost,problem.egrad,X0,options.tolgradnorm,0.1*options.maxiter,options.minstepsize,options.verbosity);
                fprintf('BFGS complete!  Cost: \t %d,  Path length: \t %d. \n',costfun(X,0),path_len(X));
            elseif it == 1
                X = X0;
            end
            X = trustregions(problem,X,options);
            fprintf('Step \t %i,  Cost: \t %d,  Path length: \t %d,  lambda: \t %d. \n',it,costfun(X,0),path_len(X),lambda);            

            % Decrease lambda.
            lambda = lambda/2;

            % If desired, plot the path for the given lambda (only for try 1).
            if tr == 1 && options.plot
                plot_path(X)
                hold on
            end
        end
        fprintf('----------------------------------------------------- \n');
        if tr == 1
            Y = X;
            if costfun(X,0) < 1e-10
                fprintf('Found new shortest path in try %i. \n',tr);
            end
        elseif costfun(X,0) < 1e-10 && path_len(X) < path_len(Y)
            Y = X;
            fprintf('Found new shortest path in try %i. \n',tr);
        end
    end
    X = Y;
    cost = costfun(X,0);
else
    fprintf('\n----------------------------------------------------- \n');
    for tr = 1:options.tries
        X0 = randn(d,m);
        fprintf('Try  %i \n----------------------------------------------------- \n',tr);

        % Define the problem cost function and its Euclidean gradient.
        problem.cost  = @(A) costfun(A);
        problem.egrad = @(A) eucgrad(A);
        problem.ehess = @(A,U) euchess(A,U);

        % Solve.
        if options.bfgs
            X = bfgs_method(problem.cost,problem.egrad,X0,options.tolgradnorm,0.1*options.maxiter,options.minstepsize,options.verbosity);
            fprintf('BFGS complete!  Cost: \t %d. \n',costfun(X));
        else
            X = X0;
        end
        X = trustregions(problem,X0,options);
        fprintf('Trust-Regions complete!  Cost: \t %d. \n',costfun(X));            

        % If desired, plot the path (only for try 1).
        if tr == 1 && options.plot
            plot_path(X)
            hold on
        end
        fprintf('----------------------------------------------------- \n');
        if costfun(X) < 1e-10
            fprintf('Found path in try %i. \n',tr);
            break
        end
    end
    cost = costfun(X);
end

end