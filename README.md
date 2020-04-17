# signature-recovery

This method is made to recover linear paths from their third order signature. If you have a signature S and you want to recover the linear transformation of the axis path, you only have to enter the number of pieces m. You can generate a random S by

Z = randn(d,m);
C = piecewise_lin_sign(m);
S = trip_prod(Z,C);

In order to recover the matrix, run

[X,cost] = recover_linear_transformation(S,m)

If m <= d, the algorithm will recover the path. Because this does not work every time, it will try 5 times by default until the cost function is small enough. If m > d, the algorithm will look for the shortest path out of the 5 tries. You can alter a lot of options, for example

options.maxiter = 100; (The manopt code will now only do 100 iterations instead of 1000)
