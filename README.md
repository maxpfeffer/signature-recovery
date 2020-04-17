# signature-recovery

This method is made to recover linear paths from their third order signature. If you have a signature S and you want to recover the linear transformation of the axis path, you only have to enter the number of pieces m. You can generate a random S by

Z = randn(d,m); <br>
C = piecewise_lin_sign(m); <br>
S = trip_prod(Z,C); <br>

In order to recover the matrix, run

[X,cost] = recover_linear_transformation(S,m)

The two outputs 'x' and 'cost' are the last reached point on the manifold and its cost. If m <= d, the algorithm will recover the path. Because this does not work every time, it will try 5 times by default until the cost function is small enough. If m > d, the algorithm will look for the shortest path out of the 5 tries. 

The method uses a trust-region solver from the manopt toolbox (cite as below) with a preliminary run of a BFGS-solver. If the number of linear pieces m is lower than or equal to the dimension d, the path will be recovered uniquely. The algorithm will be performed several times with different starting points to avoid local minima. If there are more pieces than the dimension, the algorithm will attempt to obtain the shortest such path. The initial iterate is X0 if it is provided. Otherwise, a random point is picked. To specify options whilst not specifying an initial iterate, give X0 as [] (the empty matrix).

You can alter a lot of options. The options structure is used to overwrite the default values. All options have a default value and are hence optional. To force an option value, pass an options structure with a field options.optionname, where optionname is one of the following and the default value is indicated between parentheses:

{{plainlist
   * tries (5) <br>
       Number of tries with different starting points.
   * bfgs (true) <br>
       Set to false if no preliminary BFGS is required.
   * display (true) <br>
       Set to false to disable all outputs.
   * decreases (20) <br>
       If the shortest path is to be calculated, this is the number of
       times lambda will be decreased.
   * startlambda (0.001) <br>
       If the shortest path is to be calculated, the is the initial value 
       of lambda.
   * plot (false) <br>
       Set to true if the path should be plotted after recovery (only
       possible if d = 2,3).

   Further options of the manopt toolbox (that can reasonably be altered):

   * tolgradnorm (1e-10) <br>
       The algorithm terminates if the norm of the gradient drops below
       this.
   * tolcost (1e-20) <br>
       The algorithm terminates if the value of the cost function drops
       below this. For exact recovery, it is desired that this is very
       small.
   * maxiter (1000) <br>
       The algorithm terminates if maxiter (outer) iterations were 
       executed. 
   * minstepsize (1e-20) <br>
       The algorithm terminates if the stepsize is smaller than this
       value. It is desired that even very small steps are possible.
	 maxinner (10*problem.M.dim() : the manifold's dimension)
       Maximum number of inner iterations (for tCG).
   * Delta_bar (0.1*sqrt(problem.M.dim())) <br>
       Maximum trust-region radius. If you specify this parameter but not
       Delta0, then Delta0 will be set to 1/8 times this parameter.
   * Delta0 (Delta_bar/8) <br>
       Initial trust-region radius. If you observe a long plateau at the
       beginning of the convergence plot (gradient norm VS iteration), it
       may pay off to try to tune this parameter to shorten the plateau.
       You should not set this parameter without setting Delta_bar too (at
       a larger value).
   * verbosity (0) <br>
       Integer number used to tune the amount of output the algorithm
       generates during execution (mostly as text in the command window).
       The higher, the more output. 0 means silent. 3 and above includes a
       display of the options structure at the beginning of the execution.

 Please cite the Manopt paper: <br>
       @article{manopt, <br>
         author  = {Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.}, <br>
         title   = {{M}anopt, a {M}atlab Toolbox for Optimization on Manifolds}, <br>
         journal = {Journal of Machine Learning Research}, <br>
         year    = {2014}, <br>
         volume  = {15}, <br>
         pages   = {1455--1459}, <br>
         url     = {http://www.manopt.org}, <br>
       } <br>
 Please also cite the research paper: <br>
       @article{PfefferSeigalSturmfels2019, <br>
         author  = {Max Pfeffer and Anna Seigal and Bernd Sturmfels}, <br>
         title   = {{Learning paths from signature tensors}}, <br>
         doi     = {10.1137/18M1212331}, <br>
         journal = {SIAM journal on matrix analysis and applications}, <br>
         pages   = {394--416}, <br>
         year    = {2019}, <br>
         volume  = {40}, <br>
         number  = {2}, <br>
         issn    = {0895-4798}, <br>
       }
