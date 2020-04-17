function Y = trip_prod(A,X)
     Y = nmodeproduct(X,A,1);
     Y = nmodeproduct(Y,A,2);
     Y = nmodeproduct(Y,A,3);
end