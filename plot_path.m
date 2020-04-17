% plotting a piecewise linear path from the matrix X
% of step directions

function f = plot_path(X)
   
    m = size(X,2); % 10 
    d = size(X,1); % 2 
    Y = zeros(d,m+1);
    for i = 1:m
        y = zeros(d,m+1);
        y(:,i+1:end) = repmat(X(:,i),1,m-i+1);
        Y = Y + y;
    end
    
    if d == 2 
        plot(Y(1,:),Y(2,:)); % the steps of the path 
    else plot3(Y(1,:),Y(2,:),Y(3,:));
    end
 
end