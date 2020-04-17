function C = piecewise_lin_sign(m);
C = zeros(m,m,m);

for i = 1:m
    for j = i:m
        for k = j:m
            x = [i j k];
            [a,b]=hist(x,unique(x));
            s = factorial(3);
            for l = 1:length(a)
                s = s/factorial(a(l));
            end
            C(i,j,k) = s/factorial(3);
        end
    end
end

end