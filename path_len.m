function p = path_len(A)
    m = size(A,2);
    p = 0;
    for i = 1:m
        p = p + norm(A(:,i),'fro'); 
    end
end