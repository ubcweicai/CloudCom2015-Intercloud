function [ A ] = rand_loc( r, c )

    A = zeros(r,c);
    
    for i = 1:r
        A(i,randi([1,c])) = 1;
    end

end

