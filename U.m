function [ UD ] = U( U_min, U_max, rand_dev, r, c )
%Uniform Distribution over [U_min,U_max]

    UD = ((2 * rand_dev * rand(r,c)) - rand_dev + 1) .* ... 
         ((U_max - U_min) * rand(r,c)) ... [1-rand_dev, 1+rand_dev]
       + ((2 * rand_dev * rand(r,c)) - rand_dev + 1)  * ... 
           U_min;
       
%     UD = (U_max - U_min) * rand(r,c) + U_min;  

end