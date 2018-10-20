function [ A,B ] = test( C,D,L,M,ALPHA,BETA,ROU,THETA,w,x_C,x_D )

    c = length(C);
    d = length(D);
    s = length(ALPHA);
    
    cvx_begin quiet
    
        variable A(c,s) binary
        variable B(d,s) binary
        
        p = price(A,B,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
        
        minimize(p)
        subject to
            sum(A,2) == ones(c,1)
            sum(B,2) == ones(d,1)
            sum(A,1) <= x_C*ones(1,s)
            sum(B,1) <= x_D*ones(1,s)
    cvx_end  

end

