function [ A,B ] = test2( A_old,B_old,t,C,D,L,M,ALPHA,BETA,ROU,THETA,w )

    c = length(C);
    d = length(D);
    s = length(ALPHA);
    
    cvx_begin quiet
        variable A(c,s) binary
        variable B(d,s) binary
        
        p_old = price(A_old,B_old,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
        p     = price(A    ,B    ,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
        
        r = trans_price( B_old,B,D,ROU,THETA );
        
        minimize(r + p*t)
        subject to
            sum(A,2) == ones(c,1)
            sum(B,2) == ones(d,1)
            r + p*t  <= p_old*t
    cvx_end  

end

