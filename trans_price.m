function [ r ] = trans_price( B_old,B,D,ROU,THETA )

    T = B - B_old;

    R1 = pos(T);
    r1 = D' * R1 .* ROU';
    R2 = pos(-T);
    r2 = D' * R2 .* THETA';
    r = r1 + r2;
    r = sum(r);

end

