function [ p ] = price( A,B,C,D,L,M,ALPHA,BETA,ROU,THETA,w )

    c = length(C);
    d = length(D);
    s = length(ALPHA);

    Pc = sum(C'*A .* ALPHA');
    Pd = sum(D'*B .*  BETA');
    
    Pn1 = zeros(1,s);
    for i = 1:c
        for j = 1:d
            if (L(i,j) == 1)
                temp = zeros(c,d); temp(i,j) = 1;
                cs = (temp*w*D)'*A;
                ds = (sum(temp*(-w)).*D')*B;
                csds = cs+ds;
%                 Pn1 = Pn1 + pos(csds) .* ROU' + pos(-csds).* THETA';
                Pn1 = Pn1 + max(csds, 0) .* ROU' + max(-csds, 0).* THETA';
            end
        end
    end
    
    Pn2 = zeros(1,s);
    LM = M > 0;
    for i = 1:c
        for j = 1:c
            if (LM(i,j) == 1)
                temp = zeros(c); temp(i,j) = 1;
                DD = (M(i,:))';
                c1 = (temp*DD)'*A;
                c2 = (sum(temp*(-1)).*DD')*A;
                csds = c1+c2;
%                 Pn2 = Pn2 + pos(csds) .* ROU' + pos(-csds).* THETA';
                Pn2 = Pn2 + max(csds, 0) .* ROU' + max(-csds, 0).* THETA';
            end
        end
    end

    Pn = sum(Pn1)+sum(Pn2);
    
    p = Pc+Pd+Pn;

end

