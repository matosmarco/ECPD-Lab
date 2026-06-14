% Build Matrix W
% Idea: initialize the matrix to 0
function W = build_w(A,b,C,H)
W = zeros(H,H);
for ii=1:H
    for jj = 1:ii
        if ii == jj
            W(ii,jj) = C*b;
        else
        W(ii,jj) = C*A^(ii-jj)*b;
        end

    end
end
end