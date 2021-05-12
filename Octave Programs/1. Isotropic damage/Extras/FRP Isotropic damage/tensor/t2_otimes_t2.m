function C4 = t2_otimes_t2(A,B)
%
% C4(i,j,k,l) = A(i,j) * B(k,l)
%
[n,m] = size(A);
C4 = zeros(n,n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                C4(i,j,k,l) = A(i,j)*B(k,l);
            end;
        end;
    end;
end;
