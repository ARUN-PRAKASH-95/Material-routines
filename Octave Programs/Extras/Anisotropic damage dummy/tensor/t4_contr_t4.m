function C = t4_contr_t4(A,B)
%
% double contraction 
% C(i,j) = A(i,j,k,l)*B(k,l,i,j)
%
%[n,m] = size(A);

C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l) = C(i,j,k,l) + A(i,j,k,l)*B(k,l,i,j);
            end
        end
    end
end