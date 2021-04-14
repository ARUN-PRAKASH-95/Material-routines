function C = t4_contr_t2(A,B)
%
% double contraction 
% C(i,j) = A(i,j,k,l)*B(k,l)
%
%[n,m] = size(A);

C=zeros(3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j) = C(i,j) + A(i,j,k,l)*B(k,l);
            end
        end
    end
end
