function C = t4_dot_t2(A,B)
%
% dot product of 4th order and 2nd order tensor

C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l) = C(i,j,k,l) + A(i,j,k,l)*B(l,i);
            end
        end
    end
end