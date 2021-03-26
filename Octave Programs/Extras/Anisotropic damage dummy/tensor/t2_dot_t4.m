function C = t2_dot_t4(A,B)
%


C=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C(i,j,k,l) = C(i,j,k,l) + A(k,i)*B(i,j,k,l);
            end
        end
    end
end