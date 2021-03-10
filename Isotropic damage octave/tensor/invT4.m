function [Ainv] = invT4(A)
%invT4.m inverts 4th order tensor with appropriate mapping
% coded by: B. Kiefer 9 Feb 2012
A6x6=kelvintm(A);
Ainv6x6=inv(A6x6);
Ainv=kelvinmt(Ainv6x6);
end

