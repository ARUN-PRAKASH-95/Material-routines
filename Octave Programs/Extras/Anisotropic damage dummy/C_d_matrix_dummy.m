

C  = eye(6,6);
C(1,2) = 1;
C(1,3) = 1;
C(2,1) = 1;
C(2,3) = 1;
C(3,1) = 1;
C(3,2) = 1;

C  = kelvinmt(C)


D  = zeros(3,3);
D(1,1) = 2;
D(2,2) = 3;
D(3,3) = 4;

DC = t2_dot_t4(D,C)

DC = kelvintm(DC);