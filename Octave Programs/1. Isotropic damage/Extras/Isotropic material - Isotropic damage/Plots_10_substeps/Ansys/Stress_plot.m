A=load('Stress.txt');
Stress = A(:,2);
Time  = A(:,1);


B = load('e11.txt');
strain_11 = B(:,2);


C = load('e22.txt');
strain_22 = C(:,2);

D = load('e33.txt');
strain_33 = D(:,2);


figure(1)
plot(Time,strain_11) 
xlabel('time')
ylabel('eps11')


figure(2)
plot(strain_11,Stress,'-o','color','b') 
legend('\sigma_{11}','Ref.','Location','NorthWest')
xlabel('eps11')
ylabel('sig11')


figure(3)
plot(strain_11,strain_22,'-o','color','b')
hold on
plot(strain_11,strain_33,'-o','color','g')
legend('\epsilon_{22}','\epsilon_{33}','Ref.','Location','NorthEast')
xlabel('eps11')
ylabel('eps22, eps33')