%==========================================================================
%   drive.m
%
%   driver for uni-axial isothermal stress tests 
% 
%==========================================================================
% 
%   coded by: B. Kiefer 15 Nov 2011
%   analytical solution added: S. Prueger 10 Dec 2018
%
%   comments: --> numerical tangent moduli computation for ttype = 1
%
%==========================================================================
clear all
close all
clc
addpath('tensor/'); 
addpath('analyt_sol/'); 

%--------------------------------------------------------------------------
%   define loading
%--------------------------------------------------------------------------
% 1: linear ramping of load
% 2: half cycle with linear load change
% 3: linear loading and unloading, start and end point different
% 4: full cycle with linear load change
% 5: two cycles
ltype=1;
if ltype==1
    t=[0 10];
    lam=[0 0.1];
elseif ltype==2
    t=[0 5 10];
    lam=[0 1.59155e-3];
elseif ltype==3
    t=[0 5 10];
    lam=[0 0.0015 0.0000];
elseif ltype==4
    t=[0 2.5 7.5 10];
    lam=[0 n_ampl*sigma_y0/2/mu/(1-q_el) -n_ampl*sigma_y0/2/mu/(1-q_el)];
elseif ltype==5
    t=[0 2.5 7.5 12.5 17.5 22.5 25.0];
    lam=[0 0.005 -0.005];
end

%--------------------------------------------------------------------------

% prescribed load/time step
dt=0.1;
% start and end-time of loading, time-scale, no. of steps
ta=t(1);
te=t(end);
time=ta:dt:te;
steps=size(time,2)-1;
e11=loading(ltype,dt,t,lam);

%-------------------------------------------------------------------------

% initialize strains, temperature and internal variables
epsbar=zeros(2,1);
sdv = zeros(5,steps);


% initialise quantities for post-processing
s11=zeros(1,steps);
eps22=zeros(1,steps); eps12=zeros(1,steps);

% tolerance and maximum no. of iterations for Newton iteration
tol=1e-4;
maxit=100;
ttype = 0; % 0: analytical, 1: numerical tangent moduli computation

% initialize waitbar
wb=waitbar(0,'computation in progress...');

for n=1:steps
    
    % display waitbar
    waitbar(n/steps)
    
    % display current time step
    disp(['n = ', num2str(n)])
    
    % initialize rnorm in order to enter the while-loop and set no. of
    % iterations to zero
    sbar=ones(2,1);
    iter=0;
    
    while norm(sbar) > tol % check convergence
        iter = iter+1;
        if iter > maxit
            error(['No convergence after ', num2str(maxit), ... 
                ' global iterations'])
        end
        
        % 1.) total deformation
         epsilon(1,1) = e11(n+1);
         
         epsilon(2:3,1) = epsbar;
         
        % 2.) constitutive law: algorithmic stresses and moduli 
        [s,A,sdvup]=subroutine_planestress(epsilon,sdv(:,n),ttype);
        
        % 3.) partitioning
        sbar=partition_planestress(s);
        Abar=partition_planestress(A);

        % 4.) update of lateral deformations
        epsbar=epsbar-Abar\sbar;

        % display convergence
        disp(['|sbar| = ', num2str(norm(sbar))])

    end % while
    
    
    % update of internal variables after obtaining convergence
    sdv(:,n+1) = sdvup;
    
    % store quantities for post-processing
    s11(n+1)=s(1);
    eps22(n+1)=epsilon(2); eps12(n+1)=epsilon(3); 

end % for
fprintf('********************\n')
close(wb)
%
s
%fprintf('eps11 %f\n', e11);
%fprintf('eps22 %f\n', eps22);
%fprintf('eps33 %f\n', eps33);

data = [time; e11; s11; eps22; ];

fileID = fopen('output_planestress.txt','w');
fprintf(fileID,'%2s %15s %15s %15s\n','n','Epsilon_{11}','Sigma_{11}','Epsilon_{22}');
fprintf(fileID,'%25s\n','');
fprintf(fileID,'%3.1f %15.4f %15.5f %15.4f\n',data);
fclose(fileID);

figure(1)
%subplot(2,1,1)
plot(time,e11)
xlabel('time')
ylabel('e11')


% plot stress-strain response
figure(2);
 %subplot(2,1,2)
 %plot(e11,s11, paperX,paperY, '-.rx')
plot(e11,s11,'or-')
xlabel('\epsilon_{11}')
ylabel('\sigma_{11}')
 

 
 % plot lateral strains
figure(3)
plot(e11,eps22,'ob-')
hold on
plot(e11,eps12,'og-')
xlabel('\epsilon_{11}')
ylabel('\epsilon_{22},\epsilon_{12}')

% plot damage d1
figure(4)
plot(e11,sdv(1,:),'or-')
legend('d1','Ref.','Location','West')
xlabel('\epsilon_{11}')
ylabel('d1')



% plot damage d2
figure(5)
plot(e11,sdv(2,:),'or-')
legend('d2','Ref.','Location','West')
xlabel('\epsilon_{11}')
ylabel('d2')

% plot damage d3
figure(6)
plot(e11,sdv(3,:),'or-')
legend('d3','Ref.','Location','West')
xlabel('\epsilon_{11}')
ylabel('d3')

figure(7)
plot(e11,sdv(4,:),'r-')
legend('F_f','Ref.','Location','West')
xlabel('\epsilon_{11}')
ylabel('F_f')


figure(8)
plot(e11,sdv(5,:),'r-')
legend('F_m','Ref.','Location','West')
xlabel('\epsilon_{11}')
ylabel('F_m')