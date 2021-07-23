function  UserDefinedReturn=SIMULATION(TARGETSTRAIN,PARA,SIGMA_PRE,RHO)
% ARC LENGTH METHOD TO CONTROL LOADING PATH
% LITERATURE: MASTERTHESIS ANDREAS SEUPEL 2013; 
%SIGMA_PRE: LOADING STRESS COMPONENT S11, S22==S33=RHO*S11
%numerical tangent
ttype = 0;

%% INITIATE (USER DEFINED) STATE VARIABLES
STATEV      =zeros(6,1);



%% PARAMETER OF INCREMENTATION (ARC LENGTH INCREMENT DELTAS)

DELTAS          =0.001;
DELTASMAX       =0.001;
DELTASMIN       =0.000001;

%% TOLERANCE OF THE FORCE (STRESS) EQUILIBRIUM RESIDUAL
 TOL         =1e-6;
 
 %% MAXIMALE ITERATIONS DURING NEWTON METHOD
 MAX         =20;

%%COMPONENT OF THE STRAIN TENSOR WHICH SHOULD REACH THE TARGET STRAIN TARGETSTRAIN
COMPONENT =1;

%% ESSENTIAL BOUNDARY CONDITIONS HAVE TO BE MARKED BY "1"
EPSMATRIXGES       	=[0,0.0,0.0,0,0,0]';

% ARRAY CONTAINING THE HISTORY OF STRAIN TENSOR 
EPSMATRIXSIM        	=[zeros(6,1)];
% INTERNAL VECTOR TO MARK UP THE BC'S
SHORTEN         	=[];
% REGISTRATION OF BOUNDARIES
for i=1:length(EPSMATRIXGES)
    
   if EPSMATRIXGES(i)==0
       
   else
        SHORTEN= horzcat(SHORTEN,i);
   end
end
% CURRENT STRAIN INCREMENT TENSOR
EPSMATRIX       =zeros(6,1);

%% START OF ARC LENGTH METHOD %%
%% DECLARATION OF: INTERNAL STRESS TENSOR (STRESSSIM) AND EXTERNAL; PRESCRIBED STRESS (FEXT)

STRESSSIM     =[zeros(6,1)];
STATEVSIM     =[zeros(6,1)];
FEXT            =[SIGMA_PRE,RHO,RHO,0,0,0]';

%% UPDATE OF STRAIN INCREMENT 
DELTAEPSMATRIX  =zeros(6,1);

%% LOADING PARAMETER LAMBDA
LAMBDA0         =1;
LAMBDA          =0;

%% PREDICTOR STEP, M=1
% CALCULATE INITIAL STIFFNESS K, UMAT=MATERIAL ROUTINE
% UMAT RETURNS: MATERIAL STIFFNESS K, STRESS TENSOR FOR GIVEN STRAIN INCREMENT STRESSITERATION, 
%INTERNAL OR USER DEFINED VARIABLES STATEV, CONTROL VARIABLE PNEW: PNEW==0 -> NO PROBLEM IN UMAT; PNEW==1 -> PROBLEM IN UMAT, ARC LENGTH STEP IS REDUCED
%UMAT NEEDS: STRAIN INCREMENT DELTAEPSMATRIX, OLD STRESS STATE STRESSSIM(:,end), OLD STATEV-VARIABLES, MATERIAL OR MODEL PARAMETERS PARA
PARA = inputmat();

[K,STRESSITERATION,STATEVITER,PNEW]     =  UMAT(EPSMATRIX,STRESSSIM(:,end),STATEV,PARA,ttype);

%% INCORPORATING THE BC'S IN THE SYSTEM OF EQUATIONS
[ FEXT ] 	= DeleteRow( FEXT,SHORTEN );
       
[ K ]  		= DeleteRow( K,SHORTEN );
[ K ]  		= DeleteColumn( K,SHORTEN );

%% PREDICTOR VALUE OF FIRST STRAIN INCREMENT
DUSTAR      	=linsolve(K,FEXT);
[ DUSTAR ]  	= CreateRow(DUSTAR,SHORTEN );

[DUSTARMAX,INDEX]   =max(abs(DUSTAR));
DUSTARMAX           =DUSTAR(INDEX);
%% SCALING OF THE PREDICTOR
ALPHA       =DELTAS^2/(DUSTARMAX^2);

DU          =sqrt(ALPHA)*DUSTAR;
DLAMBDA     =sqrt(ALPHA)*LAMBDA0;

%% CONTROL VARIABLE
PNEW=0;
TARGETSTRAIN = 0.2;
%% INCREMENTATION OVER LOAD STEPS
%THE LOOP IS INTERRUPTED IF TARGET STRAIN IS EXCEEDED
while max(abs(EPSMATRIX))<abs(TARGETSTRAIN)
    

% CURRENT  STRESS RESIDUAL
    ERRORGG     =2*TOL;
% ITERATION COUNTER
    COUNTER     =0;

    
%% CORRECTIONS STEPS OF NEWTON METHOD
    while ERRORGG>TOL && PNEW==0

% UMAT RETURNS: MATERIAL STIFFNESS K, STRESS TENSOR FOR GIVEN STRAIN INCREMENT STRESSITERATION, 
%INTERNAL OR USER DEFINED VARIABLES STATEV, CONTROL VARIABLE PNEW: PNEW==0 -> NO PROBLEM IN UMAT; PNEW==1 -> PROBLEM IN UMAT, ARC LENGTH STEP IS REDUCED
%UMAT NEEDS: STRAIN INCREMENT  -->DU<--, OLD STRESS STATE STRESSSIM(:,end), OLD STATEV-VARIABLES, MATERIAL OR MODEL PARAMETERS PARA

        [K,STRESSITERATION,STATEVITER,PNEW]     =UMAT(EPSMATRIX+DU,STRESSSIM(:,end),STATEV,PARA,ttype);
        
% DENOTATION AS IN ALGORITHMN 4.2 IN MASTER THESIS ANDREAS SEUPEL 2013
        [ K ]       = DeleteRow( K,SHORTEN );
        [ K ]       = DeleteColumn( K,SHORTEN );
        [ FINT ]    = DeleteRow( STRESSITERATION,SHORTEN );
        
        %SOLVE FOR DDUG
        DDUG    =linsolve(K,-1*(FINT-(LAMBDA+DLAMBDA)*FEXT));
        
        %SOLVE FOR DDUF
        DDUF    =linsolve(K,FEXT);
        
        %FIND MAXIMAL STRAIN COMPONENT
        DUMAXVEK        =DeleteRow( DU,SHORTEN );
        [DUMAX,INDEX]   =max(abs(DUMAXVEK));
        DUMAX           =DUMAXVEK(INDEX);
        DDUGMAX         =DDUG(INDEX);
        DDUFMAX         =DDUF(INDEX);
%         
%         DUMAXVEK        =DeleteRow( DU,SHORTEN );
% %         [DUMAX,INDEX]   =max(abs(DUMAXVEK));
%         INDEX=1;
%         DUMAX           =abs(DUMAXVEK(INDEX));
%         DDUGMAX         =abs(DDUG(INDEX));
%         DDUFMAX         =abs(DDUF(INDEX));
        
        
        %SOLVE CONDITION EQUATION
        C0              =DUMAX^2+2*DUMAX*DDUGMAX+DDUGMAX^2-DELTAS^2;
        C1              =2*DDUFMAX*(DUMAX+DDUGMAX);
        C2              =DDUFMAX^2;
        P               =C1/C2;
        Q               =C0/C2;
        DDLAMBDA1       =-P/2+sqrt(P^2/4-Q);
        DDLAMBDA2       =-P/2-sqrt(P^2/4-Q);
        % UPDATE AND DECISION OF FURTHER LOADING PATH
        DU1             =DU+DDUG+DDUF*DDLAMBDA1;
        G1              =DU1'*DU;
        DU2             =DU+DDUG+DDUF*DDLAMBDA2;
        G2              =DU2'*DU;
        if G1>G2
            DDLAMBDA    =DDLAMBDA1;
            [ DU ]      = CreateRow(DU1,SHORTEN );
        else
            DDLAMBDA    =DDLAMBDA2;
            [ DU ]      = CreateRow(DU2,SHORTEN );
        end
       DLAMBDA         =DLAMBDA+DDLAMBDA;
%% CALCULATE NORM (ERROR ) OF RESIDUAL
       ERRORGG         =max(norm(FINT-(LAMBDA+DLAMBDA)*FEXT),norm(DDUG+DDUF*DDLAMBDA));
        COUNTER         =COUNTER+1;
        if COUNTER>=MAX
            PNEW        =1;
            disp('ABORT BECAUSE OF DIVERGENCE')
        end
    end

%% CONTROL OF CONVERGENCE
    if PNEW==0
        STRESSSIM         =horzcat(STRESSSIM,STRESSITERATION);
        STATEV              =STATEVITER;
        STATEVSIM     =horzcat(STATEVSIM,STATEV);
        
        if DELTAS<DELTASMAX
           
            DELTAS              =DELTAS*1.25;
            
            if DELTAS > DELTASMAX
                DELTAS          =DELTASMAX;
            end

        end
        %% UPDATE AND NEW PREDICTOR STEP
            EPSMATRIX           =EPSMATRIX+DU;
            LAMBDA              =LAMBDA+DLAMBDA;
            DUSTAR              =DU;
            EPSMATRIXSIM        =horzcat(EPSMATRIXSIM,EPSMATRIX);
            [DUSTARMAX,INDEX]   =max(abs(DUSTAR));
            DUSTARMAX           =DUSTAR(INDEX);
            % SCALING
            ALPHA               =DELTAS^2/(DUSTARMAX^2);

            DU                  =sqrt(ALPHA)*DUSTAR;
            DLAMBDA             =sqrt(ALPHA)*LAMBDA0;
            

    else

            DELTAS                  =DELTAS/2;
            ALPHA                   =DELTAS^2/(DUSTARMAX^2);

            DU                      =sqrt(ALPHA)*DUSTAR;
            DLAMBDA                 =sqrt(ALPHA)*LAMBDA0;

            PNEW                    =0;                  
            disp('STEP SIZE HAS BEEN DECREASED')
    end
    
    if DELTAS<DELTASMIN
        disp('SMALLEST STEPSIZE REACHED! NO SOLUTION FOUND!')
        break  
    end


end

  
  
  
%% RETURN
%% USER DEFINED!
% STRAIN IN X1-DIRECTION
EPS1 	=EPSMATRIXSIM(1,1:end);

% STRESS IN X1-DIRECTION
S1 	=STRESSSIM(1,1:end);
D1 	=STATEVSIM(1,1:end);

% STRAIN IN X2-DIRECTION
EPS2 	=EPSMATRIXSIM(2,1:end);

% STRESS IN X2-DIRECTION
S2 	=STRESSSIM(2,1:end);
D2 	=STATEVSIM(2,1:end);

% STRAIN IN X3-DIRECTION
EPS3 	=EPSMATRIXSIM(3,1:end);

% STRESS IN X3-DIRECTION
S3 	=STRESSSIM(3,1:end);
D3 	=STATEVSIM(3,1:end);

UserDefinedReturn=horzcat(EPS1',S1');
%% PLOT
subplot(3,2,1)
plot(EPS1 ,S1,'marker','x')
title('\epsilon_{11} vs. \sigma_{11}')
subplot(3,2,3)
plot(EPS2 ,S2,'marker','x')
title('\epsilon_{22} vs. \sigma_{22}')
subplot(3,2,5)
plot(EPS3 ,S3,'marker','x')
title('\epsilon_{33} vs. \sigma_{33}')

subplot(3,2,2)
plot(EPS1 ,D1,'marker','x')
title('\epsilon_{11} vs. d_{1}')
subplot(3,2,4)
plot(EPS2 ,D2,'marker','x')
title('\epsilon_{22} vs. d_{2}')
subplot(3,2,6)
plot(EPS2 ,D3,'marker','x')
title('\epsilon_{33} vs. d_{3}')

end

