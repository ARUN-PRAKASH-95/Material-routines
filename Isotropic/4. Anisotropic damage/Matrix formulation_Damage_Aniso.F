       subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
	
	
#include "impcom.inc"
#include "ansysdef.inc"
      integer :: matId, elemId, kDomIntPt, kLayer, kSectPt, ldstep,
     &           isubst,keycut, nDirect,nShear,ncomp,nStatev,nProp


      double precision :: Time,    dTime,   Temp,    dTemp,
     &                    sedEl,   sedPl,   epseq,   epsZZ
                  
      double precision :: stress  (ncomp  ), ustatev (nStatev), 
     &                 dsdePl  (ncomp,ncomp), Strain  (ncomp  ), 
     &                 dStrain (ncomp  ), epsPl   (ncomp  ),
     &                 prop    (nProp  ), coords  (3), 
     &                 defGrad (3,3),     defGrad_t(3,3),  
     &                 tsstif  (2) 
	
	
      EXTERNAL         usermat3d
      EXTERNAL         myuserfunc
      integer :: myvar
	
	
	
      double precision :: var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8

      integer :: iott,wrinqr
      external wrinqr
c
c*************************************************************************
c
      iott = wrinqr(WR_OUTPUT)
      write(iott,*) ' ************************************************'
      write(iott,*) ' * #DEBUG# UserMatLib.dll USERMAT               *'
      write(iott,*) ' ************************************************'
      myvar = 99
      call myuserfunc(myvar)
	
	

	if (ncomp >= 4) then
	
		
		call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)


      end if
      return
	end
	
	
**************    Subroutine for 3d elements (Linear elastic isotropic model)   ********************
      subroutine usermat3d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
	
	
#include "impcom.inc"
#include "ansysdef.inc"
c
      integer :: matId, elemId, kDomIntPt, kLayer, kSectPt, ldstep,
     &           isubst,keycut, nDirect,nShear,ncomp,nStatev,nProp


      double precision :: Time,    dTime,   Temp,    dTemp,
     &                    sedEl,   sedPl,   epseq,   epsZZ
                  
      double precision :: stress  (ncomp  ), ustatev (nStatev), 
     &                 dsdePl  (ncomp,ncomp), Strain  (ncomp  ), 
     &                 dStrain (ncomp  ), epsPl   (ncomp  ),
     &                 prop    (nProp  ), coords  (3), 
     &                 defGrad (3,3),     defGrad_t(3,3),  
     &                 tsstif  (2)
	

	
	
********* USER DEFINED VARIABLES (FOR CALCULATION PURPOSE) **********	
	
      integer :: mcomp
      double precision :: HALF, ONE, TWO, ZERO,
     &                    THREE
      
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 mcomp      = 6
     &                 )
	
      EXTERNAL          vzero, vmove, get_ElmData

   
      double precision  dsdeEl(mcomp,mcomp), sigma(mcomp),
     &                  sigi (mcomp), C_T(mcomp,mcomp),
     &                  dsdeEl_d(mcomp,mcomp), T_strain(mcomp),
     &                  id(3,3), C_T_1(mcomp,mcomp),
     &                  d_C_d_d1(mcomp,mcomp), C_T_1_a(mcomp), 
     &                  C_T_1_b(mcomp), C_T_2(mcomp,mcomp), 
     &                  d_C_d_d2(mcomp,mcomp), C_T_2_a(mcomp),
     &                  C_T_2_b(mcomp),C_T_3(mcomp,mcomp), 
     &                  d_C_d_d3(mcomp,mcomp), C_T_3_a(mcomp),
     &                  C_T_3_b(mcomp)
                  

      
      integer             i,j,k,l,n,P
      double precision    young, posn, mu,  
     &                    twoG, lambda, L_c,k1,k2,k3,
     &                    xK, yield,G_c_1,G_c_2,G_c_3,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f,
     &                    sig_22_f,sig_33_f,d1,d2,d3,d1_new,d2_new
      
      double precision    d3_new,F_f,F_m,F_z,term1,term2
   
      

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

      
      integer :: iott,wrinqr
      external wrinqr
c
c*************************************************************************
c
      iott = wrinqr(WR_OUTPUT)
      write(iott,*) ' ************************************************'
      write(iott,*) ' * #DEBUG# UserMatLib.dll USERMAT1D             *'
      write(iott,*) ' ************************************************'

      keycut   = 0
	
      id = transpose(reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), 
     & shape(id)))


 
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      twoG     = young / (ONE+posn)
	yield    = 200
	G_c_1    = 10e2
      G_c_2    = 10e2
      G_c_3    = 10e2
      L_c      = 1**(1/3)
      
      
      eps_11_f_t  = yield/young 
      eps_11_f_c  = yield/young 
      eps_22_f_t  = yield/young
      eps_22_f_c  = yield/young
      eps_33_f_t  = yield/young
      eps_33_f_c  = yield/young
      eps_12_f    = yield/young
      eps_13_f    = yield/young
      eps_23_f    = yield/young
      
      
      sig_11_f =  yield
      sig_22_f =  yield
      sig_33_f =  yield
      
      
      d1   = ustatev(1)
      d2   = ustatev(2)
      d3   = ustatev(3)
	
      
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG


********   LAME'S CONSTANTS ********
      
      lambda   =  young*posn / ((ONE + posn)*(ONE - TWO*posn)) 
      mu       =  young / (TWO*(ONE+posn))
	xK       =  young / (THREE*(ONE - TWO*posn))

      

********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (lambda + TWO * mu)*(1 - d1)**2
	dsdeEl(1,2) = lambda*(1 - d1)*(1 - d2)
	dsdeEl(1,3) = lambda*(1 - d1)*(1 - d3)
	dsdeEl(1,4) = 0
	dsdeEl(1,5) = 0
	dsdeEl(1,6) = 0
	dsdeEl(2,1) = lambda*(1 - d1)*(1 - d2)
	dsdeEl(2,2) = (lambda + TWO * mu)*(1 - d2)**2
	dsdeEl(2,3) = lambda*(1 - d3)*(1 - d2)
	dsdeEl(2,4) = 0
	dsdeEl(2,5) = 0
	dsdeEl(2,6) = 0
	dsdeEl(3,1) = lambda*(1 - d3)*(1 - d1)
	dsdeEl(3,2) = lambda*(1 - d3)*(1 - d2) 
	dsdeEl(3,3) = (lambda + TWO * mu)*(1 - d3)**2
	dsdeEl(3,4) = 0
	dsdeEl(3,5) = 0
	dsdeEl(3,6) = 0
	dsdeEl(4,1) = 0
	dsdeEl(4,2) = 0 
	dsdeEl(4,3) = 0
	dsdeEl(4,4) = mu*(1 - d1)*(1 - d2)
	dsdeEl(4,5) = 0
	dsdeEl(4,6) = 0
	dsdeEl(5,1) = 0
	dsdeEl(5,2) = 0 
	dsdeEl(5,3) = 0
	dsdeEl(5,4) = 0
	dsdeEl(5,5) = mu*(1 - d1)*(1 - d3)
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = mu*(1 - d3)*(1 - d2)

  

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do
      
      
 
c****Damage initiation criteria
 
      if (T_strain(1) > 0) then
          
          
          F_f = SQRT((T_strain(1)/eps_11_f_t)**2 + (T_strain(4) 
     &       /eps_12_f)**2   + (T_strain(5)/eps_13_f)**2)
          
      else 
          
          F_f = SQRT((T_strain(1)/eps_11_f_c)**2)
  
      endif
      
      
      
      if (T_strain(2) + T_strain(3)> 0)  then
          
          F_m =   SQRT(((T_strain(2) + T_strain(3))**2/(eps_22_f_t
     &    *eps_33_f_t)) -  (T_strain(2)*T_strain(3)/eps_23_f**2) +  
     &    (T_strain(4)/eps_12_f)**2 + (T_strain(5)/eps_13_f)**2 +
     &     (T_strain(6)/eps_23_f)**2)  
          
      
      else
          
          F_m =   SQRT(((T_strain(2) + T_strain(3))**2/(eps_22_f_c
     &    *eps_33_f_c)) + ((T_strain(2)+T_strain(3)/eps_22_f_c)  
     &    *((eps_22_f_c/2*eps_12_f)  - 1)) -  
     &    (T_strain(2)*T_strain(3)/eps_23_f**2) +  
     &    (T_strain(4)/eps_12_f)**2 + (T_strain(5)/eps_13_f)**2 +
     &     (T_strain(6)/eps_23_f)**2) 
          
      
      endif
      
      
      
      if (T_strain(3) > 0 )then
      
          F_z = SQRT((T_strain(3)/eps_33_f_t)**2 + (T_strain(5) 
     &   /eps_13_f)**2  + (T_strain(6)/eps_23_f)**2) 
          
      else 
      
          F_z = SQRT((T_strain(3)/eps_33_f_c)**2 + (T_strain(5) 
     &   /eps_13_f)**2  + (T_strain(6)/eps_23_f)**2)
          
          
      endif
      
      
      
      call vzero(sigma, 6)
      call vzero(C_T_1_a, 6) 
      call vzero(C_T_2_a, 6)
      call vzero(C_T_3_a, 6) 
      if ((F_f**2 < 1) .AND. (F_m**2 < 1) .AND. (F_z**2 < 1)) then
       
          do i=1,ncomp
             do j=1,ncomp
                 sigma(i) = sigma(i) + (dsdeEl(i,j)*T_strain(j))        
             end do
          end do  
       
      
      else
      
        k1 =  (-sig_11_f*eps_11_f_t*L_c)/G_c_1
        k2 =  (-sig_22_f*eps_22_f_t*L_c)/G_c_2
        k3 =  (-sig_33_f*eps_33_f_t*L_c)/G_c_3
        
        
        if (F_f**2 >= 1) then
      
          d1_new =  1  - (exp(k2*(F_f - 1)))/F_f
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
   
        endif        

      
        if (F_m**2 >= 1) then
      
          d2_new =  1  - (exp(k2*(F_m - 1)))/F_m
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
   
        endif
      
      
        if (F_z**2 >= 1) then
      
          d3_new =  1  - (exp(k3*(F_z - 1)))/F_z    
      
          if (d3_new > d3) then
              d3 = d3_new
          else
              d3 = d3
          endif

        endif 
      
      
        dsdeEl_d(1,1) = (lambda + TWO * mu)*(1 - d1)**2
        dsdeEl_d(1,2) = lambda*(1 - d1)*(1 - d2)
        dsdeEl_d(1,3) = lambda*(1 - d1)*(1 - d3)
        dsdeEl_d(1,4) = 0
        dsdeEl_d(1,5) = 0
        dsdeEl_d(1,6) = 0
        dsdeEl_d(2,1) = lambda*(1 - d1)*(1 - d2)
        dsdeEl_d(2,2) = (lambda + TWO * mu)*(1 - d2)**2
        dsdeEl_d(2,3) = lambda*(1 - d3)*(1 - d2)
        dsdeEl_d(2,4) = 0
        dsdeEl_d(2,5) = 0
        dsdeEl_d(2,6) = 0
        dsdeEl_d(3,1) = lambda*(1 - d3)*(1 - d1)
        dsdeEl_d(3,2) = lambda*(1 - d3)*(1 - d2) 
        dsdeEl_d(3,3) = (lambda + TWO * mu)*(1 - d3)**2
        dsdeEl_d(3,4) = 0
        dsdeEl_d(3,5) = 0
        dsdeEl_d(3,6) = 0
        dsdeEl_d(4,1) = 0
        dsdeEl_d(4,2) = 0 
        dsdeEl_d(4,3) = 0
        dsdeEl_d(4,4) = mu*(1 - d1)*(1 - d2)
        dsdeEl_d(4,5) = 0
        dsdeEl_d(4,6) = 0
        dsdeEl_d(5,1) = 0
        dsdeEl_d(5,2) = 0 
        dsdeEl_d(5,3) = 0
        dsdeEl_d(5,4) = 0
        dsdeEl_d(5,5) = mu*(1 - d1)*(1 - d3)
        dsdeEl_d(5,6) = 0
        dsdeEl_d(6,1) = 0
        dsdeEl_d(6,2) = 0 
        dsdeEl_d(6,3) = 0
        dsdeEl_d(6,4) = 0
        dsdeEl_d(6,5) = 0
        dsdeEl_d(6,6) = mu*(1 - d3)*(1 - d2)
         
        
        if (d1 == 0) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_1(i,j) = 0        
             end do
          end do  
          
        else
            
          do i=1,mcomp
             do j=1,mcomp
                 d_C_d_d1(i,j) = 0        
             end do
          end do 
          
          d_C_d_d1(1,1) =  -TWO*(lambda + TWO*mu)*(1 - d1)    
          d_C_d_d1(1,2) = lambda*(d2 - 1)
          d_C_d_d1(1,3) = lambda*(d3 - 1)
          d_C_d_d1(2,1) = lambda*(d2 - 1)
          d_C_d_d1(3,1) = lambda*(d3 - 1)
          d_C_d_d1(4,4) = mu*(d2 -1)
          d_C_d_d1(5,5) = mu*(d3 -1)
          
          
          
          do i=1,mcomp
              do j=1,mcomp
                  C_T_1_a(i) = C_T_1_a(i) + (d_C_d_d1(i,j)*T_strain(j)) 
              end do
          end do 
          
          
          if (T_strain(1)> 0)  then
          
            C_T_1_b  = (/ ((1 - k1*F_f)/(F_f**3 * eps_11_f_t**2))
     &     *exp(k1*(F_f - 1))*T_strain(1), ZERO, ZERO, ((1 - k1*F_f)/ 
     &     (F_f**3 * eps_12_f**2))*exp(k1*(F_f - 1))*T_strain(4), ZERO,     
     &     ((1 - k1*F_f)/(F_f**3 * eps_13_f**2))*exp(k1*(F_f - 1))
     &     *T_strain(6)/)
            
          
          else
              
            C_T_1_b  = (/ ((1 - k1*F_f)/(F_f**3 * eps_11_f_c**2))
     &     *exp(k1*(F_f - 1))*T_strain(1),ZERO,ZERO,ZERO,ZERO,ZERO/)
            
            
          endif
      
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_1(i,j) = C_T_1_a(i)*C_T_1_b(j)
              end do
          end do
       
        endif
            

        if (d2 == 0) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_2(i,j) = 0        
             end do
          end do 
      
        else 
            
          do i=1,mcomp
             do j=1,mcomp
                 d_C_d_d2(i,j) = 0        
             end do
          end do 
          
          d_C_d_d2(1,2) = lambda*(d1 - 1)
          d_C_d_d2(2,1) = lambda*(d1 - 1)
          d_C_d_d2(2,2) = -TWO*(lambda + TWO*mu)*(1 - d2)
          d_C_d_d2(2,3) = lambda*(d3 - 1)
          d_C_d_d2(3,2) = lambda*(d3 - 1)
          d_C_d_d2(4,4) = mu*(d1 - 1)
          d_C_d_d2(6,6) = mu*(d3 - 1)
          
          
          do i=1,ncomp
              do j=1,ncomp
                  C_T_2_a(i) = C_T_2_a(i) + (d_C_d_d2(i,j)*T_strain(j)) 
              end do
          end do 
          
          
          if (T_strain(2) + T_strain(3)> 0)  then
          
              term1 = (TWO*(T_strain(2) + T_strain(3))/(eps_22_f_t
     &        *eps_33_f_t))  -  (T_strain(3)/eps_23_f**2)
          
              term2  =  (TWO*(T_strain(2) + T_strain(3))/(eps_22_f_t
     &        *eps_33_f_t))  -  (T_strain(2)/eps_23_f**2)
              
              C_T_2_b = (/ZERO, (HALF)*((1 - k2*F_m)/(F_m**3))*
     &        exp(k2*(F_m - 1))*term1, (HALF)*((1 - k2*F_m)/(F_m**3))
     &        *exp(k2*(F_m - 1))*term2, ((1 - k2*F_m)/(F_m**3 *     
     &        eps_12_f**2))*exp(k2*(F_m - 1))*T_strain(4), ((1 - k2*F_m) 
     &        /(F_m**3 * eps_23_f**2))*exp(k2*(F_m - 1))*T_strain(5),  
     &        ((1 - k2*F_m)/(F_m**3 * eps_13_f**2))*exp(k2*(F_m - 1))         
     &        *T_strain(6) /)         
     

      
          else 
                            
              term1 = ((TWO*(T_strain(2) + T_strain(3))/(eps_22_f_c
     &    *eps_33_f_c))  +  ((eps_22_f_c/(2*eps_12_f)) - 1)/eps_22_f_c  
     &     -  (T_strain(3)/eps_23_f**2))/(2*F_m)
              
              term2 = ((TWO*(T_strain(2) + T_strain(3))/(eps_22_f_c
     &    *eps_33_f_c))  +  ((eps_22_f_c/(2*eps_12_f)) - 1)/eps_22_f_c  
     &     -  (T_strain(2)/eps_23_f**2))/(2*F_m) 
              
              
              C_T_2_b  = (/ZERO, ((1 - k2*F_m)/(F_m**2))*
     &        exp(k2*(F_m - 1))*term1, ((1 - k2*F_m)/(F_m**2))
     &        *exp(k2*(F_m - 1))*term2,  ((1 - k2*F_m)/(F_m**3 *     
     &        eps_12_f**2))*exp(k2*(F_m - 1))*T_strain(4), ((1 - k2*F_m) 
     &        /(F_m**3 * eps_23_f**2))*exp(k2*(F_m - 1))*T_strain(5),  
     &        ((1 - k2*F_m)/(F_m**3 * eps_13_f**2))*exp(k2*(F_m - 1))         
     &        *T_strain(6) /)
              
              
          endif

              
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_2(i,j) = C_T_2_a(i)*C_T_2_b(j)
              end do
          end do              
              
        endif      
      
      
        if (d3 == 0) then
      
              do i=1,ncomp
                 do j=1,ncomp
                     C_T_3(i,j) = 0        
                 end do
              end do 
          
        else 
            
          do i=1,mcomp
             do j=1,mcomp
                 d_C_d_d3(i,j) = 0        
             end do
          end do 
          
         d_C_d_d3(1,3)  =  lambda*(d1 - 1) 
         d_C_d_d3(2,3)  =  lambda*(d2 - 1)     
         d_C_d_d3(3,1)  =  lambda*(d1 - 1)     
         d_C_d_d3(3,2)  =  lambda*(d2 - 1)     
         d_C_d_d3(3,3)  =  -TWO*(lambda + TWO*mu)*(1 - d3)     
         d_C_d_d3(5,5)  =  mu*(d1 - 1)     
         d_C_d_d3(6,6)  =  mu*(d2 - 1)   
         
         
          
          do i=1,ncomp
              do j=1,ncomp
                  C_T_3_a(i) = C_T_3_a(i) + (d_C_d_d3(i,j)*T_strain(j)) 
              end do
          end do 
          
          
          if (T_strain(3) > 0) then
              
              C_T_3_b =(/ZERO,ZERO,((1 - k3*F_z)/(F_z**3*eps_33_f_t**2))
     &    *exp(k3*(F_z - 1))*T_strain(3),ZERO,((1 - k3*F_z)/(F_z**3       
     &    * eps_23_f**2))*exp(k3*(F_z - 1))*T_strain(5), ((1 - k3*F_z)         
     &    /(F_z**3 * eps_13_f**2))*exp(k3*(F_z - 1))*T_strain(6) /)   
              
              
          else 
              
              C_T_3_b =(/ZERO,ZERO,((1 - k3*F_z)/(F_z**3*eps_33_f_c**2))
     &    *exp(k3*(F_z - 1))*T_strain(3),ZERO,((1 - k3*F_z)/(F_z**3       
     &    * eps_23_f**2))*exp(k3*(F_z - 1))*T_strain(5), ((1 - k3*F_z)         
     &    /(F_z**3 * eps_13_f**2))*exp(k3*(F_z - 1))*T_strain(6) /)
              
          
          endif
          
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_3(i,j) = C_T_3_a(i)*C_T_3_b(j)
              end do
          end do 
          
        endif
              
              
        do i=1,ncomp
          do j=1,ncomp
              sigma(i) = sigma(i) + (dsdeEl_d(i,j)*T_strain(j))        
          end do
        end do 
  
      
        do i=1,ncomp
           do j=1,ncomp
              C_T(i,j)   =  dsdeEl_d(i,j) +  C_T_1(i,j) + C_T_2(i,j) +
     &        C_T_3(i,j) 
           end do
        end do
      
      
      endif
      

      ustatev(1) = d1
      ustatev(2) = d2
      ustatev(3) = d3
      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f**2 < 1) .AND. (F_m**2 < 1) .AND. (F_z**2 < 1)) then
          
          do i=1,ncomp
              do j=1,ncomp
                  dsdePl(i,j) = dsdeEl(i,j)        
              end do
          end do 
      
      else
          
          do i=1,ncomp
              do j=1,ncomp
                  dsdePl(i,j) = C_T(i,j)        
              end do
          end do 
          
      endif  
      
      
      
      return
      end



              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      