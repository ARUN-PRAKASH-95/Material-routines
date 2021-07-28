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
	
      integer :: mcomp, twenty
      double precision :: HALF, ONE, TWO, ZERO,
     &                    THREE
      
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 mcomp      = 6
     &                        )
	
      EXTERNAL          vzero, vmove
   
      double precision  dsdeEl(mcomp,mcomp), sigma(mcomp),
     &                  sigi (mcomp), C_T(mcomp,mcomp),
     &                  dsdeEl_d(mcomp,mcomp), T_strain(mcomp),
     &                  id(3,3), C_T_1(mcomp,mcomp),
     &                  d_C_d_d1(mcomp,mcomp), C_T_1_a(mcomp), 
     &                  C_T_1_b(mcomp), C_T_2(mcomp,mcomp), 
     &                  d_C_d_d2(mcomp,mcomp), C_T_2_a(mcomp),
     &                  C_T_2_b(mcomp),C_T_3(mcomp,mcomp), 
     &                  d_C_d_d3(mcomp,mcomp), C_T_3_a(mcomp),
     &                  C_T_3_b(mcomp),arr(mcomp),
     &                  a(3),b(3)
                  

      
      integer             i,j,k,l,n,P1,P2,P3
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz,epsilon_f, 
     &                    epsilon_m,epsilon_z,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,d1_new,d2_new
      
      double precision    d3_new,F_f,F_m,F_z,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,V,key,ans,V2,value,L_c

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

 


      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	young_z  = prop(3)
	pr_xy    = prop(4)
	pr_yz    = prop(5)
	pr_xz    = prop(6)
	g_xy     = 5500e6
	g_yz     = 3000e6
	g_xz     = 5500e6
      
      
      sig_11_f_t = 2500e6
      sig_11_f_c = -2000e6
      sig_22_f_t = 50e6
      sig_22_f_c = -200e6
      sig_33_f_t = 50e6
      sig_33_f_c = -200e6
      sig_12_f   = 50e6
      sig_13_f   = 50e6
      sig_23_f   = 50e6

    
   
      
	pr_yx = (young_y * pr_xy) / young_x
	pr_zy = (young_z * pr_yz) / young_y
	pr_zx = (young_z * pr_xz) / young_x
	
	
	xy_yx = pr_xy*pr_yx
	yz_zy = pr_yz*pr_zy 
	zx_xz = pr_zx*pr_xz
	xyz   = TWO*pr_xy*pr_yz*pr_zx
	E_xyz = young_x*young_y*young_z
	
	delta = (ONE - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz

      V    = ustatev(1)
      d1   = ustatev(2)
      d2   = ustatev(3)
      d3   = ustatev(4)
      F_f  = ustatev(5)
      F_m  = ustatev(6)
      F_z  = ustatev(7)
	
      L_c  = V**(1/3)

      tsstif(1) = g_xz
      tsstif(2) = g_yz  
     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (ONE -yz_zy) / (young_y*young_z*delta)*(ONE - d1)**2
	dsdeEl(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)
     &              *(ONE - d1)*(ONE - d2)
	dsdeEl(1,3) = (pr_zx + pr_yx*pr_zy)/(young_y*young_z*delta)*(ONE 
     &              - d1)*(ONE - d3)
      dsdeEl(1,4) = 0
	dsdeEl(1,5) = 0
	dsdeEl(1,6) = 0
	dsdeEl(2,1) = (pr_yx + pr_zx*pr_yz) /(young_y*young_z*delta)*(ONE
     &              -d1)*(ONE - d2)
      dsdeEl(2,2) = (ONE -zx_xz)/(young_x*young_z*delta)*(ONE-d2)**2
	dsdeEl(2,3) = (pr_zy + pr_zx*pr_xy) /(young_x*young_z*delta)*(ONE 
     &             - d3)*(ONE  - d2)
      dsdeEl(2,4) = 0
	dsdeEl(2,5) = 0
	dsdeEl(2,6) = 0
	dsdeEl(3,1) = (pr_zx + pr_yx*pr_zy) /(young_y*young_z*delta)*(ONE  
     &              - d3)*(ONE - d1)
      dsdeEl(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta)*
     &              (ONE - d3)*(ONE - d2)
      dsdeEl(3,3) = (ONE -xy_yx)/(young_x*young_y*delta)*(ONE-d3)**2
	dsdeEl(3,4) = 0
	dsdeEl(3,5) = 0
	dsdeEl(3,6) = 0
	dsdeEl(4,1) = 0
	dsdeEl(4,2) = 0 
	dsdeEl(4,3) = 0
	dsdeEl(4,4) = g_xy*(ONE - d1)*(ONE - d2)
	dsdeEl(4,5) = 0
	dsdeEl(4,6) = 0
	dsdeEl(5,1) = 0
	dsdeEl(5,2) = 0 
	dsdeEl(5,3) = 0
	dsdeEl(5,4) = 0
	dsdeEl(5,5) = g_yz*(ONE - d2)*(ONE - d3)
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = g_xz*(ONE - d3)*(ONE - d1)

  
      
      eps_11_f_t = sig_11_f_t / ((ONE -yz_zy) / (young_y*young_z*delta))
      eps_11_f_c = sig_11_f_c / ((ONE -yz_zy) / (young_y*young_z*delta))
      eps_22_f_t = sig_22_f_t / ((ONE -zx_xz) / (young_x*young_z*delta))
      eps_22_f_c = sig_22_f_c / ((ONE -zx_xz) / (young_x*young_z*delta))
      eps_33_f_t = sig_33_f_t / ((ONE -xy_yx) / (young_x*young_y*delta))
      eps_33_f_c = sig_33_f_c / ((ONE -xy_yx) / (young_x*young_y*delta))
      eps_12_f   = sig_12_f / g_xy
      eps_13_f   = sig_13_f / g_xz 
      eps_23_f   = sig_23_f / g_yz

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do

 
c****Damage initiation criteria
 
      if (T_strain(1) >= 0) then          
          F_f_new = T_strain(1)/eps_11_f_t          
      else           
          F_f_new = T_strain(1)/eps_11_f_c  
      endif
      
            
      if (T_strain(2)>= 0)  then          
          F_m_new = T_strain(2)/eps_22_f_t                 
      else          
          F_m_new = T_strain(2)/eps_22_f_c                 
      endif
     
      
      if (T_strain(3) >= 0 )then      
          F_z_new = T_strain(3)/eps_33_f_t          
      else       
          F_z_new = T_strain(3)/eps_33_f_c                    
      endif
      

c******    To make sure damage initiation criteria is greater than or equal to previous step   *******   
      
      if (F_f_new >= F_f) then          
          F_f = F_f_new          
      else          
          F_f = F_f          
      endif
      
      if (F_m_new >= F_m) then          
          F_m = F_m_new          
      else          
          F_m = F_m          
      endif
          
      if (F_z_new >= F_z) then          
          F_z = F_z_new         
      else          
          F_z = F_z          
      endif

      
      
      call vzero(sigma, 6)
      call vzero(C_T_1_a, 6) 
      call vzero(C_T_2_a, 6)
      call vzero(C_T_3_a, 6) 
      if ((F_f <= 1) .AND. (F_m <= 1) .AND. (F_z <= 1)) then

c*********    Compute stress using Hooke's Law    **********
          do i=1,ncomp
             do j=1,ncomp
                 sigma(i) = sigma(i) + (dsdeEl(i,j)*T_strain(j))        
             end do
          end do  
       
      
      else
 
c*********  Terms in damage evolution equations   ***********
        if (T_strain(1) >= 0) then
          epsilon_f = eps_11_f_t
          P1 = 25
        else 
          epsilon_f = eps_11_f_c
          P1 = -25
        endif
        
        
        if (T_strain(2) >= 0) then
          epsilon_m = eps_22_f_t
          P2 = 25
        else 
          epsilon_m = eps_22_f_c
          P2 = -25
        endif
      
        
        if (T_strain(3) >= 0) then
          epsilon_z =  eps_33_f_t
          P3 = 25
        else 
          epsilon_z =  eps_33_f_c
          P3 = -25
        endif
       
        
	  
c*********   Damage evolution equations    *********** 	  
        if (F_f > 1) then
      
          d1_new =  ONE -  (exp(-P1*(T_strain(1) - epsilon_f)))
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif

        endif        

      
        if (F_m > 1) then
      
          d2_new =  ONE -  (exp(-P2*(T_strain(2) - epsilon_m)))
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif

        endif
      
      
        if (F_z > 1) then
      
          d3_new = ONE -  (exp(-P3*(T_strain(3) - epsilon_z)))
      
          if (d3_new > d3) then
              d3 = d3_new
          else
              d3 = d3
          endif
 
        endif 
      
c********      Degraded stiffness   **********      
        dsdeEl_d(1,1)=(ONE-yz_zy)/(young_y*young_z*delta)*(ONE
     &                -d1)**2
        dsdeEl_d(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)
     &              *(ONE - d1)*(ONE - d2)
        dsdeEl_d(1,3) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta)*
     &              (ONE - d1)*(ONE - d3)
        dsdeEl_d(1,4) = 0
        dsdeEl_d(1,5) = 0
        dsdeEl_d(1,6) = 0
        dsdeEl_d(2,1) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)*
     &              (ONE - d1)*(ONE - d2)
        dsdeEl_d(2,2)=(ONE -zx_xz)/(young_x*young_z*delta)*(ONE - d2)**2
        dsdeEl_d(2,3) = (pr_zy + pr_zx*pr_xy)/(young_x*young_z*delta)*
     &              (ONE - d3)*(ONE - d2)
        dsdeEl_d(2,4) = 0
        dsdeEl_d(2,5) = 0
        dsdeEl_d(2,6) = 0
        dsdeEl_d(3,1) = (pr_zx + pr_yx*pr_zy) /(young_y*young_z*delta)*
     &              (ONE - d3)*(ONE - d1)
        dsdeEl_d(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta)*
     &              (ONE - d3)*(ONE- d2)
        dsdeEl_d(3,3) = (ONE -xy_yx)/(young_x*young_y*delta)*(ONE
     &                -d3)**2
        dsdeEl_d(3,4) = 0
        dsdeEl_d(3,5) = 0
        dsdeEl_d(3,6) = 0
        dsdeEl_d(4,1) = 0
        dsdeEl_d(4,2) = 0 
        dsdeEl_d(4,3) = 0
        dsdeEl_d(4,4) = g_xy*(ONE - d1)*(ONE - d2)
        dsdeEl_d(4,5) = 0
        dsdeEl_d(4,6) = 0
        dsdeEl_d(5,1) = 0
        dsdeEl_d(5,2) = 0 
        dsdeEl_d(5,3) = 0
        dsdeEl_d(5,4) = 0
        dsdeEl_d(5,5) = g_yz*(ONE - d2)*(ONE - d3)
        dsdeEl_d(5,6) = 0
        dsdeEl_d(6,1) = 0
        dsdeEl_d(6,2) = 0 
        dsdeEl_d(6,3) = 0
        dsdeEl_d(6,4) = 0
        dsdeEl_d(6,5) = 0
        dsdeEl_d(6,6) = g_xz*(ONE - d3)*(ONE - d1)
         
        
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

c*********     First term C_T_1 ((d_C_d/d1 : eps) outerProduct (d_d1/d_epsilon))  ************

          d_C_d_d1(1,1) = -TWO*((ONE -yz_zy) / (young_y*young_z*delta))*
     &                    (ONE - d1) 
          d_C_d_d1(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*
     &                     delta))*(d2 - ONE)     
          d_C_d_d1(1,3) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*
     &                     delta))*(d3 - ONE)
          d_C_d_d1(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*
     &                     delta))*(d2 - ONE)
          d_C_d_d1(3,1) = ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*
     &                     delta))*(d3 - ONE)
          d_C_d_d1(4,4) = g_xy*(d2 -ONE)
          d_C_d_d1(6,6) = g_xz*(d3 -ONE)
          
          
c********   (d_C_d/d1 : eps)   **********          
          do i=1,mcomp
              do j=1,mcomp
                  C_T_1_a(i) = C_T_1_a(i) + (d_C_d_d1(i,j)*T_strain(j)) 
              end do
          end do 
          
 
c**********   Derivative of d1 with respect to strain (d_d1/d_epsilon)   ***********
 
c**       For Tension 
          if (T_strain(1)> 0)  then
          
            C_T_1_b  = (/ P1*exp(-P1*(T_strain(1) - eps_11_f_t)), ZERO,
     &     ZERO, ZERO, ZERO, ZERO/) 
    
            
c**       For Compression           
          else
              
            C_T_1_b  = (/ P1*exp(-P1*(T_strain(1) - eps_11_f_c)), ZERO,
     &     ZERO, ZERO, ZERO, ZERO/) 
            
            
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
          
c*******   Second term C_T_2 ((d_C_d/d2 : eps) outerProduct (d_d2/d_epsilon))    *********       
          d_C_d_d2(1,2) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*
     &                    delta))*(d1 - ONE)
          d_C_d_d2(2,1) = ((pr_yx + pr_zx*pr_yz) / (young_y*young_z*
     &                    delta))*(d1 - ONE)
          d_C_d_d2(2,2) = -TWO*((1 -zx_xz) / (young_x*young_z*
     &                    delta))*(ONE - d2)
          d_C_d_d2(2,3) =  ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*
     &                     delta))*(d3 - ONE)
          d_C_d_d2(3,2) = ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*
     &                     delta))*(d3 - ONE)
          d_C_d_d2(4,4) = g_xy*(d1 - ONE)
          d_C_d_d2(5,5) = g_yz*(d3 - ONE)
          
c*****   (d_C_d/d2 : eps)         
          do i=1,ncomp
              do j=1,ncomp
                  C_T_2_a(i) = C_T_2_a(i) + (d_C_d_d2(i,j)*T_strain(j)) 
              end do
          end do 
 
          
c******     Derivative of d2 with respect to strain (d_d2/d_epsilon)
      
c****      For Tension       
          if (T_strain(2)>= 0)  then
          
            C_T_2_b  = (/ ZERO, P2*exp(-P2*(T_strain(2) - eps_22_f_t)),
     &     ZERO, ZERO, ZERO, ZERO/)            
c**       For Compression                 
          else 
                            
            C_T_2_b  = (/ ZERO, P2*exp(-P2*(T_strain(2) - eps_22_f_c)),
     &     ZERO, ZERO, ZERO, ZERO/) 
           
                            
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

c*****  Third term C_T_3 ((d_C_d/d3 : eps) outerProduct (d_d3/d_epsilon)) 
         d_C_d_d3(1,3)  =  ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*
     &                     delta))*(d1 - ONE)
         d_C_d_d3(2,3)  =  ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*
     &                     delta))*(d2 - ONE) 
         d_C_d_d3(3,1)  =  ((pr_zx + pr_yx*pr_zy) / (young_y*young_z*
     &                     delta))*(d1 - ONE) 
         d_C_d_d3(3,2)  =  ((pr_zy + pr_zx*pr_xy) / (young_x*young_z*
     &                     delta))*(d2 - ONE) 
         d_C_d_d3(3,3)  =  -TWO*((1 -xy_yx) / (young_x*young_y*
     &                     delta))*(ONE - d3) 
         d_C_d_d3(5,5)  =  g_yz*(d2 - ONE)     
         d_C_d_d3(6,6)  =  g_xz*(d1 - ONE)  
         
         
c*****      (d_C_d/d3 : eps)        
          do i=1,ncomp
              do j=1,ncomp
                  C_T_3_a(i) = C_T_3_a(i) + (d_C_d_d3(i,j)*T_strain(j)) 
              end do
          end do 
          
c*****     Derivative of d3 with respect to strain (d_d3/d_epsilon)
          
c****     For Tension 
          if (T_strain(3) > 0) then
              
            C_T_3_b  = (/ ZERO, ZERO,
     &      P3*exp(-P3*(T_strain(3) - eps_33_f_t)), ZERO, ZERO, ZERO/) 
              
c****     For Compression              
          else 
              
            C_T_3_b  = (/ ZERO, ZERO,
     &      P3*exp(-P3*(T_strain(3) - eps_33_f_c)), ZERO, ZERO, ZERO/) 
              
          
          endif
          
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_3(i,j) = C_T_3_a(i)*C_T_3_b(j)
              end do
          end do 
          
        endif
              
c********    Compute stress              
        do i=1,ncomp
          do j=1,ncomp
              sigma(i) = sigma(i) + (dsdeEl_d(i,j)*T_strain(j))        
          end do
        end do 
  
c*********   Tangent stiffness     
        do i=1,ncomp
           do j=1,ncomp
              C_T(i,j)   =  dsdeEl_d(i,j) +  C_T_1(i,j) + C_T_2(i,j) +
     &        C_T_3(i,j) 
           end do
        end do
      
      
      endif
      

      ustatev(2) = d1
      ustatev(3) = d2
      ustatev(4) = d3
      ustatev(5) = F_f
      ustatev(6) = F_m
      ustatev(7) = F_z
      

      
 
      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f <= 1) .AND. (F_m <= 1) .AND. (F_z <= 1)) then
          
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


      