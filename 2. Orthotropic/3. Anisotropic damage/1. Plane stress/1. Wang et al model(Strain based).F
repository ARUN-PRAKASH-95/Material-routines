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
	
	
      EXTERNAL         usermatps
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
	
	

	if (ncomp == 3) then
	
		
		call usermatps (
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
      subroutine usermatps(
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
     &                 mcomp      = 3
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
     &                  C_T_3_b(mcomp),value(mcomp),arr(mcomp)
                  

      
      integer             i,j,k,l,n,P,ielem,iiter,key
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz, L_c,k1,k2,k3,
     &                    xK, yield,G_c_1,G_c_2,G_c_3,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,d1_new,d2_new
      
      double precision    d3_new,F_f,F_m,F_z,term,term2,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,V,eta,D

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

 
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	pr_xy    = prop(3)
	g_xy     = prop(4)
      pr_yx = (young_y * pr_xy) / young_x
      
      
      sig_11_f_t = 2500e6
      sig_11_f_c = -2000e6
      sig_22_f_t = 50e6
      sig_22_f_c = -150e6
      sig_12_f   = 50e6
      G_c_1      = 12.5e3
      G_c_2      = 1e3
      G_c_3      = 1e3
      L_c        = 1**(1/3)
    
   
      
      d1   = ustatev(1)
      d2   = ustatev(2)
      d3   = ONE - ((ONE - d1)*(ONE - d2))
      F_f  = ustatev(3)
      F_m  = ustatev(4)
	

********* ELASTIC STIFFNESS MATRIX *********
	D  =  ONE  -  ((ONE - d1)*(ONE - d2)*pr_xy*pr_yx)
	dsdeEl(1,1) = young_x*(ONE - d1)*(ONE/D)
	dsdeEl(1,2) = pr_xy * young_y *(ONE - d1)*(ONE - d2)*(ONE/D)
	dsdeEl(1,3) = ZERO
	dsdeEl(2,1) = pr_xy * young_y *(ONE - d1)*(ONE - d2)*(ONE/D)
      dsdeEl(2,2) = young_y*(ONE - d2)*(ONE/D)
	dsdeEl(2,3) = ZERO
	dsdeEl(3,1) = ZERO
      dsdeEl(3,2) = ZERO
      dsdeEl(3,3) = g_xy*(ONE-d3)


  
      
      eps_11_f_t = sig_11_f_t / (young_x/(ONE - pr_xy*pr_yx))
      eps_11_f_c = sig_11_f_c / (young_x/(ONE - pr_xy*pr_yx))
      eps_22_f_t = sig_22_f_t / (young_y/(ONE - pr_xy*pr_yx))
      eps_22_f_c = sig_22_f_c / (young_y/(ONE - pr_xy*pr_yx))
      eps_12_f   = sig_12_f / g_xy

      
      
c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do

 
c****Damage initiation criteria
 
      if (T_strain(1) >= 0) then
          
          
          F_f_new = SQRT((T_strain(1)/eps_11_f_t)**2 + (T_strain(3) 
     &       /eps_12_f)**2   )
          
      else 
          
          F_f_new = SQRT((T_strain(1)/eps_11_f_c)**2)
  
      endif
      
      
      
      if (T_strain(2) >= 0)  then
          
          F_m_new =   SQRT( (T_strain(2)/eps_22_f_t)**2 +
     &    (T_strain(3)/eps_12_f)**2 ) 
          
      
      else
          
          F_m_new =   SQRT( (T_strain(2)/eps_22_f_t)**2
     &    + ((T_strain(2)/eps_22_f_c)  
     &    *((eps_22_f_c/2*eps_12_f)  - ONE)) +  
     &    (T_strain(3)/eps_12_f)**2 ) 
          
      
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

      
      
      call vzero(sigma, 3)
      call vzero(C_T_1_a, 3) 
      call vzero(C_T_2_a, 3)
      call vzero(C_T_3_a, 3) 
      if ((F_f**2 <= 1) .AND. (F_m**2 <= 1) ) then
       
c*********    Compute stress using Hooke's Law    **********
          do i=1,ncomp
             do j=1,ncomp
                 sigma(i) = sigma(i) + (dsdeEl(i,j)*T_strain(j))        
             end do
          end do  
       
      
      else
          
c*********  Terms in damage evolution equations   ***********
        if (T_strain(1) >= 0) then
          k1 =   (-sig_11_f_t*eps_11_f_t*L_c)/G_c_1
        else 
          k1 =   (-sig_11_f_c*eps_11_f_c*L_c)/G_c_1
        endif
        
        
        if (T_strain(2) >= 0) then
          k2 =  (-sig_22_f_t*eps_22_f_t*L_c)/G_c_2
        else 
          k2 =  (-sig_22_f_c*eps_22_f_c*L_c)/G_c_2
        endif
      

c*********   Damage evolution equations    *********** 	  
        if (F_f**2 > 1) then
      
          d1_new =  ONE  - (exp(k1*(F_f - ONE)))/F_f
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
          
        endif        

      
        if (F_m**2 > 1) then
      
          d2_new =  ONE  - (exp(k2*(F_m - ONE)))/F_m
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
        
        endif
      
      
        d3  = 1 - ((1 -d1)*(1 - d2))
        
        
c********      Degraded stiffness   **********      
        D  =  ONE  -  ((ONE - d1)*(ONE - d2)*pr_xy*pr_yx)
        dsdeEl_d(1,1) = young_x*(ONE - d1)*(ONE/D)
        dsdeEl_d(1,2) = pr_xy * young_y *(ONE - d1)*(ONE - d2)*(ONE/D)
        dsdeEl_d(1,3) = ZERO
        dsdeEl_d(2,1) = pr_xy * young_y *(ONE - d1)*(ONE - d2)*(ONE/D)
        dsdeEl_d(2,2) = young_y*(ONE - d2)*(ONE/D)
        dsdeEl_d(2,3) = ZERO
        dsdeEl_d(3,1) = ZERO
        dsdeEl_d(3,2) = ZERO
        dsdeEl_d(3,3) = g_xy*(ONE-d3)

         
        
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
          d_C_d_d1(1,1) = -young_x/D**2
          d_C_d_d1(1,2) = ((d2 - ONE)*pr_xy*young_y)/(D**2)
          d_C_d_d1(2,1) = ((d2 - ONE)*pr_xy*young_y)/(D**2)
          d_C_d_d1(2,2) = ((d2 - ONE)**2*pr_xy*pr_yx*young_y)/(D**2)            d_C_d_d1(3,3) =  g_xy*(1 - d2)
          
          
c********   (d_C_d/d1 : eps)   **********          
          do i=1,mcomp
              do j=1,mcomp
                  C_T_1_a(i) = C_T_1_a(i) + (d_C_d_d1(i,j)*T_strain(j)) 
              end do
          end do 
          
c**********   Derivative of d1 with respect to strain (d_d1/d_epsilon)   ***********
 
c**       For Tension            
          if (T_strain(1)> 0)  then
          
            C_T_1_b  = (/ ((ONE - k1*F_f)/(F_f**3 * eps_11_f_t**2))
     &     *exp(k1*(F_f - ONE))*T_strain(1), ZERO, ((ONE - k1*F_f)/ 
     &     (F_f**3 * eps_12_f**2))*exp(k1*(F_f -ONE))*T_strain(3)/)
            
c**       For Compression           
          else
              
            C_T_1_b  = (/ ((ONE - k1*F_f)/(F_f**2 * eps_11_f_c))
     &     *exp(k1*(F_f - ONE)),ZERO,ZERO/)
            
            
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
          d_C_d_d2(1,1) = ((d1 - ONE)**2*pr_xy*pr_yx*young_x)/(D**2) 
          d_C_d_d2(1,2) = ((d1 - ONE)*pr_xy*young_y)/(D**2)
          d_C_d_d2(2,1) = ((d1 - ONE)*pr_xy*young_y)/(D**2)
          d_C_d_d2(2,2) = -young_y/D**2
          d_C_d_d2(3,3) =  g_xy*(1 - d1)
          
c*****   (d_C_d/d2 : eps)         
          do i=1,ncomp
              do j=1,ncomp
                  C_T_2_a(i) = C_T_2_a(i) + (d_C_d_d2(i,j)*T_strain(j)) 
              end do
          end do 
          
c******     Derivative of d2 with respect to strain (d_d2/d_epsilon)
      
c****      For Tension           
          if (T_strain(2) > 0)  then
          
              
              C_T_2_b=(/ZERO, ((ONE-k2*F_m)/(F_m**3))*
     &        exp(k2*(F_m -ONE))*T_strain(2),  ((ONE - k2*F_m)/
     &        (F_m**3 * eps_12_f**2))*
     &        exp(k2*(F_m - ONE))*T_strain(3) /)         
     

c**       For Compression                       
          else 
                            
              term = ( TWO*T_strain(2) )/(eps_22_f_c**2)
     &      + ((eps_22_f_c/(2*eps_12_f))- ONE)/eps_22_f_c  
     
              
 
              
              C_T_2_b  = (/ZERO, ((ONE - k2*F_m)/(F_m**2))*
     &        exp(k2*(F_m - ONE))*term, ((ONE - k2*F_m)
     &        /(F_m**3 * eps_12_f**2))*
     &        exp(k2*(F_m - ONE))*T_strain(3) /)
              
              
          endif

              
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_2(i,j) = C_T_2_a(i)*C_T_2_b(j)
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
              C_T(i,j)   =  dsdeEl_d(i,j) +  C_T_1(i,j) + C_T_2(i,j) 
           end do
        end do
      
      
      endif
      

      ustatev(1) = d1
      ustatev(2) = d2
      ustatev(3) = F_f
      ustatev(4) = F_m
  
      

      
 
      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ( (F_f**2 <= 1) .AND. (F_m**2 <= 1) ) then
          
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



              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      