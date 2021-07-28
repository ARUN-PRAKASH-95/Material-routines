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
     &                 mcomp      = 3
     &                        )
	
      EXTERNAL          vzero, vmove
   
      double precision  dsdeEl(mcomp,mcomp), sigma(mcomp),
     &                  sigi (mcomp), C_T(mcomp,mcomp),
     &                  dsdeEl_d(mcomp,mcomp), T_strain(mcomp),
     &                  a(mcomp,mcomp), b(mcomp), d(mcomp,mcomp),
     &                  c(mcomp),M(mcomp,mcomp),sigma_eff(mcomp)
                  

      
      integer             i,j,k,l,n,P1,P2,P3
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz,epsilon_f, 
     &                    epsilon_m,epsilon_z, damage,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,damage_new
      
      double precision    d1_new,d2_new,d3_new,F_f,F_m,F_z,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,V,key,ans,V2,value

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

 


      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	pr_xy    = prop(3)
	g_xy     = prop(4)
      sig_11_f_t  = prop(5)
      sig_11_f_c  = prop(6)
      sig_22_f_t  = 50e6
      sig_22_f_c  = -150e6
      pr_yx = (young_y * pr_xy) / young_x


      damage = ustatev(1)
      d1     = ustatev(2)
      d2     = ustatev(3)
      F_f    = ustatev(4)
      F_m    = ustatev(5)

     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = young_x/(ONE - pr_xy*pr_yx)
	dsdeEl(1,2) = pr_xy * young_y /(ONE - pr_xy*pr_yx)
	dsdeEl(1,3) = ZERO
	dsdeEl(2,1) = pr_xy * young_y /(ONE - pr_xy*pr_yx)
	dsdeEl(2,2) = young_y/(ONE - pr_xy*pr_yx)
	dsdeEl(2,3) = ZERO
	dsdeEl(3,1) = ZERO
	dsdeEl(3,2) = ZERO
	dsdeEl(3,3) = g_xy

      eps_11_f_t =  sig_11_f_t/dsdeEl(1,1)
      eps_11_f_c =  sig_11_f_c/dsdeEl(1,1)
      eps_22_f_t =  sig_22_f_t/dsdeEl(2,2)
      eps_22_f_c =  sig_22_f_t/dsdeEl(2,2)
      
      

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do

      
c**** Compute current stress 
      call vzero(sigma, 3)
          do i=1,ncomp
             do j=1,ncomp
                sigma(i) = sigma(i) + ((ONE-damage)
     &                     *dsdeEl(i,j)*T_strain(j))
             end do
          end do    
      
 
****** M_tensor *****
      do i=1,ncomp
          do j=1,ncomp
              if (i/=j) then
                  M(i,j)  = ZERO
              endif
          end do
      end do
      
      M(1,1) = 1/(1-d1)
      M(2,2) = 1/(1-d2)
      M(3,3) = 1/sqrt((1-d1)*(1-d2))

      
c****  Compute effective stress      
      call vzero(sigma_eff, 3)
      do i=1,ncomp
          do j=1,ncomp
              sigma_eff(i) = sigma_eff(i) + (M(i,j)*sigma(j))        
          end do
      end do         
      
      
      
      
      
c****Damage initiation criteria
 
      if (sigma_eff(1) >= 0) then                    
          F_f_new = sigma_eff(1)/sig_11_f_t          
      else           
          F_f_new = sigma_eff(1)/sig_11_f_c  
      endif
     
      
      if (sigma_eff(2)>= 0)  then          
          F_m_new = sigma_eff(2)/sig_22_f_t                
      else          
          F_m_new = sigma_eff(2)/sig_22_f_c        
      endif
            

      
c****** To make sure damage initiation criteria is greater than or equal to previous step  ************
      
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
      

          
      
      
      
      if ((F_f <= 1) .AND. (F_m <= 1) ) then
       
          do i=1,ncomp
              sigma(i) = sigma(i)       
          end do 
      
      else
      
        if (sigma_eff(1) >= 0) then
          epsilon_f = eps_11_f_t
          P1 = 25
        else 
          epsilon_f = eps_11_f_c
          P1 = -25
        endif
        
        
        if (sigma_eff(2) >= 0) then
          epsilon_m = eps_22_f_t
          P2 = 25
        else 
          epsilon_m = eps_22_f_c
          P2 = -25
        endif
      

        
	  
	  
        
        if (F_f > 1) then
      
          d1_new =  ONE -  (exp(-P1*(T_strain(1) - epsilon_f)))
          
          if (d1_new >= d1) then
             d1 = d1_new
          else
             d1 = d1
          endif
        else
          d1 = ZERO
        
        endif        

      
        if (F_m > 1) then
      
          d2_new =  ONE -  (exp(-P2*(T_strain(2) - epsilon_m)))
          
          if (d2_new >= d2) then
             d2 = d2_new
          else
             d2 = d2
          endif
        else
          d2 = ZERO       
        endif
      
      

      
       damage_new = 1  - ((1 - d1)*(1 - d2))
      
       if (damage_new >= damage) then
           damage = damage_new
       else
           damage = damage
       endif
 
       
       
       call vzero(sigma, 3)   
       do i=1,ncomp
          do j=1,ncomp
              sigma(i) = sigma(i) + ((ONE-damage)
     &            *dsdeEl(i,j)*T_strain(j))
          end do
       end do       
      
       
c*************** First term of the tangent stiffness (1-d+eta)*C **********          
       do i = 1,mcomp
          do j = 1,mcomp
                  a(i,j) = (ONE - damage)*dsdeEl(i,j) 
          end do
       end do           
                 
c*************** Second term of the tangent stiffness ********
       
        b = (/ P1*exp(-P1*(T_strain(1)-epsilon_f))*(1 - d2), 
     &     P2*exp(-P2*(T_strain(2)-epsilon_m))*(1 - d1), ZERO/)  

      
c****  Double contraction of 4th order elastic stiffness and 2nd order strain *******                    
        do i=1,ncomp
           do j=1,ncomp
               c(i)  = dsdeEl(i,j)*T_strain(j)
           end do
        end do
                       

        
c****  Dyadic product of two 2nd order tensors ******
        do i=1,ncomp
           do j=1,ncomp
                d(i,j)  = c(i)*b(j)
           end do
        end do

        
c**** Tangent stiffness *****
        do i = 1,mcomp
            do j = 1,mcomp
               C_T(i,j) = a(i,j) - d(i,j)
            end do
        end do
        
        
      endif
      
      
      ustatev(1) = damage
      ustatev(2) = d1
      ustatev(3) = d2
      ustatev(4) = F_f
      ustatev(5) = F_m
      
      
      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f <= 1) .AND. (F_m <= 1)) then
          
          do i = 1,mcomp
               do j = 1,mcomp
                   dsdePl(i,j) =  (ONE - damage)*dsdeEl(i,j)
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


      