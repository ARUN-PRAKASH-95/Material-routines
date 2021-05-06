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
     &                  a(mcomp,mcomp), b(mcomp), d(mcomp,mcomp),
     &                  c(mcomp)
                  

      
      integer             i,j,k,l,n,P1,P2,P3
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz,epsilon_f, 
     &                    epsilon_m,epsilon_z, damage,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,damage_new
      
      double precision    d3_new,F_f,F_m,F_z,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,V,key,ans,V2,value

      
      
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
      sig_22_f_c = -150e6
      sig_33_f_t = 50e6
      sig_33_f_c = -150e6
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

      
      damage = ustatev(1)


     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (ONE -yz_zy) / (young_y*young_z*delta)
	dsdeEl(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)
	dsdeEl(1,3) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta)
	dsdeEl(1,4) = 0
	dsdeEl(1,5) = 0
	dsdeEl(1,6) = 0
	dsdeEl(2,1) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)
	dsdeEl(2,2) = (ONE -zx_xz) / (young_x*young_z*delta)
	dsdeEl(2,3) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta)
	dsdeEl(2,4) = 0
	dsdeEl(2,5) = 0
	dsdeEl(2,6) = 0
	dsdeEl(3,1) = (pr_zx + pr_yx*pr_zy) / (young_y*young_z*delta)
	dsdeEl(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta)
	dsdeEl(3,3) = (ONE -xy_yx) / (young_x*young_y*delta)
	dsdeEl(3,4) = 0
	dsdeEl(3,5) = 0
	dsdeEl(3,6) = 0
	dsdeEl(4,1) = 0
	dsdeEl(4,2) = 0 
	dsdeEl(4,3) = 0
	dsdeEl(4,4) = g_xy
	dsdeEl(4,5) = 0
	dsdeEl(4,6) = 0
	dsdeEl(5,1) = 0
	dsdeEl(5,2) = 0 
	dsdeEl(5,3) = 0
	dsdeEl(5,4) = 0
	dsdeEl(5,5) = g_yz
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = g_xz

  
      
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
          
          
          F_f = T_strain(1)/eps_11_f_t
          
      else 
          
          F_f = T_strain(1)/eps_11_f_c
  
      endif
      
      
      
      if (T_strain(2)>= 0)  then
          
          F_m = T_strain(2)/eps_22_f_t 
          
      
      else
          
          F_m = T_strain(2)/eps_22_f_c  
          
      
      endif
      
      
      
      if (T_strain(3) >= 0 )then
      
          F_z = T_strain(3)/eps_33_f_t
          
      else 
      
          F_z = T_strain(3)/eps_33_f_c
          
          
      endif
      
      
      
      call vzero(sigma, 6)
      if ((F_f < 1) .AND. (F_m < 1) .AND. (F_z < 1)) then
       
          do i=1,ncomp
             do j=1,ncomp
                sigma(i) = sigma(i) + ((ONE-damage)
     &                     *dsdeEl(i,j)*T_strain(j))
             end do
          end do
      
      else
      
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
       
        
	  
	  
        
        if (F_f >= 1) then
      
          d1 =  ONE -  (exp(-P1*(T_strain(1) - epsilon_f)))
        else
          d1 = ZERO
        
        endif        

      
        if (F_m >= 1) then
      
          d2 =  ONE -  (exp(-P2*(T_strain(2) - epsilon_m)))
        else
          d2 = ZERO       
        endif
      
      
        if (F_z >= 1) then
      
          d3 =  ONE -  (exp(-P3*(T_strain(3) - epsilon_z)))
        else
          d3 = ZERO 
        endif 
      
       damage_new = 1  - ((1 - d1)*(1 - d2)*(1 - d3))
      
       if (damage_new >= damage) then
           damage = damage_new
       else
           damage = damage
       endif
           
         
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
       
        b = (/ P1*exp(-P1*(T_strain(1)-epsilon_f))*(1 - d2)*(1 - d3), 
     &     P2*exp(-P2*(T_strain(2)-epsilon_m))*(1 - d1)*(1 - d3), 
     &     P3*exp(-P3*(T_strain(3)-epsilon_z))*(1 - d2)*(1 - d1), 
     &     ZERO, ZERO, ZERO  /)  

      
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
   
      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f < 1) .AND. (F_m < 1) .AND. (F_z < 1)) then
          
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


      