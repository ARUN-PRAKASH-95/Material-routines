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
	
	
      EXTERNAL         usermatps,usermatps2
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
	
	

	if ((ncomp == 3).AND.(ustatev(2)<=0.9).AND.(ustatev(3)<=0.9)) then
	
		
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
          

	elseif ((ncomp == 3).AND.(ustatev(2)>0.9).AND.(ustatev(3)>0.9)) then
	
		
		call usermatps2 (
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
     &                  id(3,3), C_T_1(mcomp,mcomp),M_inv(mcomp,mcomp),
     &                  C_T_1_a(mcomp), C_T_0(mcomp,mcomp),
     &                  C_T_1_b(mcomp), C_T_2(mcomp,mcomp), 
     &                  C_T_2_a(mcomp),
     &                  C_T_2_b(mcomp),sigma_eff(mcomp) 

                  

      
      integer             i,j,k,l,n,P,ielem,iiter,key
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz, L_c,k1,k2,k3,
     &                    xK, yield,G_c_1,G_c_2,G_c_3,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,d1_new,d2_new
      
      double precision    d3_new,F_f,F_m,F_z,term1,term2,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,D,V,alpha,eta

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

      V    = ustatev(1) 
      d1   = ustatev(2)
      d2   = ustatev(3)
      F_f  = ustatev(4)
      F_m  = ustatev(5)  
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	pr_xy    = prop(3)
	g_xy     = prop(4)
      pr_yx = (young_y * pr_xy) / young_x
      
      
      sig_11_f_t = prop(5)
      sig_11_f_c = prop(6)
      sig_22_f_t = 500e6
      sig_22_f_c = -1500e6
      sig_12_f   = 50e6
      G_c_1      = 12.5e5
      G_c_2      = 1e5
      L_c        = V**(1.0/3.0)
      alpha = ZERO
      
     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (young_x/(ONE - pr_xy*pr_yx))
	dsdeEl(1,2) = (pr_xy * young_y /(ONE - pr_xy*pr_yx))
	dsdeEl(1,3) = ZERO
	dsdeEl(2,1) = (pr_xy * young_y /(ONE - pr_xy*pr_yx))
      dsdeEl(2,2) = (young_y/(ONE - pr_xy*pr_yx))
	dsdeEl(2,3) = ZERO
	dsdeEl(3,1) = ZERO
      dsdeEl(3,2) = ZERO
      dsdeEl(3,3) = g_xy
      tsstif(1) = g_xy
      eta = ZERO
  
      
      eps_11_f_t = sig_11_f_t / (young_x/(ONE - pr_xy*pr_yx))
      eps_11_f_c = sig_11_f_c / (young_x/(ONE - pr_xy*pr_yx))
      eps_22_f_t = sig_22_f_t / (young_y/(ONE - pr_xy*pr_yx))
      eps_22_f_c = sig_22_f_c / (young_y/(ONE - pr_xy*pr_yx))
      eps_12_f   = sig_12_f / g_xy

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do

      
c**** Compute current stress 
      call vzero(sigma_eff, 3)
      do i=1,ncomp
          do j=1,ncomp
              sigma_eff(i) = sigma_eff(i) + (dsdeEl(i,j)*T_strain(j))        
          end do
      end do          
      
     
c****Modified Hashin's Damage initiation criteria
 
      if (sigma_eff(1) >= 0) then                    
          F_f_new = sqrt( (sigma_eff(1)/sig_11_f_t)**2 +       
     & alpha*(sigma_eff(3)/sig_12_f)**2)
      else           
          F_f_new = sqrt( (sigma_eff(1)/sig_11_f_c)**2 +       
     & alpha*(sigma_eff(3)/sig_12_f)**2)
      endif
            
      
      if (sigma_eff(2)>= 0)  then         
          F_m_new = sqrt( (sigma_eff(2)/sig_22_f_t)**2 + 
     & alpha*(sigma_eff(3)/sig_12_f)**2)
      else          
          F_m_new = sqrt( (sigma_eff(2)/sig_22_f_c)**2 + 
     & alpha*(sigma_eff(3)/sig_12_f)**2)               
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
          
      
      call vzero(sigma,3)
      call vzero(C_T_1_a, 3) 
      call vzero(C_T_2_a, 3)
      if ((F_f <= 1) .AND. (F_m <= 1)) then
 
c********    Update stress
         do i=1,ncomp
            sigma(i) = sigma_eff(i)       
         end do  
        
      
      else
          
c*********  Terms in damage evolution equations   ***********
        if (sigma_eff(1) >= 0) then
          k1 =   (-sig_11_f_t*eps_11_f_t*L_c)/G_c_1
        else 
          k1 =   (-sig_11_f_c*eps_11_f_c*L_c)/G_c_1
        endif
        
        
        if (sigma_eff(2)  >= 0) then
          k2 =  (-sig_22_f_t*eps_22_f_t*L_c)/G_c_2
        else 
          k2 =  (-sig_22_f_c*eps_22_f_c*L_c)/G_c_2
        endif
      
        
       

c*********   Damage evolution equations    *********** 	  
        if (F_f > 1) then
      
          d1_new =  ONE  - (exp(k1*(F_f - ONE)))/F_f
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
          
        endif        

      
        if (F_m > 1) then
      
          d2_new =  ONE  - (exp(k2*(F_m - ONE)))/F_m
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
        
        endif
        
     
        d3  = ONE - ((ONE - d1)*(ONE - d2))

c*****      Inverse of the damage effect tensor  ********
        
        do i = 1,ncomp
            do j= 1,ncomp
                M_inv(i,j) = ZERO
            end do
        end do
        
        M_inv(1,1)  = (ONE - d1)
        M_inv(2,2)  = (ONE - d2)
        M_inv(3,3)  = (ONE - d3)
      
c*******   Find nominal stress from effective stress ********
        do i=1,ncomp
           do j=1,ncomp
               sigma(i) = sigma(i) + (M_inv(i,j)*sigma_eff(j))        
           end do
        end do  
        
        do i = 1,ncomp
            do j= 1,ncomp
                C_T_0(i,j) = ZERO
            end do
        end do

c*******   Degraded stiffness  ********
        do i=1,ncomp
           do j=1,ncomp
               do k=1,ncomp
                   C_T_0(i,j) = C_T_0(i,j) + (M_inv(i,k)*dsdeEl(k,j))
               end do
           end do
        end do

        
        
        dsdeEl_d(1,1) = (young_x/(ONE - d3*pr_xy*pr_yx))*(ONE+eta- d1)
        dsdeEl_d(1,2) =(pr_yx*young_x/(ONE-d3*pr_xy*pr_yx))*(ONE+eta-d1)
     &                   *(ONE+eta-d2)
        dsdeEl_d(1,3) = ZERO
        dsdeEl_d(2,1) =(pr_xy*young_y/(ONE-d3*pr_xy*pr_yx))*(ONE+eta-d2)
     &                   *(ONE+eta-d1)
        dsdeEl_d(2,2) = (young_y/(ONE - d3*pr_xy*pr_yx))*(ONE+eta- d2)
        dsdeEl_d(2,3) = ZERO
        dsdeEl_d(3,1) = ZERO
        dsdeEl_d(3,2) = ZERO
        dsdeEl_d(3,3) = g_xy*(ONE+eta-d3)
        
        
        if ((d1 == ZERO)) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_1(i,j) = 0        
             end do
          end do  
          
        else
            
          C_T_1_a = (/-ONE*sigma_eff(1),ZERO,(ONE - d2)*sigma_eff(3)/)

          
c**********   Derivative of d1 with respect to strain (d_d1/d_epsilon)   ***********
 
c**       For Tension          
          if (sigma_eff(1)> 0)  then
              
            term1 = ((ONE-k1*F_f)*exp(k1*(F_f - ONE)))/(F_f**2)
          
            C_T_1_b  = (/(term1*sigma_eff(1)*dsdeEl(1,1))
     &   /(F_f*sig_11_f_t**2), (term1*sigma_eff(1)*dsdeEl(1,2))
     &   /(F_f*sig_11_f_t**2), (alpha*sigma_eff(3)*term1*dsdeEl(3,3))
     &    /(F_f*sig_12_f**2)/)        
c**       For Compression           
          else
              
            term1 = ((ONE-k1*F_f)*exp(k1*(F_f - ONE)))/(F_f**2)
          
            C_T_1_b  = (/(term1*sigma_eff(1)*dsdeEl(1,1))
     &   /(F_f*sig_11_f_c**2), (term1*sigma_eff(1)*dsdeEl(1,2))
     &   /(F_f*sig_11_f_c**2), (alpha*sigma_eff(3)*term1*dsdeEl(3,3))
     &    /(F_f*sig_12_f**2)/)     
            
          endif
      
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_1(i,j) = C_T_1_a(i)*C_T_1_b(j)
              end do
          end do
       
        endif
            
        

        if ((d2 == ZERO)) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_2(i,j) = 0        
             end do
          end do 
      
        else 
            
          C_T_2_a = (/ZERO,-ONE*sigma_eff(2),(ONE - d1)*sigma_eff(3)/)

          
c******     Derivative of d2 with respect to strain (d_d2/d_epsilon)
      
c****      For Tension           
          if (sigma_eff(2)>= 0)  then
          
            term2 = ((ONE-k2*F_m)*exp(k2*(F_m - ONE)))/(F_m**2)
          
            C_T_2_b  = (/(term2*sigma_eff(2)*dsdeEl(2,1))
     &   /(F_m*sig_22_f_t**2), (term2*sigma_eff(2)*dsdeEl(2,2))
     &   /(F_m*sig_22_f_t**2), (alpha*sigma_eff(3)*term2*dsdeEl(3,3))
     &    /(F_m*sig_12_f**2)/) 
c**       For Compression                       
          else                                         
              
            term2 = ((ONE-k2*F_m)*exp(k2*(F_m - ONE)))/(F_m**2)
          
            C_T_2_b  = (/(term2*sigma_eff(2)*dsdeEl(2,1))
     &   /(F_m*sig_22_f_c**2), (term2*sigma_eff(2)*dsdeEl(2,2))
     &   /(F_m*sig_22_f_c**2), (alpha*sigma_eff(3)*term2*dsdeEl(3,3))
     &    /(F_m*sig_12_f**2)/) 
          endif

              
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_2(i,j) = C_T_2_a(i)*C_T_2_b(j)
              end do
          end do              
              
        endif      
      

        do i = 1,ncomp
            do j= 1,ncomp
                C_T(i,j) = ZERO
            end do
        end do
        
        do i=1,ncomp
            do j=1,ncomp
                  C_T(i,j)   = dsdeEl_d(i,j)  +  C_T_1(i,j) + C_T_2(i,j) 
            end do
        end do
           
       

        open(10,file='data.txt')
        write(10,*)  ldstep
        write(10,*)  isubst
        write(10,*)  '        '
        write(10,*)  'C_matrix'
        do i=1,ncomp
             write(10,*) (dsdeEl(i,j), j = 1,ncomp)
        enddo
        
        write(10,*)  '        '
        write(10,*) d1
        write(10,*) d2
        write(10,*) d3
        
        write(10,*)  '        '
        write(10,*)  'Stress'
        do i=1,ncomp
             write(10,*) sigma(i)
        enddo  
        
       
        write(10,*)  '        '
        write(10,*)  'C_T_0 matrix'
        do i=1,ncomp
             write(10,*) (C_T_0(i,j), j = 1,ncomp)
        enddo 
        
        write(10,*)  '        '
        write(10,*)  'C_T_2_a'
        do i=1,ncomp
             write(10,*) C_T_2_a(i)
        enddo
        
        write(10,*)  '        '
        write(10,*)  'C_T_2_b'
        do i=1,ncomp
             write(10,*) C_T_2_b(i)
        enddo
        
        write(10,*)  '        '
        write(10,*)  'C_T_2 matrix'
        do i=1,ncomp
             write(10,*) (C_T_2(i,j), j = 1,ncomp)
        enddo 
        
        write(10,*)  '        '
        write(10,*)  'C_T matrix'       
        do i=1,ncomp
             write(10,*) (C_T(i,j), j = 1,ncomp)
        enddo 
                  
      endif
      

      ustatev(2) = d1
      ustatev(3) = d2
      ustatev(4) = F_f
      ustatev(5) = F_m


      
 
      call vzero(stress,3)      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f <= 1) .AND. (F_m <= 1)) then
          
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


**************    Subroutine for 3d elements (Linear elastic isotropic model)   ********************
      subroutine usermatps2(
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
     &                  id(3,3), C_T_1(mcomp,mcomp),M_inv(mcomp,mcomp),
     &                  C_T_1_a(mcomp), C_T_0(mcomp,mcomp),
     &                  C_T_1_b(mcomp), C_T_2(mcomp,mcomp), 
     &                  C_T_2_a(mcomp),
     &                  C_T_2_b(mcomp),sigma_eff(mcomp) 

                  

      
      integer             i,j,k,l,n,P,ielem,iiter,key
      double precision    young_x, young_y, young_z,  
     &                    pr_xy, pr_yz, pr_xz, L_c,k1,k2,k3,
     &                    xK, yield,G_c_1,G_c_2,G_c_3,
     &                    eps_11_f_t,eps_11_f_c,eps_22_f_t,
     &                    eps_22_f_c,eps_33_f_t,eps_33_f_c,
     &                    eps_12_f,eps_13_f,eps_23_f,sig_11_f_t,
     &                    sig_22_f_t,sig_33_f_t,d1,d2,d3,d1_new,d2_new
      
      double precision    d3_new,F_f,F_m,F_z,term1,term2,g_xy, g_yz, 
     &                    g_xz,delta,pr_yx, pr_zy, pr_zx,xy_yx, yz_zy,
     &                    zx_xz,xyz, E_xyz, sig_11_f_c, sig_22_f_c,
     &                    sig_33_f_c,sig_12_f,sig_23_f,sig_13_f,F_f_new,
     &                    F_m_new,F_z_new,D,V,alpha,eta

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

      V    = ustatev(1) 
      d1   = ustatev(2)
      d2   = ustatev(3)
      F_f  = ustatev(4)
      F_m  = ustatev(5)  
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	pr_xy    = prop(3)
	g_xy     = prop(4)
      pr_yx = (young_y * pr_xy) / young_x
      
      
      sig_11_f_t = prop(5)
      sig_11_f_c = prop(6)
      sig_22_f_t = 500e6
      sig_22_f_c = -1500e6
      sig_12_f   = 50e6
      G_c_1      = 12.5e5
      G_c_2      = 1e5
      L_c        = V**(1.0/3.0)
      alpha = ZERO
      
     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (young_x/(ONE - pr_xy*pr_yx))
	dsdeEl(1,2) = (pr_xy * young_y /(ONE - pr_xy*pr_yx))
	dsdeEl(1,3) = ZERO
	dsdeEl(2,1) = (pr_xy * young_y /(ONE - pr_xy*pr_yx))
      dsdeEl(2,2) = (young_y/(ONE - pr_xy*pr_yx))
	dsdeEl(2,3) = ZERO
	dsdeEl(3,1) = ZERO
      dsdeEl(3,2) = ZERO
      dsdeEl(3,3) = g_xy
      tsstif(1) = g_xy
      eta = ZERO
  
      
      eps_11_f_t = sig_11_f_t / (young_x/(ONE - pr_xy*pr_yx))
      eps_11_f_c = sig_11_f_c / (young_x/(ONE - pr_xy*pr_yx))
      eps_22_f_t = sig_22_f_t / (young_y/(ONE - pr_xy*pr_yx))
      eps_22_f_c = sig_22_f_c / (young_y/(ONE - pr_xy*pr_yx))
      eps_12_f   = sig_12_f / g_xy

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) 
      end do
     
c**** Compute current stress 
      call vzero(sigma_eff, 3)
      do i=1,ncomp
          do j=1,ncomp
              sigma_eff(i) = sigma_eff(i) + (dsdeEl(i,j)*T_strain(j))        
          end do
      end do         
      
      
      call vzero(sigma,3)
      call vzero(C_T_1_a, 3) 
      call vzero(C_T_2_a, 3)
      if ((F_f <= 1) .AND. (F_m <= 1)) then
 
c********    Update stress
         do i=1,ncomp
            stress(i) = sigma_eff(i)       
         end do  
        
      
      else
          
c*********  Terms in damage evolution equations   ***********
        if (sigma_eff(1) >= 0) then
          k1 =   (-sig_11_f_t*eps_11_f_t*L_c)/G_c_1
        else 
          k1 =   (-sig_11_f_c*eps_11_f_c*L_c)/G_c_1
        endif
        
        
        if (sigma_eff(2)  >= 0) then
          k2 =  (-sig_22_f_t*eps_22_f_t*L_c)/G_c_2
        else 
          k2 =  (-sig_22_f_c*eps_22_f_c*L_c)/G_c_2
        endif
      
        
       

c*********   Damage evolution equations    *********** 	  
        if (F_f > 1) then
      
          d1_new =  d1
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
          
        endif        

      
        if (F_m > 1) then
      
          d2_new =  d2
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
        
        endif


        
        d3  = ONE - ((ONE - d1)*(ONE - d2))
        
c*****      Inverse of the damage effect tensor  ********
        
        do i = 1,ncomp
            do j= 1,ncomp
                M_inv(i,j) = ZERO
            end do
        end do
        
        M_inv(1,1)  = (ONE - d1)
        M_inv(2,2)  = (ONE - d2)
        M_inv(3,3)  = (ONE - d3)
      
c*******   Find nominal stress from effective stress ********
        do i=1,ncomp
           do j=1,ncomp
               sigma(i) = sigma(i) + (M_inv(i,j)*sigma_eff(j))        
           end do
        end do 
        
        
      
        dsdeEl_d(1,1) = (young_x/(ONE - d3*pr_xy*pr_yx))*(ONE+eta- d1)
        dsdeEl_d(1,2) =(pr_yx*young_x/(ONE-d3*pr_xy*pr_yx))*(ONE+eta-d1)
     &                   *(ONE+eta-d2)
        dsdeEl_d(1,3) = ZERO
        dsdeEl_d(2,1) =(pr_xy*young_y/(ONE-d3*pr_xy*pr_yx))*(ONE+eta-d2)
     &                   *(ONE-d1)
        dsdeEl_d(2,2) = (young_y/(ONE - d3*pr_xy*pr_yx))*(ONE+eta- d2)
        dsdeEl_d(2,3) = ZERO
        dsdeEl_d(3,1) = ZERO
        dsdeEl_d(3,2) = ZERO
        dsdeEl_d(3,3) = g_xy*(ONE+eta-d3)
        
        
        if ((d1 == ZERO)) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_1(i,j) = 0        
             end do
          end do  
          
        else
            
          C_T_1_a = (/-ONE*sigma_eff(1),ZERO,(ONE - d2)*sigma_eff(3)/)

          
c**********   Derivative of d1 with respect to strain (d_d1/d_epsilon)   ***********
 
c**       For Tension          
          if (sigma_eff(1)> 0)  then
              
            term1 = ((ONE-k1*F_f)*exp(k1*(F_f - ONE)))/(F_f**2)
          
            C_T_1_b  = (/(term1*sigma_eff(1)*dsdeEl(1,1))
     &   /(F_f*sig_11_f_t**2), (term1*sigma_eff(1)*dsdeEl(1,2))
     &   /(F_f*sig_11_f_t**2), (alpha*sigma_eff(3)*term1*dsdeEl(3,3))
     &    /(F_f*sig_12_f**2)/)        
c**       For Compression           
          else
              
            term1 = ((ONE-k1*F_f)*exp(k1*(F_f - ONE)))/(F_f**2)
          
            C_T_1_b  = (/(term1*sigma_eff(1)*dsdeEl(1,1))
     &   /(F_f*sig_11_f_c**2), (term1*sigma_eff(1)*dsdeEl(1,2))
     &   /(F_f*sig_11_f_c**2), (alpha*sigma_eff(3)*term1*dsdeEl(3,3))
     &    /(F_f*sig_12_f**2)/)     
            
          endif
      
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_1(i,j) = C_T_1_a(i)*C_T_1_b(j)
              end do
          end do
       
        endif
            
        

        if ((d2 == ZERO)) then
          
          do i=1,mcomp
             do j=1,mcomp
                 C_T_2(i,j) = 0        
             end do
          end do 
      
        else 
            
          C_T_2_a = (/ZERO,-ONE*sigma_eff(2),(ONE - d1)*sigma_eff(3)/)

          
c******     Derivative of d2 with respect to strain (d_d2/d_epsilon)
      
c****      For Tension           
          if (sigma_eff(2)> 0)  then
          
            term2 = ((ONE-k2*F_m)*exp(k2*(F_m - ONE)))/(F_m**2)
          
            C_T_2_b  = (/(term2*sigma_eff(2)*dsdeEl(2,1))
     &   /(F_m*sig_22_f_t**2), (term2*sigma_eff(2)*dsdeEl(2,2))
     &   /(F_m*sig_22_f_t**2), (alpha*sigma_eff(3)*term2*dsdeEl(3,3))
     &    /(F_m*sig_12_f**2)/) 
c**       For Compression                       
          else                                         
              
            term2 = ((ONE-k2*F_m)*exp(k2*(F_m - ONE)))/(F_m**2)
          
            C_T_2_b  = (/(term2*sigma_eff(2)*dsdeEl(2,1))
     &   /(F_m*sig_22_f_c**2), (term2*sigma_eff(2)*dsdeEl(2,2))
     &   /(F_m*sig_22_f_c**2), (alpha*sigma_eff(3)*term2*dsdeEl(3,3))
     &    /(F_m*sig_12_f**2)/) 
          endif

              
          do i = 1,mcomp
              do j = 1,mcomp
                  C_T_2(i,j) = C_T_2_a(i)*C_T_2_b(j)
              end do
          end do              
              
        endif      
      

        do i = 1,ncomp
            do j= 1,ncomp
                C_T(i,j) = ZERO
            end do
        end do
        
        do i=1,ncomp
            do j=1,ncomp
                  C_T(i,j)   = dsdeEl_d(i,j)  +  C_T_1(i,j) + C_T_2(i,j) 
            end do
        end do
           
       

        open(10,file='data2.txt')
        
        write(10,*)  ldstep
        write(10,*)  isubst
        write(10,*)  '        '
        write(10,*)  'C_matrix'
        do i=1,ncomp
             write(10,*) (dsdeEl(i,j), j = 1,ncomp)
        enddo
        
        write(10,*)  '        '
        write(10,*) d1
        write(10,*) d2
        write(10,*) d3
        
        write(10,*)  '        '
        write(10,*)  'Stress'
        do i=1,ncomp
             write(10,*) sigma(i)
        enddo  
        
       
        write(10,*)  '        '
        write(10,*)  'C_T_0 matrix'
        do i=1,ncomp
             write(10,*) (C_T_0(i,j), j = 1,ncomp)
        enddo 

        write(10,*)  '        '
        write(10,*)  'C_T_2_a'
        do i=1,ncomp
             write(10,*) C_T_2_a(i)
        enddo
        
        write(10,*)  '        '
        write(10,*)  'C_T_2_b'
        do i=1,ncomp
             write(10,*) C_T_2_b(i)
        enddo
        
        write(10,*)  '        '
        write(10,*)  'C_T_2 matrix'
        do i=1,ncomp
             write(10,*) (C_T_2(i,j), j = 1,ncomp)
        enddo 
                 
      endif
      

      ustatev(2) = d1
      ustatev(3) = d2
      ustatev(4) = F_f
      ustatev(5) = F_m


      call vzero(stress,3)
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      
      write(10,*)  '        '
      write(10,*)  'Stress'       
      do i=1,ncomp
          write(10,*) (Stress(i))
      enddo       

      
      if ((F_f <= 1) .AND. (F_m <= 1)) then
          
          do i=1,ncomp
              do j=1,ncomp
                  dsdePl(i,j) = C_T(i,j) 
              end do
          end do 
      
      else
          
          do i=1,ncomp
              do j=1,ncomp
                  dsdePl(i,j) = C_T(i,j)  
              end do
          end do 
          
      endif  
      
      write(10,*)  '        '
      write(10,*)  'C_T matrix'       
      do i=1,ncomp
          write(10,*) (dsdePl(i,j), j = 1,ncomp)
      enddo       
      
      return
      end



              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      