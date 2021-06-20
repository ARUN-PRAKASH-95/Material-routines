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
     &           isubst,keycut, nDirect,nShear,ncomp,nStatev,nProp,ttype


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
	    ttype = 1
		
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
     &                   var6, var7, var8, ttype)


      end if
      return
	end
	
	
**************    Subroutine for 3d elements (Linear elastic isotropic model)   ********************
      recursive subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8, ttype)
	
	
#include "impcom.inc"
#include "ansysdef.inc"
c
      integer :: matId, elemId, kDomIntPt, kLayer, kSectPt, ldstep,
     &           isubst,keycut, nDirect,nShear,ncomp,nStatev,nProp,ttype


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
   
      double precision  dsdeEl(ncomp,ncomp), sigma(ncomp),
     &                  C_T(ncomp,ncomp), C_T_num(ncomp,ncomp),
     &                  dsdeEl_d(ncomp,ncomp), T_strain(ncomp),
     &                  id(3,3), C_T_1(ncomp,ncomp),M_inv(ncomp,ncomp),
     &                  C_T_1_a(ncomp),sigma_eff(ncomp),
     &                  C_T_1_b(ncomp), C_T_2(ncomp,ncomp), 
     &                  C_T_2_a(ncomp),C_T_2_b(ncomp),C_T_0(ncomp,ncomp) 
                  

      
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
     &                    F_m_new,F_z_new,D,V,eta,hper,dmax,multiplier

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

 
      V    = ustatev(1) 
      d1   = ustatev(2)
      d2   = ustatev(3)
      d3   = ONE - ((ONE - d1)*(ONE - d2))
      F_f  = ustatev(4)
      F_m  = ustatev(5) 
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
	young_x  = prop(1)
	young_y  = prop(2)
	pr_xy    = prop(3)
	g_xy     = prop(4)
      pr_yx = (young_y * pr_xy) / young_x
      
      
      sig_11_f_t = 2500e6
      sig_11_f_c = -2000e6
      sig_22_f_t = 500e6
      sig_22_f_c = -1500e6
      sig_12_f   = 50e6
      G_c_1      = 12.5e7
      G_c_2      = 1e7
      G_c_3      = 1e7
      L_c        = V**(1.0/3.0)
	eta        =  0
      dmax = 1
      multiplier = 1e-3
     
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
        if (F_f > ONE) then
      
          d1_new = ONE  - (exp(k1*(F_f - ONE)))/F_f 
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
          
        endif        

      
        if (F_m > ONE) then
      
          d2_new = ZERO
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
        
        endif
      
        d3  = ONE - ((ONE -d1)*(ONE - d2))
          
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
        
c*******   Degraded stiffness  ********
        do i=1,ncomp
           do j=1,ncomp
               do k=1,ncomp
                   C_T_0(i,j) = C_T_0(i,j) + (M_inv(i,k)*dsdeEl(k,j))
               end do
           end do
        end do
        
        dsdeEl_d(1,1) = (young_x/(ONE - pr_xy*pr_yx))*(ONE + eta - d1)
        dsdeEl_d(1,2) = (pr_xy *young_y/(ONE-pr_xy*pr_yx))*(ONE+eta-d1)
        dsdeEl_d(1,3) = ZERO
        dsdeEl_d(2,1) = (pr_xy *young_y/(ONE-pr_xy*pr_yx))*(ONE+ eta-d2)
        dsdeEl_d(2,2) = (young_y/(ONE - pr_xy*pr_yx))*(ONE+ eta - d2)
        dsdeEl_d(2,3) = ZERO
        dsdeEl_d(3,1) = ZERO
        dsdeEl_d(3,2) = ZERO
        dsdeEl_d(3,3) = g_xy*(ONE+eta-d3)

       
        if (d1 == ZERO) then
          
          do i=1,ncomp
             do j=1,ncomp
                 C_T_1(i,j) = 0   
             end do
          end do  
          
        else
            
            
            
          C_T_1_a = (/-ONE*sigma_eff(1),ZERO,(ONE - d2)*sigma_eff(3)/)
            
c**********   Derivative of d1 with respect to strain (d_d1/d_epsilon)   ***********
 
c**       For Tension          
          if (sigma_eff(1)> 0)  then
          
            C_T_1_b  = (/((ONE - k1*F_f)/(F_f**2 * sig_11_f_t*(ONE-d1)))
     &     *exp(k1*(F_f - ONE))*dsdeEl(1,1),  ((ONE - k1*F_f)/(F_f**2 
     &     * sig_11_f_t*(ONE-d1)))*exp(k1*(F_f -ONE))*dsdeEl(1,2),ZERO/)
            
c**       For Compression           
          else
              
            C_T_1_b  = (/((ONE - k1*F_f)/(F_f**2 * sig_11_f_c*(ONE-d1)))
     &     *exp(k1*(F_f - ONE))*dsdeEl(1,1),  ((ONE - k1*F_f)/(F_f**2 
     &     * sig_11_f_c*(ONE-d1)))*exp(k1*(F_f -ONE))*dsdeEl(1,2),ZERO/)
            
          endif
      
          do i = 1,ncomp
              do j = 1,ncomp
                  C_T_1(i,j) = C_T_1_a(i)*C_T_1_b(j)
              end do
          end do
       
        endif
            
        

        if (d2 == ZERO) then
          
          do i=1,ncomp
             do j=1,ncomp
                 C_T_2(i,j) = ZERO     
             end do
          end do 
      
        else 
            
          C_T_2_a = (/ZERO,-ONE*sigma_eff(2),(ONE - d1)*sigma_eff(3)/)

          
c******     Derivative of d2 with respect to strain (d_d2/d_epsilon)
      
c****      For Tension           
          if (sigma_eff(2)> 0)  then
          
            C_T_2_b  = (/((ONE - k2*F_m)/(F_m**2 * sig_22_f_t*(ONE-d2)))
     &     *exp(k2*(F_m - ONE))*dsdeEl(2,1),  ((ONE - k2*F_m)/(F_m**2 
     &     *sig_22_f_t*(ONE-d2)))*exp(k2*(F_m - ONE))*dsdeEl(2,2),ZERO/)
     
c**       For Compression                       
          else                                         
              
            C_T_2_b  = (/((ONE - k2*F_m)/(F_m**2 * sig_22_f_c*(ONE-d2)))
     &     *exp(k2*(F_m - ONE))*dsdeEl(2,1),  ((ONE - k2*F_m)/(F_m**2 
     &     *sig_22_f_c*(ONE-d2)))*exp(k2*(F_m - ONE))*dsdeEl(2,2),ZERO/) 
              
          endif

              
          do i = 1,ncomp
              do j = 1,ncomp
                  C_T_2(i,j) = C_T_2_a(i)*C_T_2_b(j)
              end do
          end do              
              
        endif   
        
************   Algorithmic tangent   ************
         do i=1,ncomp
            do j=1,ncomp
               C_T(i,j) = dsdeEl_d(i,j) + C_T_1(i,j) + C_T_2(i,j)
            end do
         end do
            

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
      

      if (ttype == 0) then
          
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
          
      elseif (ttype == 1) then 
          
          
          if ((F_f <= 1) .AND. (F_m <= 1)) then
          
              do i=1,ncomp
                  do j=1,ncomp
                      dsdePl(i,j) = dsdeEl(i,j)        
                  end do
              end do 
      
          else
          
              do i = 1,ncomp
                  hper = 1e-12
                  dStrain(i) = dStrain(i) + hper
                  ttype = 0
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
     &                   var6, var7, var8, ttype)
      
      
              
                  do j = 1,ncomp
                  
                      C_T_num(j,i) = (stress(j) - sigma(j))/hper
                  end do
              end do
              
              if (d1<=dmax) then
                  do i=1,ncomp
                      do j=1,ncomp
                          dsdePl(i,j) = C_T_num(i,j)       
                      end do
                  end do 
              else
                  do i=1,ncomp
                      do j=1,ncomp
                          dsdePl(i,j) = multiplier*C_T_num(i,j)       
                      end do
                  end do               
              endif
              
              
          endif
          
      endif
                  
                  
      
      return
      end



              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      