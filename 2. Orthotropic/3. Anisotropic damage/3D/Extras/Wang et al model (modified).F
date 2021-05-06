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
     &                  C_T_3_b(mcomp),arr(mcomp),value(11)
                  

      
      integer             i,j,k,l,n,P,ielem,iiter
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
     &                    F_m_new,F_z_new,V,eta

      
      
      double precision    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	

      
      integer :: iott,wrinqr
      external wrinqr
c

	
      id = transpose(reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), 
     & shape(id)))


 
      
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
      G_c_1      = 12.5e3
      G_c_2      = 1e3
      G_c_3      = 1e3

      
      L_c        = 1**(1/3)
      
   
      
	pr_yx = (young_y * pr_xy) / young_x
	pr_zy = (young_z * pr_yz) / young_y
	pr_zx = (young_z * pr_xz) / young_x
	
	
	xy_yx = pr_xy*pr_yx
	yz_zy = pr_yz*pr_zy 
	zx_xz = pr_zx*pr_xz
	xyz   = TWO*pr_xy*pr_yz*pr_zx
	E_xyz = young_x*young_y*young_z
	
	delta = (ONE - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz

      
      d1   = ustatev(1)
      d2   = ustatev(2)
      d3   = ustatev(3)
      F_f  = ustatev(4)
      F_m  = ustatev(5)
      F_z  = ustatev(6)
	


     
********* ELASTIC STIFFNESS MATRIX *********
	
	dsdeEl(1,1) = (ONE -yz_zy) / (young_y*young_z*delta)*(ONE - d1)**2
	dsdeEl(1,2) = (pr_yx + pr_zx*pr_yz) / (young_y*young_z*delta)
     &              *(ONE - d1)*(ONE - d2)
	dsdeEl(1,3) = (pr_zx + pr_yx*pr_zy) /(young_y*young_z*delta)*(ONE - d1)
     &              *(ONE - d3)
      dsdeEl(1,4) = 0
	dsdeEl(1,5) = 0
	dsdeEl(1,6) = 0
	dsdeEl(2,1) = (pr_yx + pr_zx*pr_yz) /(young_y*young_z*delta)*(ONE - d1)
     &              *(ONE - d2)
      dsdeEl(2,2) = (ONE -zx_xz) / (young_x*young_z*delta)*(ONE - d2)**2
	dsdeEl(2,3) = (pr_zy + pr_zx*pr_xy) /(young_x*young_z*delta)*(ONE - d3)
     &             *(ONE - d2)
      dsdeEl(2,4) = 0
	dsdeEl(2,5) = 0
	dsdeEl(2,6) = 0
	dsdeEl(3,1) = (pr_zx + pr_yx*pr_zy) /(young_y*young_z*delta)*(ONE - d3)
     &              *(ONE - d1)
      dsdeEl(3,2) = (pr_zy + pr_zx*pr_xy) / (young_x*young_z*delta)*
     &              (ONE - d3)*(ONE - d2)
      dsdeEl(3,3) = (ONE -xy_yx) / (young_x*young_y*delta)*(ONE - d3)**2
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
	dsdeEl(5,5) = g_yz*(ONE - d1)*(ONE - d3)
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = g_xz*(ONE - d3)*(ONE - d2)

  
      
      eps_11_f_t = sig_11_f_t / ((ONE -yz_zy) / (young_y*young_z*delta))
      eps_11_f_c = sig_11_f_c / ((ONE -yz_zy) / (young_y*young_z*delta))
      eps_22_f_t = sig_22_f_t / ((ONE -zx_xz) / (young_x*young_z*delta))
      eps_22_f_c = sig_22_f_c / ((ONE -zx_xz) / (young_x*young_z*delta))
      eps_33_f_t = sig_33_f_t / ((ONE -xy_yx) / (young_x*young_y*delta))
      eps_33_f_c = sig_33_f_c / ((ONE -xy_yx) / (young_x*young_y*delta))
      eps_12_f   = sig_12_f / g_xy
      eps_13_f   = sig_13_f / g_yz
      eps_23_f   = sig_23_f / g_xz

c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do

 
c****Damage initiation criteria
 
      if (T_strain(1) >= 0) then
          
          
          F_f_new = SQRT((T_strain(1)/eps_11_f_t)**2 + (T_strain(4) 
     &       /eps_12_f)**2   + (T_strain(5)/eps_13_f)**2)
          
      else 
          
          F_f_new = SQRT((T_strain(1)/eps_11_f_c)**2)
  
      endif
      
      
      
      if (T_strain(2) + T_strain(3)>= 0)  then
          
          F_m_new =   ((T_strain(2) + T_strain(3))**2/(eps_22_f_t
     &    *eps_33_f_t)) -  (T_strain(2)*T_strain(3)/eps_23_f**2) +  
     &    (T_strain(4)/eps_12_f)**2 + (T_strain(5)/eps_13_f)**2 +
     &     (T_strain(6)/eps_23_f)**2  
          
      
      else
          
          F_m_new =   ((T_strain(2) + T_strain(3))**2/(eps_22_f_c
     &    *eps_33_f_c)) + ((T_strain(2)+T_strain(3)/eps_22_f_c)  
     &    *((eps_22_f_c/2*eps_12_f)  - ONE)) -  
     &    (T_strain(2)*T_strain(3)/eps_23_f**2) +  
     &    (T_strain(4)/eps_12_f)**2 + (T_strain(5)/eps_13_f)**2 +
     &     (T_strain(6)/eps_23_f)**2 
          
      
      endif
      
      
      
      if (T_strain(3) >= 0 )then
      
          F_z_new = SQRT((T_strain(3)/eps_33_f_t)**2 + (T_strain(5) 
     &   /eps_13_f)**2  + (T_strain(6)/eps_23_f)**2) 
          
      else 
      
          F_z_new = SQRT((T_strain(3)/eps_33_f_c)**2 + (T_strain(5) 
     &   /eps_13_f)**2  + (T_strain(6)/eps_23_f)**2)
          
          
      endif
      
      
      
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
      if ((F_f**2 < 1) .AND. (F_m < 1) .AND. (F_z**2 < 1)) then
       
          do i=1,ncomp
             do j=1,ncomp
                 sigma(i) = sigma(i) + (dsdeEl(i,j)*T_strain(j))        
             end do
          end do  
       
      
      else
      
        if (T_strain(1) > 0) then
          k1 =  -0.01
        else 
          k1 =  -1
        endif
        
        
        if (T_strain(2) > 0) then
          k2 =  (-sig_22_f_t*eps_22_f_t*L_c)/G_c_2
        else 
          k2 =  (-sig_22_f_c*eps_22_f_c*L_c)/G_c_2
        endif
      
        
        if (T_strain(3) > 0) then
          k3 =  (-sig_33_f_t*eps_33_f_t*L_c)/G_c_3
        else 
          k3 =  (-sig_33_f_c*eps_33_f_c*L_c)/G_c_3
        endif
       
        
        if (F_f**2 >= 1) then
      
          d1_new =  ONE  - (exp(k1*(F_f - ONE)))/F_f
      
          if (d1_new > d1) then
              d1 = d1_new
          else
              d1 = d1
          endif
      
   
        endif        

      
        if (F_m >= 1) then
      
          d2_new =  ONE   - (exp(k2*(sqrt(F_m) - ONE)))/sqrt(F_m)
      
          if (d2_new > d2) then
              d2 = d2_new
          else
              d2 = d2
          endif
          
  
        endif
      
      
        if (F_z**2 >= 1) then
      
          d3_new =  ONE   - (exp(k3*(F_z - ONE)))/F_z   
      
          if (d3_new > d3) then
              d3 = d3_new
          else
              d3 = d3
          endif

        endif 
      
      
        dsdeEl_d(1,1)=(ONE -yz_zy)/(young_y*young_z*delta)*(ONE - d1)**2
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
        dsdeEl_d(3,3) = (ONE -xy_yx)/(young_x*young_y*delta)*(ONE-d3)**2
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
        dsdeEl_d(5,5) = g_yz*(ONE - d1)*(ONE - d3)
        dsdeEl_d(5,6) = 0
        dsdeEl_d(6,1) = 0
        dsdeEl_d(6,2) = 0 
        dsdeEl_d(6,3) = 0
        dsdeEl_d(6,4) = 0
        dsdeEl_d(6,5) = 0
        dsdeEl_d(6,6) = g_xz*(ONE - d3)*(ONE - d2)
         
        
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
          d_C_d_d1(5,5) = g_yz*(d3 -ONE)
          
          
          
          do i=1,mcomp
              do j=1,mcomp
                  C_T_1_a(i) = C_T_1_a(i) + (d_C_d_d1(i,j)*T_strain(j)) 
              end do
          end do 
          
          
          if (T_strain(1)> 0)  then
          
            C_T_1_b  = (/ ((ONE - k1*F_f)/(F_f**3 * eps_11_f_t**2))
     &     *exp(k1*(F_f - ONE))*T_strain(1), ZERO, ZERO,((ONE - k1*F_f)/ 
     &     (F_f**3 * eps_12_f**2))*exp(k1*(F_f -ONE))*T_strain(4), ZERO,     
     &     ((ONE - k1*F_f)/(F_f**3 * eps_13_f**2))*exp(k1*(F_f - ONE))
     &     *T_strain(6)/)
            
          
          else
              
            C_T_1_b  = (/ ((ONE - k1*F_f)/(F_f**2 * eps_11_f_c))
     &     *exp(k1*(F_f - ONE)),ZERO,ZERO,ZERO,ZERO,ZERO/)
            
            
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
          d_C_d_d2(6,6) = g_xz*(d3 - ONE)
          
          
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
              
              C_T_2_b=(/ZERO,(HALF)*((ONE-k2*sqrt(F_m))/(sqrt(F_m)**3))*
     &        exp(k2*(sqrt(F_m) -ONE))*term1, (HALF)*((1 - k2*sqrt(F_m))
     &        /(sqrt(F_m)**3))
     &        *exp(k2*(sqrt(F_m) - ONE))*term2, ((ONE - k2*sqrt(F_m))/
     &        (sqrt(F_m)**3 *
     &        eps_12_f**2))*exp(k2*(sqrt(F_m) - ONE))*T_strain(4), 
     &        ((ONE - k2*sqrt(F_m))/(sqrt(F_m)**3*eps_23_f**2))*
     &        exp(k2*(sqrt(F_m)-1))*T_strain(5),    
     &        ((ONE - k2*sqrt(F_m))/(sqrt(F_m)**3 * eps_13_f**2))*
     &        exp(k2*(sqrt(F_m) - ONE)) 
     &        *T_strain(6) /)         
     

      
          else 
                            
              term1 = ((TWO*(T_strain(2) + T_strain(3))/(eps_22_f_c
     &    *eps_33_f_c))  +  ((eps_22_f_c/(2*eps_12_f)) - ONE)/eps_22_f_c  
     &     -  (T_strain(3)/eps_23_f**2))/(2*F_m)
              
              term2 = ((TWO*(T_strain(2) + T_strain(3))/(eps_22_f_c
     &    *eps_33_f_c))  +  ((eps_22_f_c/(2*eps_12_f)) - ONE)/eps_22_f_c  
     &     -  (T_strain(2)/eps_23_f**2))/(2*F_m) 
              
              
              C_T_2_b  = (/ZERO, ((ONE - k2*sqrt(F_m))/(sqrt(F_m)**2))*
     &        exp(k2*(sqrt(F_m) - ONE))*term1, ((ONE - k2*sqrt(F_m))
     &        /(sqrt(F_m)**2))
     &        *exp(k2*(sqrt(F_m) - ONE))*term2,  ((ONE - k2*sqrt(F_m))
     &        /(sqrt(F_m)**3 * 
     &        eps_12_f**2))*exp(k2*(sqrt(F_m) - ONE))*T_strain(4), 
     &        ((ONE - k2*sqrt(F_m)) 
     &        /(sqrt(F_m)**3 * eps_23_f**2))*exp(k2*(sqrt(F_m) - ONE))
     &        *T_strain(5),
     &        ((ONE - k2*sqrt(F_m))/(sqrt(F_m)**3 * eps_13_f**2))*
     &        exp(k2*(sqrt(F_m) - ONE))
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
         d_C_d_d3(5,5)  =  g_yz*(d1 - ONE)     
         d_C_d_d3(6,6)  =  g_xz*(d2 - ONE)  
         
         
          
          do i=1,ncomp
              do j=1,ncomp
                  C_T_3_a(i) = C_T_3_a(i) + (d_C_d_d3(i,j)*T_strain(j)) 
              end do
          end do 
          
          
          if (T_strain(3) > 0) then
              
              C_T_3_b =(/ZERO,ZERO,((ONE-k3*F_z)/(F_z**3*eps_33_f_t**2))
     &    *exp(k3*(F_z - ONE))*T_strain(3),ZERO,((ONE - k3*F_z)/(F_z**3       
     &    * eps_23_f**2))*exp(k3*(F_z -ONE))*T_strain(5),((ONE - k3*F_z)         
     &    /(F_z**3 * eps_13_f**2))*exp(k3*(F_z - ONE))*T_strain(6) /)   
              
              
          else 
              
              C_T_3_b =(/ZERO,ZERO,((ONE-k3*F_z)/(F_z**3*eps_33_f_c**2))
     &    *exp(k3*(F_z - ONE))*T_strain(3),ZERO,((ONE - k3*F_z)/(F_z**3       
     &    * eps_23_f**2))*exp(k3*(F_z -ONE))*T_strain(5),((ONE - k3*F_z)         
     &    /(F_z**3 * eps_13_f**2))*exp(k3*(F_z - ONE))*T_strain(6) /)
              
          
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
      ustatev(4) = F_f
      ustatev(5) = F_m
      ustatev(6) = F_z
      

      

      
c *** Update stress i
      do i=1,ncomp
         stress(i) = sigma(i)
      end do
      

      
      if ((F_f**2 < 1) .AND. (F_m < 1) .AND. (F_z**2 < 1)) then
          
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



              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              

      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      