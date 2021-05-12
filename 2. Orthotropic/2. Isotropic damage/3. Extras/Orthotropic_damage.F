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
	
	
	
      integer :: mcomp
      DOUBLE PRECISION    HALF, ONE, TWO, ZERO
      
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 mcomp      = 6
     &                 )
	

      EXTERNAL          vzero, vmove, get_ElmData,
     &                  fourth_order

********* USER DEFINED VARIABLES (FOR CALCULATION PURPOSE) **********	


      DOUBLE PRECISION  dsdeEl(mcomp,mcomp), C(3,3,3,3),
     &                  ii(mcomp), jj(mcomp), a(3,3,3,3),
     &                  Strain_tensor(3,3),b2(3,3),b(3,3),
     &                  c_three(3,3),d(3,3,3,3), C_T(3,3,3,3)

      
	DOUBLE PRECISION    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	
	
      INTEGER             i, j,k,l,n,P
      DOUBLE PRECISION    sigy0, twoG, young_x, young_y, 
     &                    young_z, pr_xy, pr_yz, pr_xz, 
     &                    g_xy, g_yz, g_xz, delta,
     &                    pr_yx, pr_zy, pr_zx,
     &                    xy_yx, yz_zy, zx_xz,
     &                    xyz, E_xyz,damage,eta,
     &                    damage_new,yield,epsilon_f,b1
 
	

	
	
	
	integer :: iott,wrinqr
      external wrinqr
c
c*************************************************************************
c
      iott = wrinqr(WR_OUTPUT)
      write(iott,*) ' ************************************************'
      write(iott,*) ' * #DEBUG# UserMatLib.dll USERMAT1D             *'
      write(iott,*) ' ************************************************'

	
      ii = (/1,2,3,1,2,1/)
      jj = (/1,2,3,2,3,3/)

	
	
      keycut   = 0
	young_x  = prop(1)
	young_y  = prop(2)
	young_z  = prop(3)
	pr_xy    = prop(4)
	pr_yz    = prop(5)
	pr_xz    = prop(6)
	g_xy     = 4700e6
	g_yz     = 4700e6
	g_xz     = 4700e6
	
	
      yield    = 1250e6
	n        = 3 
      damage   = ustatev(1)
	eta      = 1e-5
	
	pr_yx = (young_y * pr_xy) / young_x
	pr_zy = (young_z * pr_yz) / young_y
	pr_zx = (young_z * pr_xz) / young_x
	
	
	xy_yx = pr_xy*pr_yx
	yz_zy = pr_yz*pr_zy 
	zx_xz = pr_zx*pr_xz
	xyz   = TWO*pr_xy*pr_yz*pr_zx
	E_xyz = young_x*young_y*young_z
	
	delta = (ONE - (xy_yx) - (yz_zy) - (zx_xz) - (xyz)) / E_xyz
	
	
	
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
	dsdeEl(4,4) = TWO * g_yz
	dsdeEl(4,5) = 0
	dsdeEl(4,6) = 0
	dsdeEl(5,1) = 0
	dsdeEl(5,2) = 0 
	dsdeEl(5,3) = 0
	dsdeEl(5,4) = 0
	dsdeEl(5,5) = TWO * g_xz
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = TWO * g_xy

	
	
c**********  Calculate the failure strain ************         
         epsilon_f = yield / dsdeEl(1,1)
         P = 100000

      

c *** Calculate current strain
      do i=1,ncomp
         Strain(i) =  Strain(i) + dStrain(i) 
	end do
	
	
	
C********  Restore strain tensor ***********
      Strain_tensor = transpose(reshape((/ Strain(1), Strain(4)/TWO, 
     &     Strain(6)/TWO, Strain(4)/TWO, Strain(2), Strain(5)/TWO,
     &                Strain(6)/TWO	, Strain(5)/TWO, Strain(3) /), 
     &                shape(Strain_tensor)))

	
	

	
***** Calculate the current stress ********	
      if (Strain(1) < epsilon_f) then
          
          do i=1,ncomp
             do j=1,ncomp
                stress(i) = stress(i) + (1+eta-damage)
     &            *dsdeEl(j,i)*dStrain(j)
             end do
          end do
          
      else
          
          damage_new = 1  - exp(-P*(Strain(1)-epsilon_f))
          
          if (damage_new >= damage) then
              damage = damage_new
          else
              damage = damage
          endif
          
          
          
          do i=1,ncomp
             do j=1,ncomp
                stress(i) = stress(i) + (1+eta-damage)
     &            *dsdeEl(j,i)*dStrain(j)
             end do
		end do          
    

c***************  Tangent stiffness *************		

	call fourth_order(dsdeEl,C)
	

c*************** First term of the tangent stiffness (1-d+eta)*C **********          
          do i = 1,n
		       do j = 1,n
			       do k=1,n
				       do l=1,n
                            a(i,j,k,l) = (1 + eta- damage)*C(i,j,k,l) 
				       end do
                    end do
		       end do
		end do   
		

		
		
c*************** Second term of the tangent stiffness ********
          b1 = P*exp(-P*(Strain(1)-epsilon_f))
      
          b2 = transpose(reshape((/ ONE, ZERO, ZERO, ZERO, ZERO, ZERO, 
     &    ZERO,ZERO,ZERO /),shape(b2)))  
                  
      
          do i = 1,3
              do j = 1,3
                  b(i,j) = b1*b2(i,j)
              end do
		end do
		

		
c****  Double contraction of 4th order elastic stiffness and 2nd order strain *******
          call t4_contr_t2(C,Strain_tensor,c_three) 

c****  Dyadic product of two 2nd order tensors ******
          call t2_dyad_t2(ONE,c_three,b,d)

c**** Tangent stiffness *****
          do i = 1,n
		       do j = 1,n
			       do k=1,n
				       do l=1,n
                            C_T(i,j,k,l) = a(i,j,k,l) - d(i,j,k,l)
				       end do
                    end do
		       end do
          end do 

      endif
      
      ustatev(1) = damage
	
	
	
**********Restore stiffness tensor as 6*6 matrix***********
       if (Strain(1) < epsilon_f) then      
          do i = 1,mcomp
               do j = 1,mcomp
                   dsdePl(i,j) = (1+eta - damage)*dsdeEl(i,j)
  
               end do
          end do    
       else 
          do i = 1,mcomp
               do j = 1,mcomp
                   dsdePl(i,j) = C_T(ii(i),jj(i),ii(j),jj(j))
               end do
          end do
      endif    
      


c
      return
	end
	
	
	
	
	
c************  Subroutine for dyadic product of two second order tensors   **********      

      
      subroutine t2_dyad_t2(s,ST1,ST2,FT1)
        double precision s
        double precision ST1(3,3)
	  double precision ST2(3,3)
        double precision FT1(3,3,3,3)
        integer i,j,k,l,n
        n = 3
        


	   
         do i = 1,n
             do j = 1,n
			   do k=1,n
				   do l=1,n
						FT1(i,j,k,l) = s*ST1(i,j)*ST2(k,l)
				   end do
                end do
		   end do
	   end do
	  
	  return
	end
	
	

	
c************  Subroutine for dyadic product of 4th and 2nd tensors   **********      
	subroutine t4_contr_t2(C,Strain_tensor,c_three)
        
        double precision C(3,3,3,3)
        double precision Strain_tensor(3,3)
	  double precision c_three(3,3)
        integer i,j,k,l,n
        n = 3
        


	   
         do i = 1,n
             do j = 1,n
			   do k=1,n
				   do l=1,n
						c_three(i,j) = C(i,j,k,l) * Strain_tensor(k,l)
				   end do
                end do
		   end do
	   end do
	  
	  return
	end	
	
	
	
c************  Subroutine for 6*6 elastic stiffnes matrix into 4th order tensor   **********      
	subroutine fourth_order(ST3,FT2)

	  double precision ST3(6,6), nn(3,3)
        double precision FT2(3,3,3,3)
        integer i,j,k,l,n,m
        n = 3
	  m = 6
	  
	  
 
	  
	  nn   = transpose(reshape((/1,4,6,4,2,5,6,5,3/),shape(nn)))
	  

	  
	  
	  
         do i = 1,n
           do j = 1,n
              do k=1,n
                  do l=1,n
                   FT2(i,j,k,l) = ST3(nn(i,j),nn(k,l))
                  end do
              end do
           end do
	   end do
	  
      return
      end