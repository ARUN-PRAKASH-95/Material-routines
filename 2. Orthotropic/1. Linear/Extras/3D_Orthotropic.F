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


      DOUBLE PRECISION  dsdeEl(mcomp,mcomp), C4(3,3,3,3),
     &                  C(mcomp,mcomp), ii(mcomp), jj(mcomp),
     &                  sigma(mcomp), T_strain(mcomp)

      
	DOUBLE PRECISION    var0, var1, var2, var3, var4, var5,
     &                    var6, var7, var8
	
	
      INTEGER             i, j
      DOUBLE PRECISION    sigy0, twoG, young_x, young_y, 
     &                    young_z, pr_xy, pr_yz, pr_xz, 
     &                    g_xy, g_yz, g_xz, delta,
     &                    pr_yx, pr_zy, pr_zx,
     &                    xy_yx, yz_zy, zx_xz,
     &                    xyz, E_xyz
	
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
	g_xy     = 5000000000
	g_yz     = 3846200000
	g_xz     = 5000000000
	
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
	dsdeEl(4,4) = 2*g_yz
	dsdeEl(4,5) = 0
	dsdeEl(4,6) = 0
	dsdeEl(5,1) = 0
	dsdeEl(5,2) = 0 
	dsdeEl(5,3) = 0
	dsdeEl(5,4) = 0
	dsdeEl(5,5) = 2*g_xz
	dsdeEl(5,6) = 0
	dsdeEl(6,1) = 0
	dsdeEl(6,2) = 0 
	dsdeEl(6,3) = 0
	dsdeEl(6,4) = 0
	dsdeEl(6,5) = 0
	dsdeEl(6,6) = 2*g_xy

	
      
c *** Calculate current strain
      do i=1,ncomp
         T_strain(i) =  Strain(i) + dStrain(i) 
      end do
  
      
	call fourth_order(dsdeEl,C4)
	
	
**********Restore stiffness tensor as 6*6 matrix***********
         do i = 1,mcomp
           do j = 1,mcomp
               C(i,j) = C4(ii(i),jj(i),ii(j),jj(j))
           end do
	   end do  
	

	call vzero(sigma, 6) 	   
***** Calculate the current stress ********	
      
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(i,j) = C(i,j)
            sigma(i) = sigma(i) + C(i,j)*T_strain(j)
         end do
	end do

	
	
c *** Update stress 
      do i=1,ncomp
         stress(i) = sigma(i)
      end do	
	
      return
	end
	
	
	
	
	
	
c************  Subroutine for 6*6 elastic stiffnes matrix into 4th order tensor   **********      
	subroutine fourth_order(ST1,FT1)

	  double precision ST1(6,6), nn(3,3)
        double precision FT1(3,3,3,3)
        integer i,j,k,l,n,m
        n = 3
	  m = 6
	  
	  
 
	  
	  nn   = transpose(reshape((/1,4,6,4,2,5,6,5,3/),shape(nn)))
	  

	  
	  
	  
         do i = 1,n
           do j = 1,n
              do k=1,n
                  do l=1,n
                   FT1(i,j,k,l) = ST1(nn(i,j),nn(k,l))
                  end do
              end do
           end do
	   end do
	  
      return
      end