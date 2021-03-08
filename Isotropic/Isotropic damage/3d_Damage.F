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
	
      EXTERNAL          vzero, vmove, get_ElmData,
     &                  t2_dyad_t2, getP4sym,
     &                  t4_contr_t2
   
      double precision  dsdeEl(mcomp,mcomp), sigma(mcomp),
     &                  sigi (mcomp), strains(mcomp),
     &                  id(3,3), c2(3,3,3,3), c1(3,3,3,3),
     &                  C(3,3,3,3), ii(mcomp), jj(mcomp),
     &                  c3(3,3), Strain_tensor(3,3),
     &                  a(3,3,3,3),b2(3,3),b(3,3),
     &                  d(3,3,3,3), C_T(3,3,3,3)
                  

      
      integer             i,j,k,l,n,P
      double precision    young, posn, mu,  
     &                    twoG, lambda,
     &                    xK, yield,epsilon_f,damage,
     &                    damage_new,b1
   
      

      
      
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
      ii = (/1,2,3,1,2,1/)
      jj = (/1,2,3,2,3,3/)

 
      
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      twoG     = young / (ONE+posn)
	yield    = 200
	n        = 3 
	
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG


********   LAME'S CONSTANTS ********
      
      lambda   =  young*posn / ((ONE + posn)*(1 - TWO*posn)) 
      mu       =  young / (TWO*(1+posn))
	xK    =  young / (THREE*(1 - TWO*posn))

      
      
c****** First term of the stiffness tensor(c1)
      call getP4sym(mu,id,c1)
      
      
      
c****** Second term of the stiffness tensor(c2)
      call t2_dyad_t2(xK,id,id,c2)
      
      
********* ELASTIC STIFFNESS TENSOR *********
        
      
      do i = 1,n
		   do j = 1,n
			   do k=1,n
				   do l=1,n
                        C(i,j,k,l) = c1(i,j,k,l) + c2(i,j,k,l)
				   end do
                end do
		   end do
      end do 
      
      
**********Restore stiffness tensor as 6*6 matrix***********
         do i = 1,mcomp
           do j = 1,mcomp
               dsdeEl(i,j) = C(ii(i),jj(i),ii(j),jj(j))
           end do
         end do  
	   

c**********  Calculate the failure strain ************         
         epsilon_f = yield / dsdeEl(1,1)
         P = 2500

      

c *** Calculate current strain
      do i=1,ncomp
         strains(i) = Strain(i) + dStrain(i) 
	end do
      
      Strain_tensor = transpose(reshape((/ strains(1), strains(4)/TWO, 
     & strains(6)/TWO, strains(4)/TWO, strains(2), strains(5)/TWO,
     & strains(6)/TWO	, strains(5)/TWO, strains(3) /), 
     & shape(Strain_tensor)))

      
      
      
***** Calculate the current stress ********	
      if (strains(1) < epsilon_f) then
          
          do i=1,ncomp
             do j=1,ncomp
                stress(i) = stress(i) + dsdeEl(j,i)*dStrain(j)
             end do
          end do
          
      else
          
          damage_new = 1  - exp(-P*(strains(1)-epsilon_f))
          
          if (damage_new >= damage) then
              damage = damage_new
          else
              damage = damage
          endif
      
          do i=1,ncomp
             do j=1,ncomp
                stress(i) = stress(i)+(1-damage)*dsdeEl(j,i)*dStrain(j)
             end do
          end do          
          
c*************** First term of the tangent stiffness (1-d+eta)*C **********          
          do i = 1,n
		       do j = 1,n
			       do k=1,n
				       do l=1,n
                            a(i,j,k,l) = (1 - damage)*C(i,j,k,l) 
				       end do
                    end do
		       end do
          end do           
          

      
c*************** Second term of the tangent stiffness ********
          b1 = P*exp(-P*(strains(1)-epsilon_f))
      
          b2 = transpose(reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 0 /), 
     &         shape(b2)))      
      
          do i = 1,3
              do j = 1,3
                  b = b1*b2(i,j)
              end do
          end do

c****  Double contraction of 4th order elastic stiffness and 2nd order strain *******
          call t4_contr_t2(C,Strain_tensor,c3) 

          call t2_dyad_t2(ONE,b,c3,d)

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
      
**********Restore stiffness tensor as 6*6 matrix***********
       if (strains(1) < epsilon_f) then      
          do i = 1,mcomp
               do j = 1,mcomp
                   dsdePl(i,j) = C(ii(i),jj(i),ii(j),jj(j))
               end do
          end do    
       else 
          do i = 1,mcomp
               do j = 1,mcomp
                   dsdePl(i,j) = C_T(ii(i),jj(i),ii(j),jj(j))
               end do
          end do
      endif    
      


      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i))
      END DO
      sedEl    = sedEl * HALF

	
c
      return
      end
 
      

c************  Subroutine for dyadic product of two second order tensors   **********      
	subroutine t2_dyad_t2(s,A,B,C)
        double precision s
        double precision C(3,3,3,3)
        double precision A(3,3)
	  double precision B(3,3)
        integer i,j,k,l,n
        n = 3
        


	   
         do i = 1,n
             do j = 1,n
			   do k=1,n
				   do l=1,n
						C(i,j,k,l) = s*A(i,j)*B(k,l)
				   end do
                end do
		   end do
	   end do
	  
	  return
      end

      
c************  Subroutine for P4sym term  **********      
	subroutine getP4sym(mu,A,C)
        double precision mu
        double precision C(3,3,3,3),I4sym(3,3,3,3)
        double precision A(3,3)
        integer i,j,k,l,n
        n = 3
                 


        
        do i = 1,n
           do j = 1,n
               do k=1,n
                   do l=1,n
                      I4sym(i,j,k,l) = (0.5)*(A(i,k)*A(j,l)
     &                 +A(i,l)*A(j,k))
                  end do
              end do
          end do
        end do
       
        do i = 1,n
          do j = 1,n
              do k=1,n
                  do l=1,n
                      C(i,j,k,l) = 2*mu*(I4sym(i,j,k,l)-
     &                      (0.33)*(A(i,j)*A(k,l)))
                  end do
              end do
           end do
        end do      

	  
	  return
	  end


c************  Subroutine for dyadic product of 4th and 2nd tensors   **********      
	subroutine t4_contr_t2(C,A,B)
        
        double precision C(3,3,3,3)
        double precision A(3,3)
	  double precision B(3,3)
        integer i,j,k,l,n
        n = 3
        


	   
         do i = 1,n
             do j = 1,n
			   do k=1,n
				   do l=1,n
						B(i,j) = B(i,j) + C(i,j,k,l)*A(k,l)
				   end do
                end do
		   end do
	   end do
	  
	  return
      end	

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      