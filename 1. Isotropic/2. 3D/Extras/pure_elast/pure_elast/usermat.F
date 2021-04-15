*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements or plane elements
c       in plane strain or axisymmetric stress state. The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variable
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - statev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,l)            three time of shear moduli
c
c --- temperary variables for solution purpose
c      i, j
c      threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3)
c*************************************************************************
c
      keycut   = 0
      cutFactor = 0.d0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)

      twoG     = young / (ONE+posn)
      elast1= young*posn/((1.0D0+posn)*(1.0D0-TWO*posn))
      elast2= HALF*twoG
	  
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = elast2
c
c *** calculate elastic stiffness matrix (3d)
c
      dsdeEl(1,1)=(elast1+TWO*elast2)*G(1)*G(1)
      dsdeEl(1,2)=elast1*G(1)*G(2)+elast2*TWO*G(4)*G(4)
      dsdeEl(1,3)=elast1*G(1)*G(3)+elast2*TWO*G(5)*G(5)
      dsdeEl(1,4)=elast1*G(1)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(1,5)=elast1*G(1)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(1,6)=elast1*G(1)*G(6)+elast2*TWO*G(4)*G(5)
      dsdeEl(2,2)=(elast1+TWO*elast2)*G(2)*G(2)
      dsdeEl(2,3)=elast1*G(2)*G(3)+elast2*TWO*G(6)*G(6)
      dsdeEl(2,4)=elast1*G(2)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(2,5)=elast1*G(2)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(2,6)=elast1*G(2)*G(6)+elast2*TWO*G(2)*G(6)
      dsdeEl(3,3)=(elast1+TWO*elast2)*G(3)*G(3)
      dsdeEl(3,4)=elast1*G(3)*G(4)+elast2*TWO*G(5)*G(6)
      dsdeEl(3,5)=elast1*G(3)*G(5)+elast2*TWO*G(5)*G(3)
      dsdeEl(3,6)=elast1*G(3)*G(6)+elast2*TWO*G(6)*G(3)
      dsdeEl(4,4)=elast1*G(4)*G(4)+elast2*(G(1)*G(2)+G(4)*G(4))
      dsdeEl(4,5)=elast1*G(4)*G(5)+elast2*(G(1)*G(6)+G(5)*G(4))
      dsdeEl(4,6)=elast1*G(4)*G(6)+elast2*(G(4)*G(6)+G(5)*G(2))
      dsdeEl(5,5)=elast1*G(5)*G(5)+elast2*(G(1)*G(3)+G(5)*G(5))
      dsdeEl(5,6)=elast1*G(5)*G(6)+elast2*(G(4)*G(3)+G(5)*G(6))
      dsdeEl(6,6)=elast1*G(6)*G(6)+elast2*(G(2)*G(3)+G(6)*G(6))
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
c
c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)
c
c *** calculate the trial stress and
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         strainEl(i) = Strain(i) + dStrain(i) - epsPl(i)
      end do
      call vzero(sigElp, 6)
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i)+dsdeEl(j,i)*strainEl(j)
         end do
         sigElp(i) = sigElp(i) + sigi(i)
      end do
	  
c *** current yield stress
	  sigy     = sigy0 + dsigdep * pleq
	  
c *** Update stress for elastic/unloading -- There is no plastic strain here!!
      do i=1,ncomp
         stress(i) = sigElp(i)
      end do
	  
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
c
!      cutFactor = 0.5d0
      return
      end