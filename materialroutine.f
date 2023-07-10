       subroutine rmap(params,F,STATEV,dt,sig)
       implicit none 
       real*8 muA,lamlA,F(3,3),B(3,3),J,sigA(3,3), STATEV(32),
     +        muB,kapB,m,dgam0,Finv(3,3),alpha,
     +        FiBOld(3,3),FiBOldinv(3,3),ID(3,3),dgamb1,dgamb0,
     +        tauhat,FiB0inv(3,3),FeB0(3,3),BeB0(3,3),
     +        sigB0(3,3),sigBvm0,NiBold(3,3),f0,tau,FiB1(3,3),dt,
     +        FiB1inv(3,3),FeB1(3,3),BeB1(3,3),sigB1(3,3),det,
     +        sigBvm1,f1,NiB1(3,3),dgamB,tol,sigB(3,3),res,
     +        tau0,tauss,h,tauOld,tr,dgambOld,NiB0(3,3),tm1(3,3),
     +        tm2(3,3),tm3(3,3),sigTot(3,3),sigHtot,
     +        params(32),sig(3,3), array(10)
       integer i, maxit
       real*8 one,three,third
       parameter(one=1.,three=3.,third=one/three)
C      Identity matrix
       ID(1,1)=1.
       ID(2,2)=1.
       ID(3,3)=1.
       ID(1,2)=0.
       ID(2,1)=0.
       ID(1,3)=0.
       ID(3,1)=0.
       ID(2,3)=0.
       ID(3,2)=0.
C======================================================================
C           Material params
C======================================================================
      muA    = params(1)         !Initial shear modulus part A 
      lamLA  = params(2)         !Locking stretch part A
      muB    = params(3)         !Initial shear modulus part B
      kapB   = params(4)         !Bulk modulus part B
      alpha  = params(5)         !Pressure sensitivity parameter
      dgam0  = params(6)         !Reference shear strain rate
      m      = params(7)         !Strain rate sensitivity param
      tau0   = params(8)         !Initial shear strength
      tauss  = params(9)         !Steady state shear strength
      h      = params(10)        !Softening modulus
      maxit  = int(params(12))   !Maximum number of iterations
      tol    = params(13)        !Convergence criteria

C######################################################################
C
C         PART A
C
C######################################################################

C     Caculating the jacobian, det(F)
      J = det(F)

C     Calculating the  isochoric left Cauchy-Green def.grad.
C     B = J^(-2/3)FF^T
C     Transposing F and storing in temp matrix
      call mattrans(F,tm1)
      call matmult(F,tm1,B)
C     Calculating Cauchy stress from the Eight-Chain model
      call ECStressA(muA,lamLA,B,J,sigA)
C######################################################################
C
C         PART B
C
C######################################################################

C     Retrieving and calulating variables at time tn used in the 
C     semi-implicit update scheme, but not updated in the iterations

C     Old inelastic deformation gradient 
C     initialising        
      if(STATEV(1) .eq. 0.0 ) then
            FiBold = ID
      else
            FiBold(1,1) = STATEV(1)
            FiBold(2,2) = STATEV(2)
            FiBold(3,3) = STATEV(3)
            FiBold(1,2) = STATEV(4)
            FiBold(2,3) = STATEV(5)
            FiBold(3,1) = STATEV(6)
            FiBold(2,1) = STATEV(7)
            FiBold(3,2) = STATEV(8)
            FiBold(1,3) = STATEV(9)
      end if

C     Inverting some deformation gradients
      call matinv(FiBold,FiBoldinv)
      call matinv(F,Finv)

      NiBold(1,1) = STATEV(10)
      NiBold(2,2) = STATEV(11)
      NiBold(3,3) = STATEV(12)
      NiBold(1,2) = STATEV(13)
      NiBold(2,3) = STATEV(14)
      NiBold(3,1) = STATEV(15)
      NiBold(2,1) = STATEV(16)
      NiBold(3,2) = STATEV(17)
      NiBold(1,3) = STATEV(18)
        
      if(STATEV(22).eq.0.0) then
            tauOld = tau0
      else
            tauOld = STATEV(22)
      end if


C=======================================================================
C
C        Quasi-implicit iteration scheme 
C
C=======================================================================
      do i= 1,maxit

C     Initial guesses
      if(i .eq. 1) then
            if(STATEV(23) .eq. 0.) then
                  dgamBOld = 1.d-12
            else
                  dgamBOld = 1.*STATEV(23)
            end if
      end if

      if(i .eq. 1) then
C=======================================================================
C        Increment 0, Starting the iteration scheme
C=======================================================================
      dgamB0 = dgamBOld

C     New inelastic deformation gradient
C     FiB1inv= (I-dgamb0*dt*F^-1 NiB0 F^-1) * FiBoldinv 
      call matmult(Finv,NiBold,tm1)
      call matmult(tm1,F,tm2)
      tm3 = ID-dgamb0*dt*tm2

      call matmult(tm3,FiBoldinv,FiB0inv)

C     Calculating the elastic deformation gradient FeB=F*FiBinv
      call matmult(F,FiB0inv,FeB0)

C     Transposing FeB1 and storing in temp matrix
      call mattrans(FeB0,tm1)
C     Calculating the right Cauchy Green def.grad B=FeB*FeB^T
      call matmult(FeB0,tm1,BeB0)
C     Calculating Cauchy stress from the neo-Hookean model
      call NHStressB(muB,kapB,BeB0,J,sigB0)

      tau = (tauOld+dt*dgamB0*h)/(1.+dt*dgamB0*h/tauss)
      
      call sigeq(sigB0,sigBvm0)

C     Total hydrostatic stress
      sigtot = sigA+sigB0
      sightot = tr(sigtot)*third

      call tauSub(tau,sightot,alpha,tauHat)

C     Residual function
      f0 = dgamb0-dgam0*(sigBvm0/tauhat)**m
C     Direction of inelastic flow
      call PlasticFlowDir(sigB0,NiB0)
C     First estimate of dgamb1
      dgamb1 = dgam0*(sigBvm0/tauhat)**m
      if(dgamb1 .gt. 10.*dgamb0) then
            dgamb1 = 10*dgamb0
      end if

      end if

C=======================================================================
C        Increment I
C=======================================================================
C     New inelastic deformation gradient
C     FiB1= (I-dgamb1*dt*NiB0)^-1 * FiBold 
      call matmult(Finv,NiB0,tm1)
      call matmult(tm1,F,tm2)
      tm3 = ID-dgamb1*dt*tm2

      call matmult(tm3,FiBoldinv,FiB1inv)

C     Calculating the elastic deformation gradient FeB=F*FiBinv
      call matmult(F,FiB1inv,FeB1)

C     Transposing FeB1 and storing in temp matrix
      call mattrans(FeB1,tm1)
C     Calculating the right Cauchy Green def.grad B=FeB*FeB^T
      call matmult(FeB1,tm1,BeB1)
C     Calculating Cauchy stress from the neo-Hookean model
      call NHStressB(muB,kapB,BeB1,J,sigB1)
C     Direction of inelastic flow
      call PlasticFlowDir(sigB1,NiB1)

      tau = (tauOld+dt*dgamB1*h)/(1.+dt*dgamB1*h/tauss)

      call sigeq(sigB1,sigBvm1)

C     Total hydrostatic stress
      sigtot = sigA+sigB1
      sightot = tr(sigtot)*third 

      call tauSub(tau,sightot,alpha,tauHat)
C     Residual functions
      f1 = dgamb1-dgam0*(sigBvm1/tauHat)**m

C=============================================================================
C     Update effective strain rate, check for convergence
C=============================================================================

C     Secant update of effective shear strain rate
      dgamB = dgamB1-(dgamB1-dgamB0)/(f1-f0)*f1  
      if(dgamb .gt. 10.*dgamb1) then
            dgamb = 10*dgamb1
      else if(dgamb .lt. 0.1*dgamb1) then
            dgamb = 0.1*dgamb1
      end if
      res = abs(f1)
      array(i)=dgamb1

      if((res .le. tol) .and. (i .ge. 2)) then
            STATEV(23) = dgamB1
            STATEV(22) = tau
            STATEV(31) = 1
        exit 
      end if    

C     Updating eff. shear strain rates for next iteration
      dgamB0 = max(dgamB1,0.0)
      dgamB1 = max(dgamB,0.0)
      f0 = f1
      NiB0 = NiB1

      if(i .eq. maxit) then
            write(*,*) 'MATERIAL ROUTINE DID NOT CONVERGE'
            write(*,*) 'J=', J
            write(*,*) 'damage=', STATEV(30)
            write(*,*) 'sigh=', sightot
            write(*,*) 'sigbvm0=',sigbvm0
            write(*,*) 'tauHat=',tauHat
            write(*,*) 'dgambold=', dgambold
            write(*,*) 'dgamb=',array
            write(*,*) 'sigA=',sigA
            write(*,*) 'sigB=',sigB1
            write(*,*) 'F11=', F(1,1)
            write(*,*) 'F22=', F(2,2)
            write(*,*) 'F33=', F(3,3)
            write(*,*) 'F12=', F(1,2)
            write(*,*) 'F23=', F(2,3)
            write(*,*) 'F13=', F(1,3)
            write(*,*) 'F21=', F(2,1)
            write(*,*) 'F32=', F(3,2)
            write(*,*) 'F31=', F(3,1)
            STATEV(31) = 0 
            CALL XPLB_EXIT 
      end if
      end do 

      sigB = sigB1
      sig = sigA + sigB


      call matinv(FiB1inv,FiB1)

C======================================================================
C       Storing history variables
C======================================================================
C     New inverse inelastic deformation gradient
      STATEV(1) = FiB1(1,1)
      STATEV(2) = FiB1(2,2)
      STATEV(3) = FiB1(3,3)
      STATEV(4) = FiB1(1,2)
      STATEV(5) = FiB1(2,3)
      STATEV(6) = FiB1(3,1)
      STATEV(7) = FiB1(2,1)
      STATEV(8) = FiB1(3,2)
      STATEV(9) = FiB1(1,3)

C     New inelastic flow direction
      STATEV(10) = NiB1(1,1)
      STATEV(11) = NiB1(2,2)
      STATEV(12) = NiB1(3,3)
      STATEV(13) = NiB1(1,2)
      STATEV(14) = NiB1(2,3)
      STATEV(15) = NiB1(3,1)
      STATEV(16) = NiB1(2,1)
      STATEV(17) = NiB1(3,2)
      STATEV(18) = NiB1(1,3)

      STATEV(19) = i

      STATEV(20) = res
      STATEV(21) = tauhat

      end
C#######################################################################
C                                                                      #
C                             FUNCTIONS                                #
C                                                                      #
C#######################################################################
c
C ----------------------------------------------------------------------
C  det: Calculate determinant of 3x3 matrix
C  Input: A (3x3 matrix)
C ----------------------------------------------------------------------
      function det(A)
          implicit none
          real*8 det,A,zero
          dimension A(3,3)
          parameter(zero=0.d0)
c
          det = A(1,1)*A(2,2)*A(3,3)
     1        - A(1,2)*A(2,1)*A(3,3)
     2        - A(1,1)*A(2,3)*A(3,2)
     3        + A(1,3)*A(2,1)*A(3,2)
     4        + A(1,2)*A(2,3)*A(3,1)
     5        - A(1,3)*A(2,2)*A(3,1)
          return
      end
C ----------------------------------------------------------------------
C  kdel: Kronecker delta
C        1 if m=n
C        0 else
C ----------------------------------------------------------------------
      function kdel(m,n)
          implicit none
          integer m,n
          real*8 kdel, one, zero
          parameter(one=1.d0,zero=0.d0)
c
          if (m.eq.n) then
              kdel = one
          else
              kdel = zero
          endif
          return
      end
C ----------------------------------------------------------------------
C  tr: Calculate trace of matrix
C  Input: stress (3x3 matrix)
C  Output: stress_kk
C ----------------------------------------------------------------------
      function tr(stress)
          implicit none
          real*8 stress, tr
          dimension stress(3,3)
c
          tr = stress(1,1)+stress(2,2)+stress(3,3)
          return
      end
C ----------------------------------------------------------------------
C  dev_sub: Calculate deviatoric tensor
C  Input: stress (3x3 matrix)
C  Output: dstress (3x3 matrix)
C ----------------------------------------------------------------------
      subroutine dev_sub(stress,dstress)
          implicit none
          real*8 tr, one, three, third, dstress, stress
          dimension dstress(3,3), stress(3,3)
          parameter(one=1., three=3., third=one/three)
c
          dstress(1,1) = stress(1,1)-tr(stress)*third
          dstress(2,2) = stress(2,2)-tr(stress)*third
          dstress(3,3) = stress(3,3)-tr(stress)*third
          dstress(1,2) = stress(1,2)
          dstress(2,1) = stress(2,1)
          dstress(1,3) = stress(1,3)
          dstress(3,1) = stress(3,1)
          dstress(2,3) = stress(2,3)
          dstress(3,2) = stress(3,2)
      end
C ----------------------------------------------------------------------
C  invl: Calculate inverse langevin
C  Input: Scalar, x
C  Output: invl
C ----------------------------------------------------------------------
      function invl(x)
          implicit none
          real*8 invl,x,y
          y = min(x,0.999)
          invl = y*(3.d0-2.6d0*y+0.7d0*y**2)/
     1              ((1.d0-y)*(1.d0+0.1d0*y))
          return
      end
c
C ----------------------------------------------------------------------
C  ECStress: Compute stress from Eight-Chain potential in part A
C  Input: mu - Shear modulus
C       lamL - Locking stretch
C         K  - Bulk modulus
C         Be - Elastic left Cauchy-Green deformation tensor
C         Je - Elastic jacobi
C  Output: stress
C ----------------------------------------------------------------------
      subroutine ECStressA(mu,Laml,Be,Je,stress)
          implicit none
          real*8 mu,Be,Je,tr,kdel,stress,Laml,invl,trBe,lambdaEc,two,three,
     1           one, Biso, frac, third
          dimension Be(3,3), stress(3,3), Biso(3,3)
          integer i,j
          parameter(two=2.d0,three=3.d0,one=1.d0,third=one/three)
          Biso = Je**(-two*third)*Be
          trBe = tr(Biso)
          lambdaEc = sqrt(trBe*third)
          frac = Laml*mu/(three*Je*lambdaEc)*invl(lambdaEc/Laml)
          do i=1,3
              do j=1,3
              stress(i,j) = frac*(Biso(i,j)-lambdaEc**2*kdel(i,j))
              enddo
          enddo
      end
C ----------------------------------------------------------------------
C  ECStressA: Compute stress from neo-Hookean potential in part B
C  Input: mu - Shear modulus
C         K  - Bulk modulus
C         Be - Elastic left Cauchy-Green deformation tensor
C         Je - Elastic jacobi
C  Output: stress
C ----------------------------------------------------------------------
      subroutine NHStressB(mu,K1,Be,Je,stress)
          implicit none
          real*8 mu,Be,Je,tr,kdel,stress,trBe,two,three,
     1           one,K1,Biso,third,frac
          dimension Be(3,3), stress(3,3),Biso(3,3)
          integer i,j
          parameter(two=2.d0,three=3.d0,one=1.d0,third=one/three)
          Biso = Je**(-two*third)*Be
          trBe = tr(Biso)
          frac = mu/Je
          do i=1,3
              do j=1,3
              stress(i,j) = frac*(Biso(i,j)-trBe*third*kdel(i,j))
     1                     +K1*(Je-1.)*kdel(i,j)
              enddo
          enddo
      end
C ----------------------------------------------------------------------
C  ddot: Calculate double dot product
C  Input: t1 (3x3) and t2 (3x3)
C  Output: t1_ij*t1_ij
C ----------------------------------------------------------------------
      function ddot(t1,t2)
          implicit none
          real*8 t1,t2,ddot,zero
          dimension t1(3,3),t2(3,3)
          integer i,j
          parameter(zero=0.d0)
          ddot = zero
          do i=1,3
              do j=1,3
                  ddot=ddot+t1(i,j)*t2(i,j)
              enddo
          enddo
      end
C ----------------------------------------------------------------------
C  matmult: Matrix multiplaction
C  Input: a,b (3x3 matrices)
C  Output: c (3x3 matrix)
C ----------------------------------------------------------------------
      subroutine matmult(a,b,c)
          implicit none
          real*8 a,b,c,tmp,zero
          dimension a(3,3),b(3,3),c(3,3)
          integer i,j,k
          parameter(zero=0.d0)
c
          tmp = zero
          do i=1,3
              do j=1,3
                  do k=1,3
                      tmp = tmp+a(i,k)*b(k,j)
                  enddo
                  c(i,j) = tmp
                  tmp = zero
              enddo
          enddo
      end
C ----------------------------------------------------------------------
C matinv: compute inverse of 3x3 matrix
C input: A (3x3)
C output: B (3x3)
C ----------------------------------------------------------------------
      subroutine matinv(a,b)
        implicit none
        real*8 a(3,3), b(3,3), det
        b(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
        b(1,2) = a(1,3)*a(3,2)-a(1,2)*a(3,3)
        b(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
        b(2,1) = a(2,3)*a(3,1)-a(2,1)*a(3,3)
        b(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
        b(2,3) = a(1,3)*a(2,1)-a(1,1)*a(2,3)
        b(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
        b(3,2) = a(1,2)*a(3,1)-a(1,1)*a(3,2)
        b(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)
        b = b/det(a)
      end
C ----------------------------------------------------------------------
C mattrans: transpose 3x3 matrix
C input: A (3x3)
C output: B (3x3)
C ----------------------------------------------------------------------
      subroutine mattrans(a,b)
          implicit none
          real*8 a(3,3), b(3,3)
          integer i,j
          do i=1,3
              do j=1,3
                  b(i,j) = a(j,i)
              enddo
          enddo
      end
c
C ----------------------------------------------------------------------
C PlasticFlowDir: Direction of plastic flow
C Input:  beta     : parameter (amount of plastic vol. strain)
C         sigma    : Cauchy stress tensor
C Output: Ngrad    : Direction of plastic flow
C ----------------------------------------------------------------------
c
      subroutine PlasticFlowDir(sigma,Ngrad)
          implicit none
          real*8 sigma, Ngrad, devsigma, ID, J2, ddot, one, two, three,
     1           zero,half
          dimension sigma(3,3), Ngrad(3,3), ID(3,3), devsigma(3,3)
          parameter(one=1.d0,two=2.d0,three=3.d0,zero=0.d0,half=0.5)
          ! Construct identity matrix
          ID(1,1) = one
          ID(2,2) = one
          ID(3,3) = one
          ID(1,2) = zero
          ID(2,1) = zero
          ID(3,1) = zero
          ID(1,3) = zero
          ID(2,3) = zero
          ID(3,2) = zero
          ! Calculate deviatoric Cauchy stress and second deviatoric invariant
          call dev_sub(sigma,devsigma)
          J2 = ddot(devsigma,devsigma)*half
c
          if (J2.eq.zero) then
            Ngrad = zero*ID
          else
          ! Calculate direction of plastic flow
            Ngrad = three*half*devsigma/sqrt(three*J2) 
          endif
      end
c
C ----------------------------------------------------------------------
C sigeq: von Mises quivalent stress
C Input: stress   : 3x3 stress tensor
C Output: sigeq
C ----------------------------------------------------------------------
c
      subroutine sigeq(stress,eqstress)
          implicit none
          real*8 eqstress,stress,two,three,J2,ddot,devstress
          dimension stress(3,3),devstress(3,3)
          parameter(two=2.d0,three=3.d0)
c
          call dev_sub(stress,devstress)
          J2 = ddot(devstress,devstress)*0.5
          eqstress = sqrt(three*J2)
      end
C ----------------------------------------------------------------------
C tauHat: Pressure modified shear strength
C Input:  tau, sigH, alpha
C Output: tauHat
C ----------------------------------------------------------------------
c
      subroutine tauSub(tau,sigH,alpha,tauHat)
          implicit none
          real*8 tau,sigH,alpha,tauHat
c
          tauHat = max(tau-alpha*sigH,1.)
      end