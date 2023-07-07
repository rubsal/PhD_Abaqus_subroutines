      INCLUDE './materialroutine.f'
C
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),jInfoArray(*)
C
      character*80 cmname

C DECLARE VARIABLES
      integer i,q,k,sflag,intPt  !Integers variables
      integer conFlag            !Convergence flag. 1=converged,0=not
      real*8 STATEVTmp(nstatev)  !temporary state variable vector
      real*8 muA,lamLA           !Elastic params part A
      real*8 muB,kapB            !Elastic params part B
      real*8 dgam0,m,tau         !Viscous params part B
      real*8 alpha               !Pressure sensitive parameters
      real*8 tau0,tauss,h        !Softening params part B
      real*8 epstr, depseq, epscr      
      real*8 rodInc(ndir,ndir)   !Incremental rate of def. tensor 
      real*8 epsdev(ndir,ndir)
      real*8 sigA(ndir,ndir)     !Cauchy stress tensor part A
      real*8 sigB(ndir,ndir)     !Cauchy stress tensor part B
      real*8 sig(ndir,ndir)      !Total Cauchy stress tensor
      real*8 s1,s2,s3,seq        !Principal stresses and failure stress
      real*8 sigc,Gf,sigc_min    !Fracture params
      real*8 sigc_m,sigc_std
      real*8 damage,eta         !Damage scalar
      real*8 epseq,epseq_u       !Equivalent strains
      real*8 epseqInc, epseq_c   !------"----------
      real*8 sigCo(ndir,ndir)    !Corotated stress
      real*8 sigTmp1(ndir,ndir)  !Temporary stress tensor
      real*8 sigTmp2(ndir,ndir)  !-"-
      real*8 tm1(3,3), tm2(3,3)  !Temporary 3x3 matrices
      real*8 tm3(3,3), tm4(3,3)  !-"-
      real*8 t1, t2, t3          !Temporary floats
      real*8 v1(3)               !Temporary 3x1 vector
      real*8 tol, maxit, res     !Convergence tolerance, max iterations and residual
      real*8 ID(3,3)             !identity matrix
      real*8 Nsub, chi           !Number of substeps and scaling factor
      real*8 Fq(3,3), L(3,3)     !Substep deformation gradient and vel. gradient
      real*8 Fqold(3,3)
      real*8 F(3,3), U(3,3)      !New deformation gradient and stretch tensor 
      real*8 Uinv(3,3), R(3,3)   !Inverse stretch tensor and rotation tensor
      real*8 Fold(3,3)           !Old deformation gradient
      real*8 deps(nshr+ndir),depsv 
      real*8 depsdev(nshr+ndir),tr
      real*8 stressOldtemp(ndir+nshr)
      real*8 De(3,3),Di(3,3),NiB(3,3)
      real*8 one,three,third
      parameter(one=1.,three=3.,third=one/three)

C     Common cluster between blocks
      integer numprocesses, kprocessnum    !Process numbers
      integer Ms,Ns, iseed, seeds(100)
      integer iflag,ix,iy,iz,ic,nx,ny,nz
      data iflag/0/
c
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,cluster(1)
      integer lc
      pointer(ptr_cluster,cluster)
      integer, parameter :: id_clust = 1
c
      call vgetnumcpus(numprocesses)
      call vgetrank(kprocessnum)

      ID(1,1)=1.
      ID(2,2)=1.
      ID(3,3)=1.
      ID(1,2)=0.
      ID(2,1)=0.
      ID(1,3)=0.
      ID(3,1)=0.
      ID(2,3)=0.
      ID(3,2)=0.

      conFlag=1

      intPt  = jInfoArray(2)

C======================================================================
C           Material params
C======================================================================
      muA    = props(1)	        !Initial shear modulus part A 
      lamLA  = props(2)	        !Locking stretch part A
      muB    = props(3)	        !Initial shear modulus part B
      kapB   = props(4)         !Bulk modulus part B
      alpha  = props(5)         !Pressure sensitivity parameter
      dgam0  = props(6)         !Reference shear strain rate
      m      = props(7)         !Strain rate sensitivity param
      tau0   = props(8)         !Initial shear strength
      tauss  = props(9)         !Steady state shear strength
      h      = props(10)        !Softening modulus
      chi    = props(11)        !Substepping factor
      maxit  = int(props(12))   !Maximum number of iterations
      tol    = props(13)        !Convergence tolerance
      GF     = props(14)        !Fracture energy
      dtc    = props(15)        !Critical stress duration
      sigc_m = props(16)        !Mean critical stress
      sigc_std = props(17)      !Std critical stress
      sigc_min = props(18)      !Minimum critical stress
      sigc_max = props(19)      !Maximum critical stress
      iseed  = props(20)        !Seed flag
      sflag  = props(21)        !Softening flag 1=linear, 2=exponential
      xmin   = props(22)        !1st x coordinate of assignment mesh
      xmax   = props(23)        !2nd x coordinate of assignment mesh
      ymin   = props(24)        !1st y coordinate of assignment mesh
      ymax   = props(25)        !2nd y coordinate of assignment mesh
      zmin   = props(26)        !1st z coordinate of assignment mesh
      zmax   = props(27)        !2nd z coordinate of assignment mesh
      nx     = props(28)        !Number of assignment els in x dir
      ny     = props(29)        !Number of assignment els in x dir
      nz     = props(30)        !Number of assignment els in x dir

C======================================================================
C           Assign stochastic failure strains
C====================================================================== 
      if (totaltime.eq.dt) then
            if (kprocessnum .eq. 0) then
                  if (iflag.eq.0) then
C                    Check if seed is specified by user
                  if (iseed.gt.0) call random_seed()
C                 Assign random failure stress to assignment mesh
                  lc = nx*ny*nz
                  write(*,*) 'Allocate cluster array with length =',lc
                  ptr_cluster = SMAFloatarrayCreate(id_clust,lc,0.0)
                  write(*,*) 'Finished allocating'
                  write(*,*) 'Fill cluster array'
                  do i = 1,lc
                     cluster(i)=gd(sigc_min,sigc_m,sigc_max,sigc_std)
                  end do
                  write(*,*) 'Finished filling cluster array'
                  write(*,*) cluster(1),cluster(lc)
                  iflag = 1
                  end if
            end if
      elseif(totaltime.gt.1.5*dt .and. totaltime.lt.2.5*dt) then
         ptr_cluster = SMAFloatarrayAccess(id_clust)
            do i=1,nblock
C                 Element numbers in assignment mesh
                  ix = nx*(coordmp(i,1)-xmin)/(xmax-xmin)
                  iy = ny*(coordmp(i,2)-ymin)/(ymax-ymin)
                  iz = nz*(coordmp(i,3)-zmin)/(zmax-zmin)
C                 Assignment element number corresponding to element k
                  ic = iz*nx*ny+iy*nx+ix+1
                  if (ic.lt.1 .or. ic.gt.nx*ny*nz) then
                     write(*,*) 'Element outside cluster'
                     write(*,*) ic,nx*ny*nz
                     write(*,*) coordmp(i,1),coordmp(i,2),coordmp(i,3)
                     stop
                  endif
c                 Assign failure stress to element from assignment mesh
                  stateNew(i,25) = cluster(ic)
            end do
      elseif (totaltime.gt.2.5*dt .and. totaltime.lt.3.5*dt) then
         if (kprocessnum .eq. 0) then
            if (iflag.eq.1) then
               write(*,*) 'Delete cluster array'
               call SMAFloatArrayDelete(id_clust)
               write(*,*) 'Finished deleting cluster array'
               iflag = 2
            endif
         endif
      end if
      do k = 1,nblock

C======================================================================
C           Dummy elastic step
C======================================================================
        if(stepTime .eq. 0.0 ) then
         t1 = mua+muB
         t2 = kapB
       
         deps(1) = strainInc(k,1)
         deps(2) = strainInc(k,2)
         deps(3) = strainInc(k,3)
         deps(4) = strainInc(k,4)
         if(nshr .eq. 3) then
         deps(5) = strainInc(k,5)
         deps(6) = strainInc(k,6)
         end if 
         depsv = deps(1)+deps(2)+deps(3)
         depsdev(1) = deps(1) - depsv*third
         depsdev(2) = deps(2) - depsv*third
         depsdev(3) = deps(3) - depsv*third
         depsdev(4) = deps(4)
         if(nshr .eq. 3) then
         depsdev(5) = deps(5)
         depsdev(6) = deps(6)
         end if 

C      Calculating Cauchy stress tensor of part A
         stressNew(k,1) = 2*t1*depsdev(1) + t2*depsv
         stressNew(k,2) = 2*t1*depsdev(2) + t2*depsv
         stressNew(k,3) = 2*t1*depsdev(3) + t2*depsv
         stressNew(k,4) = 2*t1*depsdev(4)
         if (nshr .eq. 3) then
         stressNew(k,5) = 2*t1*depsdev(5)
         stressNew(k,6) = 2*t1*depsdev(6)
         end if
        else

C======================================================================
C           Ordering tensors
C======================================================================      
C      Ordering the new def.grad. in a matrix       
         F(1,1) = defgradNew(k,1)
         F(2,2) = defgradNew(k,2)
         F(3,3) = defgradNew(k,3)
         F(1,2) = defgradNew(k,4)
         if(nshr .eq. 1) then
         F(2,1) = defgradNew(k,5)
         F(2,3) = 0.
         F(3,1) = 0.
         F(3,2) = 0.
         F(1,3) = 0.
         else
            F(2,3) = defgradNew(k,5)
            F(3,1) = defgradNew(k,6)
            F(2,1) = defgradNew(k,7)
            F(3,2) = defgradNew(k,8)
            F(1,3) = defgradNew(k,9)
         end if 

C      Ordering the old def.grad. in a matrix       
         Fold(1,1) = defgradOld(k,1)
         Fold(2,2) = defgradOld(k,2)
         Fold(3,3) = defgradOld(k,3)
         Fold(1,2) = defgradOld(k,4)
         if(nshr .eq. 1) then
         Fold(2,1) = defgradOld(k,5)
         Fold(2,3) = 0.
         Fold(3,1) = 0.
         Fold(3,2) = 0.
         Fold(1,3) = 0.
         else
            Fold(2,3) = defgradOld(k,5)
            Fold(3,1) = defgradOld(k,6)
            Fold(2,1) = defgradOld(k,7)
            Fold(3,2) = defgradOld(k,8)
            Fold(1,3) = defgradOld(k,9)
         end if

C      Ordering the new right stretch tensor in a matrix  
         U(1,1) = stretchNew(k,1)
         U(2,2) = stretchNew(k,2)
         U(3,3) = stretchNew(k,3)
         U(1,2) = stretchNew(k,4)
         if(nshr .eq. 1) then
         U(2,3) = 0.
         U(3,1) = 0.
         U(3,2) = 0.
         U(1,3) = 0.
         U(2,1) = U(1,2)
         else
         U(2,3) = stretchNew(k,5)
         U(3,1) = stretchNew(k,6)
         U(2,1) = stretchNew(k,4)
         U(3,2) = stretchNew(k,5)
         U(1,3) = stretchNew(k,6)
         end if

         stressOldtemp = stressOld(k,:)
C======================================================================
C           Calculating equivalent strain inc. and neccessary substeps
C======================================================================
      epstr = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
      epsdev(1,1) = strainInc(k,1)-epstr*third
      epsdev(1,2) = strainInc(k,4)
      epsdev(2,1) = epsdev(1,2)
c
      if (NTENS.gt.4) then
          epsdev(1,3) = strainInc(k,6)
          epsdev(2,3) = strainInc(k,5)
      else
          epsdev(1,3) = 0.
          epsdev(2,3) = 0.
      endif
c
      epsdev(2,2) = strainInc(k,2)-epstr*third
      epsdev(3,1) = epsdev(1,3)
      epsdev(3,2) = epsdev(2,3)
      epsdev(3,3) = strainInc(k,3)-epstr*third

      depseq = sqrt(2.*ddot(epsdev,epsdev)*third)

      epscr = chi*tau0/muB

      Nsub = nint(max(depseq/epscr,1.))

      call matinv(Fold,tm1)
      call matmult(F-Fold,tm1,L)
      L = L/dt
      call mattrans(L,tm1)
      rodInc = 0.5*(L+tm1)*dt
C======================================================================
C           Substep loop
C======================================================================
      STATEVtmp = stateOld(k,:)
      Fqold = Fold
      do q=1,Nsub !start substep loop
      tm2 = ID+dt*L/real(Nsub)
      call matmult(tm2,FqOld,Fq)


C     Calculating stresses from part A and B
      call rmap(props,Fq,STATEVtmp,dt/real(Nsub),sig)
      FqOld = Fq
      end do !end substep loop
      do i=1,23
            stateNew(k,i) = STATEVtmp(i)
      end do
      conFlag = STATEVtmp(31)
      stateNew(k,31) = conFlag
      stateNew(k,24) = Nsub

      damage = STATEVtmp(30)


C======================================================================
C       Rotation of stresses
C======================================================================

C     Inverting the right stretch tensor 
      call matinv(U,Uinv)

C     Calculate the rotation tensor R=FUinv
      call matmult(F,Uinv,R)

C     Rotate stresses to corotational configuration 
C     1st multiplication: sig*R
      call matmult(sig,R,tm1)

C     2nd multiplication: R^T*sig*R
      call mattrans(R,tm2)
      call matmult(tm2,tm1,sigCo)

C======================================================================
C           Principal stresses and check for failure
C======================================================================
       if(totalTime.ge.3*dt) then
            stateNew(k,25) = stateOld(k,25)
            sigc = stateNew(k,25)
       else
            sigc=1000000000.
       end if

      if(conFlag.eq.1) then
C     Calcualte principal stresses
            call princStress(sig,v1)
            s1 = v1(1)
            s2 = v1(2)
            s3 = v1(3)

            call principalNorm(s1,s2,s3,seq)
C     Calculate equivalent strain
            if(seq.eq.0.) then
                  epseqInc=0.
            else
                  epseqInc = ddot(sig,rodInc)/seq
            end if
            if(epseqInc.ge.0.) then
                  epseq = stateOld(k,26)+epseqInc
            else
                  epseq = stateOld(k,26)
            end if
            if(sflag.eq.1) then
                  if(seq.ge.sigc) then
                        duration = stateOld(k,29) +dt
                        stateNew(k,29) = duration
                        if(duration.ge.dtc) then
                              if(stateOld(k,27).eq.0.) then
                                    epseq_c = epseq
                                    epseq_u = epseq_c+2.*Gf/(seq*charLength(k))
                              else
                                    epseq_c = stateOld(k,27)
                                    epseq_u = stateOld(k,28)
                              end if
                              damage = 1.-(epseq_u-epseq)/(epseq_u-epseq_c)

                              stateNew(k,27) = epseq_c
                              stateNew(k,28) = epseq_u
                        end if
                  else
                        damage = STATEVtmp(30)
                        duration = 0.
                  end if
            else if(sflag.eq.2) then
                  if(seq.ge.sigc) then
                        duration = stateOld(k,29) +dt
                        stateNew(k,29) = duration
                        if(duration.ge.dtc) then
                              if(stateOld(k,27).eq.0.) then
                                    epseq_c = epseq
                              else
                                    epseq_c = stateOld(k,27)
                              end if
                              eta = sigc*charLength(k)/Gf
                              damage = 1.-exp(-eta*(epseq-epseq_c))
            
                              stateNew(k,27) = epseq_c
                              !stateNew(k,28) = epseq_u
                        end if
                  else
                        damage = STATEVtmp(30)
                        duration = 0.
                  end if
            else
                  write(*,*) "Not a valid softening flag"
            end if

            if (damage.ge.0.9) then
                  failFlag = 0
            else
                  failFlag = 1
            end if

      else
            failFlag = 0
      end if
      stateNew(k,31) = conFlag
      stateNew(k,32) = failFlag
      stateNew(k,26) = epseq
      stateNew(k,30) = damage

C======================================================================
C           Update stress
C======================================================================

         sigCo = (1.-damage)*sigCo

C       Update stressNew
         stressNew(k,1) = sigCo(1,1)
         stressNew(k,2) = sigCo(2,2)
         stressNew(k,3) = sigCo(3,3)
         stressNew(k,4) = sigCo(1,2)
         if(nshr .eq. 3) then
         stressNew(k,5) = sigCo(2,3)
         stressNew(k,6) = sigCo(3,1)
         end if
         end if
      end do
      return
      end
C#######################################################################
C                                                                      #
C                             FUNCTIONS                                #
C                                                                      #
C#######################################################################

C ----------------------------------------------------------------------
C sigeq: Equivalent principal stress
C Input: 3x1 principal stress vector
C Output: sigeq
C ----------------------------------------------------------------------
c
      subroutine principalNorm(s1,s2,s3,eqstress)
          implicit none
          real*8 eqstress,s1,s2,s3   
c
          s1 = (s1+abs(s1))*0.5
          s2 = (s2+abs(s2))*0.5
          s3 = (s3+abs(s3))*0.5
          eqstress = sqrt(s1**2+s2**2+s3**2)
      end
C ----------------------------------------------------------------------
C princStress: Calculate principal stresses of a stress tensor
C Input:   stress
C Output:  eigvals:   eigenvalues ordered
C ----------------------------------------------------------------------
c
      subroutine princStress(stress,eigvals)
          implicit none
          real*8 stress(3,3),eigvals(3), devstress(3,3),
     +           J2,J3,theta, det, hydstress, tr, pi, ddot, 
     +           one,two,three,third,half,frac
          parameter(one=1.,two=2.,three=3.,third=one/three,half=0.5)
          pi=4.d0*datan(1.d0)
          call dev_sub(stress,devstress)
          hydstress = tr(stress)*third
          J2 = ddot(devstress,devstress)*half
          J3 = det(devstress)
          if (J2.gt.1d-6) then
            theta = dacos(3.*sqrt(3.)*half*(J3/J2**(3./2.)))*third
          else
            theta = 0.
          end if
          frac = 2./sqrt(3.)*sqrt(J2)
          eigvals(1) = hydstress+frac*dcos(theta)
          eigvals(2) = hydstress+frac*dcos(2.*pi*third-theta)
          eigvals(3) = hydstress+frac*dcos(2.*pi*third+theta)
 
          return
      end

!-----------------------------------------------------------------------
!     Normal distribution function
!-----------------------------------------------------------------------
      real*8 function gd(dmin,d1m,dmax,dis)
      implicit none
      real*8 dmin,dmax,d1m,dis
      integer iset
      real*8 x,fac,gset,rsq,v1,v2,gd1,gd2
      save iset, gset
      data iset/0/
      if (iset.eq.0) then
    1    call random_number(x)
         v1 = 2.0*x-1.0
         call random_number(x)
         v2 = 2.0*x-1.0
         rsq=v1**2+v2**2
         if (rsq.ge.1..or.rsq.eq.0.) goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gd1=d1m+dis*v2*fac
         gd2=d1m+dis*v1*fac
         if(gd1.ge.dmin.and.gd2.ge.dmin) then
            gd=d1m+dis*v2*fac
         else 
            goto 1
         end if
         iset=1
      else
         gd=d1m+dis*gset
         iset=0
      endif
      return
      end
c
      SUBROUTINE eig(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )

*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2

*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)

*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

c      PRINT *, "DSYEVJ3: No convergence."
c      print *, "Matrix = ", A

      END SUBROUTINE
* End of subroutine DSYEVJ3
