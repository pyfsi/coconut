      subroutine homogstressAtIntPoint_v5(
! Read only -
     1          Fnew, Fold,
     2          props, statOld,
! Write only -
     3          statNew, sNew)
      implicit none
!
!
!Declare input/out variables
!
      real*8, intent(in)        :: Fnew(3,3), Fold(3,3), props(11)
      real*8, intent(in)        :: statOld(247)
      real*8, intent(out)       :: statNew(247), sNew(6)
!
!
!     Local variables
      real*8, dimension(3,3) :: FFo, FFn, ROTORI0, ROTORI, iFFn, iFFo,
     1                          ROTLN, LL1, ROTORI00, rotref, ROT2, 
     2                          FFnref, FForef, iFForef
      real*8, dimension(6,6) :: invA66
      real*8, dimension(6,6) :: Q, C, AS66, QT
      real*8, dimension(7)   :: voce
      real*8, dimension(6)   :: Gammaferr, GammaCem, SGM, Bferr, BCem
      real*8, dimension(6)   :: Kferr, KCem
      real*8, dimension(6)   :: RR, defarrtot
      real*8, dimension(6)   :: sferr0, sCem0, sfinal
      real*8, dimension(6)   :: strhatferr, strhatCem, sCferr, sCCem
      real*8, dimension(6)   :: sChatferr, sChatCem
      real*8, dimension(4)   :: axang0
      real*8, dimension(2)   :: JJ, ff, ldot
      real*8, dimension(3)   :: vect
      real*8, dimension(6)   :: phihat
      real*8, dimension(6)   :: shatferrtro, shatCemtro, defhattot
      real*8, dimension(6)   :: shatferr0, shatCem0
      real*8, dimension(6)   :: erateplhatCem, erateplCem
      real*8, dimension(6)   :: erateplhatferr, erateplferr
      real*8 enoutferr, enoutCem
      real*8  YCem, ferrat, frac, sintlam, lameff, peeq
      real*8  sY, sY33, sYc, num, Cperc, flocferr, flocCem, colsize
      real*8  N3, D3, N6, D6, N5, D5, phi3, phi5, phi6
      real*8  omegaferr, omegaCem
      real*8  dlferr, dlCem, lambda, mu, E, nu, hard
      real*8  seqferr, seqCem, nldot
      real*8 axint(3,21), wtint(21), wtnew(21), detF
      integer*4 j1, j2
      !
      ! Parameters
      !
      real*8, parameter :: zero = 0.d0, one = 1.d0, two = 2.d0,
     1      fracini=8.8d-1, sqrt3=1.732050808d0, dlmin=1.d-08,
     2      sqrt2=1.414213562373d0,  half=5.d-1,  rad2deg=5.7295780d1,
     3      C1=5.31315d1, C2=1.6039d-6, C3=-3.925d-1,
     4      sqrt26=8.1649658092773d-1, three=3.d0, four=4.d0,
     5      isqrt3=5.7735026919d-1, isqrt2=7.07106781187d-1,
     6      isqrt6=4.08248290464d-1, s2=7.071067811865476d-1,
     7      s3=3.87907304067d-1, cos3=8.360955967489268d-1,
     8      w1=2.65214244093d-2, w2=1.99301476312d-2,
     9      w3=2.50712367487d-2, eps=1.d-16
      integer, parameter :: nvar=11
      !
      ! Code
      !
      E       = props(1)        !Young's modulus
      nu      = props(2)        !Poisson's ratio
      voce(1) = props(3)        !Tau Peierls ferrite
      voce(2) = props(4)        !Pre-factor tau(s)
      voce(3) = props(5)        !Additive factor = 0.68
      voce(4) = props(6)        !Voce saturation stress
      voce(5) = props(7)        !Voce initial disloc density
      voce(6) = props(8)        !Voce exponential factor
      YCem    = props(9)        !Yield stress cementite
      voce(7) = props(10)       !Burgers modulus
      Cperc   = props(11)       !C content (wt%)
!
! Initiate state variables
      do j1=1, 247
        statNew(j1)=statOld(j1)
      enddo
!
!
      axang0(1) = statOld(1)
      axang0(2) = statOld(2)
      axang0(3) = statOld(3)
      axang0(4) = statOld(4)
      call axang2rot(axang0,ROTORI00)   !Initial colony orientation (Colony->Ref)
!
!   Update Stored deformation gradient
!
      FForef(1,1)=statOld(5)
      FForef(1,2)=statOld(6)
      FForef(1,3)=statOld(7)
      FForef(2,1)=statOld(8)
      FForef(2,2)=statOld(9)
      FForef(2,3)=statOld(10)
      FForef(3,1)=statOld(11)
      FForef(3,2)=statOld(12)
      FForef(3,3)=statOld(13)
      call inverse3(Fold, iFForef)
      FFnref=matmul(Fnew, matmul(iFForef, FForef))
      statNew(5)=FFnref(1,1)
      statNew(6)=FFnref(1,2)
      statNew(7)=FFnref(1,3)
      statNew(8)=FFnref(2,1)
      statNew(9)=FFnref(2,2)
      statNew(10)=FFnref(2,3)
      statNew(11)=FFnref(3,1)
      statNew(12)=FFnref(3,2)
      statNew(13)=FFnref(3,3)
!
! Bazant's weights and axes (surface integration unit sphere)
!
      axint(1:3,1)=(/one, zero, zero/)
      axint(1:3,2)=(/zero, one, zero/)
      axint(1:3,3)=(/zero, zero, one/)
      axint(1:3,4)=(/s2, s2, zero/)
      axint(1:3,5)=(/s2, -s2, zero/)
      axint(1:3,6)=(/s2, zero, s2/)
      axint(1:3,7)=(/s2, zero, -s2/)
      axint(1:3,8)=(/zero, s2, s2/)
      axint(1:3,9)=(/zero, s2, -s2/)
      axint(1:3,10)=(/s3, s3, cos3/)
      axint(1:3,11)=(/s3, s3, -cos3/)
      axint(1:3,12)=(/s3, -s3, cos3/)
      axint(1:3,13)=(/s3, -s3, -cos3/)
      axint(1:3,14)=(/s3, cos3, s3/)
      axint(1:3,15)=(/s3, cos3, -s3/)
      axint(1:3,16)=(/s3, -cos3, s3/)
      axint(1:3,17)=(/s3, -cos3, -s3/)
      axint(1:3,18)=(/cos3, s3, s3/)
      axint(1:3,19)=(/cos3, s3, -s3/)
      axint(1:3,20)=(/cos3, -s3, s3/)
      axint(1:3,21)=(/cos3, -s3, -s3/)
      !weights
      wtint(1:21)=(/w1,w1,w1,w2,w2,w2,w2,w2,w2,
     1       w3,w3,w3,w3,w3,w3,w3,w3,w3,w3,w3,w3/)
      wtint=wtint*two
!
! Calculate the thickness ratio of ferrite due to composition
!
      call thrat(Cperc,ferrat)     !ferrat is the ferrite volume fraction relative to eutectoid
      frac=fracini*ferrat
! 
!   Elasticity ! tensors
! 
      lambda=E*nu/(one-two*nu)/(one+nu)
      mu=E/two/(one+nu)
      !Initialize diagonal elastic stiffness matrix
      SGM=zero
      SGM(1)=three*lambda+two*mu
      SGM(2)=two*mu
      SGM(3)=two*mu
      SGM(4)=mu
      SGM(5)=mu
      SGM(6)=mu
      !elastic stiffness matrix
      C=zero
      C(1:3,1:3)=lambda
      do j1=1,6
        C(j1,j1)=C(j1,j1)+mu
      enddo
      do j1=1,3
        C(j1,j1)=C(j1,j1)+mu
      enddo
      !
      !
      ! Initialize Q and QT
      !
      !
      Q(1,:)=[isqrt3, -isqrt2, -isqrt6, zero, zero, zero]
      Q(2,:)=[isqrt3, isqrt2, -isqrt6, zero, zero, zero]
      Q(3,:)=[isqrt3, zero, two*isqrt6, zero, zero, zero]
      Q(4,:)=[zero, zero, zero, one, zero, zero]
      Q(5,:)=[zero, zero, zero, zero, one, zero]
      Q(6,:)=[zero, zero, zero, zero, zero, one]
      QT=transpose(Q)
      !
      !Initialize stress and tangent matrix
      !
      sNew=zero
      enoutferr = zero
      enoutCem  = zero
      !
      ! Iteration over all the orientations
      !
      do j2 = 1,21
          axang0(1)=axint(1,j2)
          axang0(2)=axint(2,j2)
          axang0(3)=axint(3,j2)
          axang0(4)= zero
          call axang2rot(axang0,rotref)   !Initial reference colony orientation
          ROTORI0=matmul(ROTORI00,rotref)
          !
          ! Deformation gradient at the beginning of increment at colony frame
          !
          FFo=matmul(transpose(ROTORI0), matmul(FForef, ROTORI0))
          !
          ! Deformation gradient at the end of increment at colony frame
          !
          FFn=matmul(transpose(ROTORI0), matmul(FFnref, ROTORI0))
            
          ! Inverse of F-s
          call inverse3(FFo, iFFo)
          call inverse3(FFn, iFFn)
          !
          ! New weights
          !
          call determ(FFn, detF)
          vect=detF*matmul(transpose(iFFn),ROTORI0(:,3))
          !
          call norm3(vect, wtnew(j2))
          wtnew(j2)=wtint(j2)*wtnew(j2)
          !
          ! Corresponding colony rotation at the end of the increment
          !
          call calcrot(FFn,iFFn, ROTLN)
          !
          !
          !
          ! New colony orientation
          !
          ROTORI=matmul(ROTORI0,ROTLN)
          !
          !!Rotate FROM reference system to colony system (at the end of the increment)
          !
          ROT2=transpose(ROTORI)
          !
          ! ! get 6x6 rotation matrices
          call rot2A(ROT2, AS66,1)
          call rot2A(transpose(ROT2), invA66,0)
          !
          !! Rotate LL1 and get strains at colony system
          !
          LL1=(matmul(FFn,iFFo)-matmul(FFo,iFFn))*half  !DSTRAIN-ekin konfirmatua
          LL1=matmul(transpose(ROTLN), matmul(LL1, ROTLN))
          defarrtot(1)=LL1(1,1)
          defarrtot(2)=LL1(2,2)
          defarrtot(3)=LL1(3,3)
          defarrtot(4)=(LL1(1,2)+LL1(2,1))
          defarrtot(5)=(LL1(1,3)+LL1(3,1))
          defarrtot(6)=(LL1(3,2)+LL1(2,3))
          
    !
    !Stress initialization
    !
          sferr0(1)=statOld((j2-1)*nvar+17)
          sferr0(2)=statOld((j2-1)*nvar+18)
          sferr0(3)=statOld((j2-1)*nvar+19)
          sferr0(4)=statOld((j2-1)*nvar+20)
          sferr0(5)=statOld((j2-1)*nvar+21)
          sferr0(6)=statOld((j2-1)*nvar+22)
          sCem0(1)=statOld((j2-1)*nvar+23)
          sCem0(2)=statOld((j2-1)*nvar+24)
          sCem0(3)=statOld((j2-1)*nvar+19)
          sCem0(4)=statOld((j2-1)*nvar+25)
          sCem0(5)=statOld((j2-1)*nvar+21)
          sCem0(6)=statOld((j2-1)*nvar+22)
          
    !            
    !            
    !            
          call norm3(iFFn(3,1:3),sintlam)
          !call norm3(iFFo(3,1:3), sintlamo)
          lameff=statOld(14)
          lameff=lameff/sintlam
          peeq=statOld(15+(j2-1)*nvar)
          colsize=C1*(peeq+C2)**C3
          call gettau2(voce,lameff*ferrat,peeq, sY, hard)           !gettau-rekin konparatu dut
          call gettau2(voce,lameff*sqrt2*ferrat,peeq, sY33, hard)
          call gettau2(voce,colsize,peeq, sYc, hard)

          !STATEV(1+(j2-1)*nvar+4)=colsize
          !STATEV(2+(j2-1)*nvar+4)=lameff      
          !STATEV(7+(j2-1)*nvar+4)=sY
    !
          RR(1)=one
          RR(2)=RR(1)
          RR(3)=sY33/sY
          RR(4)=RR(1)
          RR(5)=sYc/sY
          RR(6)=RR(5)
          !STATEV(9+(j2-1)*nvar+4)=RR(3)*sY
          !STATEV(10+(j2-1)*nvar+4)=RR(5)*sY
          ! Calculate Gamma for both phases
          call getgamma(RR, Gammaferr)   !Ferrite
          RR=one !Isotropic cementite
          call getgamma(RR, GammaCem)    !Cementite
    !
    !
    !

          !
          !Initialize variables
          !
          nldot=one
          num=zero
          !
          ! Transformed initial stress and strains
          !
          defhattot=matmul(QT,defarrtot)
          shatferr0=matmul(QT,sferr0)
          shatCem0=matmul(QT,sCem0)
          !
          !Initialize shatferrtro and shatCemtro
          !
          shatferrtro=shatferr0+SGM*defhattot
          shatCemtro=shatCem0+SGM*defhattot
          !
          ! Initialize B and K
          !
          Bferr=one
          BCem=one
          Kferr=Gammaferr
          KCem=GammaCem
          !
          ! Initialize phi
          !
          N3=sqrt3*((shatCemtro(1)/BCem(1)-shatferrtro(1)/Bferr(1))+
     1           sqrt2*(shatCemtro(3)/BCem(3)-shatferrtro(3)/Bferr(3)))
          D3=SGM(1)*((one-frac)/Bferr(1)+frac/BCem(1))+
     1           two*SGM(3)*((one-frac)/Bferr(3)+frac/BCem(3))
          N5=(shatCemtro(5)-shatferrtro(5))
          D5=SGM(5)
          N6=(shatCemtro(6)-shatferrtro(6))
          D6=SGM(6)
          phi3=N3/D3
          phi5=N5/D5
          phi6=N6/D6
          phihat=[phi3/sqrt3,zero,phi3*sqrt26,zero,phi5,phi6]
          !
          ! Initialize strials
          !
          strhatferr=shatferrtro+SGM*phihat*(one-frac)
          strhatCem=shatCemtro-SGM*phihat*frac
          !
          ! Initialize Newton array
          !
          seqferr=sqrt(sum(half*strhatferr*Kferr*strhatferr))
          dlferr=zero
          seqCem=sqrt(sum(half*strhatCem*KCem*strhatCem))
          dlCem=zero
          !
          call gettau2(voce,lameff*ferrat,peeq+dlferr, sY, hard) ! YCem is constant
          !
          flocferr=seqferr-sY
          flocCem=seqCem-YCem
          ! Calculate sY and hardening rate (ferrite)
          !
          !
          if ((flocferr>zero).or.(flocCem>zero)) then
          !
          ! Start Newton-Raphson
          !
           do while (nldot>dlmin)
            !
            ! Construct f1-f2
            !
            seqferr=sqrt(half*sum(strhatferr*Kferr*strhatferr))
            seqCem=sqrt(half*sum(strhatCem*KCem*strhatCem))
            ff(1)=((seqferr-sY)+abs(seqferr-sY))/two
            ff(2)=((seqCem-YCem)+abs(seqCem-YCem))/two
            !!
            !!=== Construct Jacobian ===
            !!
            !
            omegaferr=sum(Kferr*Kferr*Bferr*SGM*strhatferr*strhatferr)
            omegaCem=sum(KCem*KCem*BCem*SGM*strhatCem*strhatCem)
            !
            !!J ferr
            JJ(1)=-(one-hard/sY*dlferr)/(four*sY*(ff(1)+sY))*omegaferr
            JJ(2)=-(one)/(four*YCem*(ff(2)+YCem))*omegaCem
            !!
            ! Calculate new dl values from Newton iteration
            ldot=-ff/JJ
            nldot=sqrt(sum(ldot*ldot))
            dlferr=dlferr+ldot(1)
            dlCem=dlCem+ldot(2)
            !
            !
            !Update variables for new iteration
            !
            !
            call gettau2(voce,lameff*ferrat,peeq+dlferr, sY, hard) ! YCem is constant
            Bferr=(one+dlferr/(two*sY)*SGM*Gammaferr)
            Kferr=Gammaferr/Bferr/Bferr
            BCem=(one+dlCem/(two*YCem)*SGM*GammaCem)
            KCem=GammaCem/BCem/BCem
            ! Update phi
            N3=sqrt3*((shatCemtro(1)/BCem(1)-shatferrtro(1)/Bferr(1))+
     1           sqrt2*(shatCemtro(3)/BCem(3)-shatferrtro(3)/Bferr(3)))
            D3=SGM(1)*((one-frac)/Bferr(1)+frac/BCem(1))+
     1           two*SGM(3)*((one-frac)/Bferr(3)+frac/BCem(3))
            N5=(Bferr(5)*shatCemtro(5)-BCem(5)*shatferrtro(5))
            D5=SGM(5)*(BCem(5)*(one-frac)+Bferr(5)*frac)
            N6=(Bferr(6)*shatCemtro(6)-BCem(6)*shatferrtro(6))
            D6=SGM(6)*(BCem(6)*(one-frac)+Bferr(6)*frac)
            phi3=N3/D3
            phi5=N5/D5
            phi6=N6/D6
            phihat=[phi3/sqrt3,zero,phi3*sqrt26,zero,phi5,phi6]
            !
            !Update trial stresses
            !
            strhatferr=shatferrtro+SGM*phihat*(one-frac)
            strhatCem=shatCemtro-SGM*phihat*frac
            !
            !Verify
            !
            !
            !
            !
            ! Initialize strials
            num=num+one
           enddo
          endif
          !
          ! Calculate ferr/Cem Cauchy stresses
          !
          sChatferr=strhatferr/Bferr
          erateplhatferr=dlferr/two/(seqferr+eps)*Gammaferr*sChatferr
          sCferr=matmul(Q,sChatferr)
          erateplferr=matmul(Q,erateplhatferr)
          sChatCem=strhatCem/BCem
          sCCem=matmul(Q,sChatCem)
          erateplhatCem=dlCem/two/(seqCem+eps)*GammaCem*sChatCem
          erateplCem=matmul(Q,erateplhatCem)
          !
          ! Calculate the out-of-plane plastic energy dissipated 
          !
          erateplferr(3)=(erateplferr(3)+abs(erateplferr(3)))
     1       /two
          erateplCem(3)=(erateplCem(3)+abs(erateplCem(3)))
     1       /two
          enoutferr=enoutferr+wtnew(j2)*
     1       (sCferr(3)+abs(sCferr(3)))/two*erateplferr(3)
          enoutCem=enoutCem+wtnew(j2)*
     1       (sCCem(3)+abs(sCCem(3)))/two*erateplCem(3)
          !
          !
          ! Calculate element Cauchy stress
          !
          sfinal=sCferr*frac+sCCem*(one-frac)
          ! print *,'defhattot: ', j2,', ', defhattot
!
          sNew=sNew+wtnew(j2)*matmul(invA66,sfinal)
    !
    ! Update SDV
    !
          statNew(15+(j2-1)*nvar)=dlferr+statOld(15+(j2-1)*nvar)
          statNew(16+(j2-1)*nvar)=dlCem+statOld(16+(j2-1)*nvar)
          do j1=1,6
                statNew(16+j1+(j2-1)*nvar)=sCferr(j1)
          enddo
          statNew(23+(j2-1)*nvar)=sCCem(1)
          statNew(24+(j2-1)*nvar)=sCCem(2)
          statNew(25+(j2-1)*nvar)=sCCem(4)
      enddo
      ! print *,'wtnew2: ', (wtnew)
      ! print *,'STATEV: ', STATEV(1:25)
      sNew=sNew/sum(wtnew)
      wtnew=wtnew/sum(wtnew)
      statNew(246)=statOld(246)+enoutferr/sum(wtnew)
      statNew(247)=statOld(247)+enoutCem/sum(wtnew)
      !do j2=1,21
      !    STATEV(22+(j2-1)*nvar+4)=wtnew(j2)
      !enddo
      RETURN
      END subroutine homogstressAtIntPoint_v5
!
!
!     
!=============     
! FUNCTIONS  |
!=============
!
      !
      !
      subroutine determ(s, J)
      implicit none
      real*8, intent(in)        :: s(3,3)
      real*8, intent(out)       :: J
      J=s(1,1)*s(2,2)*s(3,3)+
     1  s(2,1)*s(3,2)*s(1,3)+
     2  s(1,2)*s(2,3)*s(3,1)-
     3  s(1,3)*s(2,2)*s(3,1)-
     4  s(1,2)*s(2,1)*s(3,3)-
     5  s(1,1)*s(2,3)*s(3,2)
      return
      end subroutine determ
      !
      !
      subroutine thrat(Cperc2, ferrat)
      implicit none
      real*8, intent(in)  :: Cperc2
      real*8, intent(out) :: ferrat
      real*8 Vfer,Cperc
      real*8, parameter :: mat12=-1.3988006d1,
     1      Vref=8.845923538d-1, hund=1.d2, one=1.d0
      !
      Cperc=Cperc2/hund
      Vfer=(one-Cperc)+mat12*Cperc
      ferrat=Vfer/Vref
      return
      end subroutine thrat
      !
      !
      subroutine gettau2(voce,lameff2,peeq, tau, hard)
      implicit none
      real*8, intent(in)  :: voce(7), lameff2, peeq
      real*8, intent(out)  :: tau, hard
      real*8, parameter :: one=1.d0, M=1.d0, sqrt3=1.732050808d0
      real*8 lameff
            if (lameff2.le.voce(7)) then
                  lameff=voce(7)
            else
                  lameff=lameff2
            endif
            tau=sqrt3*(voce(1)+voce(2)/lameff*(log(lameff/voce(7))+
     1            voce(3))+voce(5)+(voce(4)-voce(5))*
     2            (one-exp(-voce(6)*M*peeq)))
            hard=sqrt3*voce(6)*M*(voce(4)-voce(5))*exp(-voce(6)*M*peeq)
      return
      end subroutine gettau2
!
!
      subroutine axang2rot(axang, rot)
      implicit none
      real*8, intent(in)    :: axang(4)
      real*8, intent(out)   :: rot(3,3)
      real*8, dimension(3)  :: R3, uvect, vvect
      real*8 nr, n1, n2, n3, mxval, cth, sth
      integer i1, ord(3), imin
      real*8, parameter :: m1=2.d0
      ! Order of normal to the lamella
      R3(1:3)=abs(axang(1:3))
      mxval=m1
      do i1=1,3
            if (R3(i1).lt.mxval) then
                  mxval=R3(i1)
                  imin=i1
            endif
      enddo
      ord(3)=imin
      ord(1)=imin+1
      ord(2)=imin+2
      do i1=1,2
            if (ord(i1).gt.3) then
                  ord(i1)=ord(i1)-3
            endif
      enddo
      n1=axang(ord(1))
      n2=axang(ord(2))
      n3=axang(ord(3))
      nr=sqrt(1.d0-n3*n3)
      
      uvect(ord(1))=-n2/nr
      uvect(ord(2))=n1/nr
      uvect(ord(3))=0.d0
      vvect(ord(1)) = -n1*n3/nr
      vvect(ord(2)) = -n2*n3/nr
      vvect(ord(3)) = nr
      cth=cos(axang(4))
      sth=sin(axang(4))
      rot(ord(1),1)=uvect(ord(1))*cth+vvect(ord(1))*sth
      rot(ord(1),2)=vvect(ord(1))*cth-uvect(ord(1))*sth
      rot(ord(2),1)=uvect(ord(2))*cth+vvect(ord(2))*sth
      rot(ord(2),2)=vvect(ord(2))*cth-uvect(ord(2))*sth
      rot(ord(3),1)=nr*sth
      rot(ord(3),2)=nr*cth
      rot(1:3,3)=axang(1:3)
      return
      end subroutine axang2rot
      !
      !
      !
      SUBROUTINE rot2axang(R, axang)
      implicit none
      real*8, intent(in)    :: R(3,3)
      real*8, intent(out)   :: axang(4)
      real*8, dimension(3)  :: R3
      real*8 mxval
      integer i1, imin
      real*8, parameter :: m1=2.d0
      axang(1:3)=R(1:3,3)
      ! Order of normal to the lamella
      R3(1:3)=abs(R(1:3,3))
      imin=0
      mxval=m1
      do i1=1,3
            if (R3(i1).lt.mxval) then
                  mxval=R3(i1)
                  imin=i1
            endif
      enddo
      axang(4)=atan2(R(imin,1),R(imin,2))
      return
      end subroutine rot2axang
!
!
      subroutine calcrot(FF,iFF,ROT)
      implicit none
      real*8, intent(in)     :: FF(3,3),iFF(3,3)
      real*8, intent(out)    :: ROT(3,3)
      real*8, dimension(3)   :: u, v, w2
      real*8, dimension(3,3) :: G
      real*8 upar, nu, nv, nw2
      integer*4 i1
      G=transpose(iFF)  !F^(-T)
      u(1)=FF(1,1)
      u(2)=FF(2,1)
      u(3)=FF(3,1)
      call norm3(u,nu)
      u=u/nu
      upar=0.d0
      do i1=1,3
            upar=upar+u(i1)*FF(i1,2)
      enddo
      do i1=1,3
            v(i1)=FF(i1,2)-upar*u(i1)
      enddo
      call norm3(v,nv)
      v=v/nv
      w2=G(:,3)
      call norm3(w2,nw2)
      w2=w2/nw2
      do i1=1,3
            ROT(i1,:)=[u(i1),v(i1),w2(i1)]
      enddo
      return
      end subroutine calcrot
!
!
      subroutine norm3(vect, norm)
            implicit none
            real*8, intent(in)  :: vect(3)
            real*8, intent(out) :: norm
            norm=sqrt(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))
            return
      end subroutine norm3
!
!
      !
      subroutine getgamma(R,Gamma)
      implicit none
      real*8, intent(in)     :: R(6)
      real*8, intent(out)    :: Gamma(6)
      real*8 N, M
      real*8, parameter :: zero=0.d0, six=6.d0, one=1.d0
      real*8, parameter :: four=4.d0, three=3.d0
      !
      N=one/R(3)/R(3)
      M=one/R(5)/R(5)
      !
      !
      Gamma=[zero,four-N, three*N, six, six*M, six*M]
      return
      end subroutine getgamma
      !
      subroutine inverse3(M, invM)
      implicit none
      real*8, intent(in)    ::  M(3,3)
      real*8, intent(out)   ::  invM(3,3)
      !
      real*8 detM
      !
      detM=M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-M(2,1)*M(1,2)
     1      *M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*M(1,2)*M(2,3)-M(3,1)
     2      *M(1,3)*M(2,2)
      !
      invM(1,1)=(M(2,2)*M(3,3)-M(2,3)*M(3,2))/detM
      invM(1,2)=-(M(1,2)*M(3,3)-M(1,3)*M(3,2))/detM
      invM(1,3)=(M(1,2)*M(2,3)-M(1,3)*M(2,2))/detM
      !__
      invM(2,1)=-(M(2,1)*M(3,3)-M(2,3)*M(3,1))/detM
      invM(2,2)=(M(1,1)*M(3,3)-M(1,3)*M(3,1))/detM
      invM(2,3)=-(M(1,1)*M(2,3)-M(1,3)*M(2,1))/detM
      !__
      invM(3,1)=(M(2,1)*M(3,2)-M(2,2)*M(3,1))/detM
      invM(3,2)=-(M(1,1)*M(3,2)-M(1,2)*M(3,1))/detM
      invM(3,3)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))/detM
      !
      return
      end subroutine inverse3
      !
      !
      !
      subroutine rot2A(R,A,labl)
      implicit none
      real*8, intent(in)     :: R(3,3)
      integer*4, intent(in)  :: labl
      real*8, intent(out)    :: A(6,6)
      integer*4 i1, i2
      real*8 fact
      real*8, parameter      :: two=2.d0, half=5.d-1, one=1.d0
      fact=one
      if (labl==1) then
        fact=half
      endif
      do i1=1,3
        do i2=1,3
            A(i1,i2)=R(i1,i2)*R(i1,i2)
        enddo
      enddo
      A(1,4)=two*R(1,1)*R(1,2)*fact
      A(1,5)=two*R(1,1)*R(1,3)*fact
      A(1,6)=two*R(1,3)*R(1,2)*fact
      A(2,4)=two*R(2,1)*R(2,2)*fact
      A(2,5)=two*R(2,1)*R(2,3)*fact
      A(2,6)=two*R(2,2)*R(2,3)*fact
      A(3,4)=two*R(3,1)*R(3,2)*fact
      A(3,5)=two*R(3,1)*R(3,3)*fact
      A(3,6)=two*R(3,2)*R(3,3)*fact
      !
      A(4,1)=R(2,1)*R(1,1)/fact
      A(4,2)=R(2,2)*R(1,2)/fact
      A(4,3)=R(2,3)*R(1,3)/fact
      A(5,1)=R(1,1)*R(3,1)/fact
      A(5,2)=R(1,2)*R(3,2)/fact
      A(5,3)=R(1,3)*R(3,3)/fact
      A(6,1)=R(2,1)*R(3,1)/fact
      A(6,2)=R(2,2)*R(3,2)/fact
      A(6,3)=R(2,3)*R(3,3)/fact
      !
      !
      A(4,4)=R(1,1)*R(2,2)+R(1,2)*R(2,1)
      A(4,5)=R(1,1)*R(2,3)+R(1,3)*R(2,1)
      A(4,6)=R(1,2)*R(2,3)+R(1,3)*R(2,2)
      A(5,4)=R(1,1)*R(3,2)+R(1,2)*R(3,1)
      A(5,5)=R(1,1)*R(3,3)+R(1,3)*R(3,1)
      A(5,6)=R(1,2)*R(3,3)+R(1,3)*R(3,2)
      A(6,4)=R(2,1)*R(3,2)+R(2,2)*R(3,1)
      A(6,5)=R(2,1)*R(3,3)+R(2,3)*R(3,1)
      A(6,6)=R(2,2)*R(3,3)+R(2,3)*R(3,2)
      return
      end subroutine rot2A
      !
