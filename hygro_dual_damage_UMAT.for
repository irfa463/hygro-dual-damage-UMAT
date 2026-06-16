      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &   RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     &   TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &   NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     &   DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,
     &   NPT,LAYER,KSPT,KSTEP,KINC)

      IMPLICIT NONE

      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT
     &       ,LAYER,KSPT,KSTEP,KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
     &                ,DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SSE,SPD,SCD,RPL
     &                ,DDSDDT(NTENS,NTENS),DRPLDE(NTENS),DRPLDT
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS),TIME(2)
      DOUBLE PRECISION DTIME,TEMP,DTEMP,PREDEF(*),DPRED(*)
      DOUBLE PRECISION PROPS(NPROPS),COORDS(*)
      DOUBLE PRECISION DROT(3,3),PNEWDT,CELENT
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)

      INTEGER I,J
      DOUBLE PRECISION E1,E2,E3,NU12,NU13,NU23,G12,G13,G23
      DOUBLE PRECISION beta_h,moisrc,sec2day
      DOUBLE PRECISION kP(3),kH(3),mP,mH
      DOUBLE PRECISION bP(3),bH(3)
      DOUBLE PRECISION aXT_P,aXT_H,XT0,RRES
      DOUBLE PRECISION Csat,tau_diff_days,beta_trap,kappa,tau_P_days
      DOUBLE PRECISION P(3),H(3),Cm,Cb,C_total
      DOUBLE PRECISION ETOT(6),epsm(6),eps_h(6)
      DOUBLE PRECISION sig(6),Cij(6,6)
      DOUBLE PRECISION Eeff(3),Geff(3),nueff12,nueff13,nueff23
      DOUBLE PRECISION P_target,H_target,Pbar,Hbar
      DOUBLE PRECISION t_days,dtime_days
      DOUBLE PRECISION XT,sig1eq,expm,fT,softf
      DOUBLE PRECISION CKA(2,2),CKB(2),DET,Cm_new,Cb_new
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL=1.0D-30)

C=====================================================================
C  STEP 1: Read properties
C=====================================================================
      E1   = PROPS(1);  E2   = PROPS(2);  E3   = PROPS(3)
      NU12 = PROPS(4);  NU13 = PROPS(5);  NU23 = PROPS(6)
      G12  = PROPS(7);  G13  = PROPS(8);  G23  = PROPS(9)
      beta_h       = PROPS(10)
      moisrc       = PROPS(11)
      sec2day      = PROPS(12)
      kP(1)=PROPS(13); kP(2)=PROPS(14); kP(3)=PROPS(15)
      kH(1)=PROPS(16); kH(2)=PROPS(17); kH(3)=PROPS(18)
      mP   =PROPS(19); mH   =PROPS(20)
      bP(1)=PROPS(21); bP(2)=PROPS(22); bP(3)=PROPS(23)
      bH(1)=PROPS(24); bH(2)=PROPS(25); bH(3)=PROPS(26)
      aXT_P        = PROPS(27)
      aXT_H        = PROPS(28)
      XT0          = PROPS(29)
      RRES         = PROPS(30)
      Csat         = PROPS(31)
      tau_diff_days= PROPS(32)
      beta_trap    = PROPS(33)
      kappa        = PROPS(34)
      tau_P_days   = PROPS(35)

C=====================================================================
C  STEP 2: Time conversion
C=====================================================================
      t_days     = TIME(2) * sec2day
      dtime_days = DTIME   * sec2day

      DO I=1,6
        ETOT(I) = STRAN(I) + DSTRAN(I)
      END DO

C=====================================================================
C  STEP 3: Initialise STATEV on first increment
C=====================================================================
      IF (KSTEP.EQ.1 .AND. KINC.EQ.1) THEN
        DO I=1,NSTATV
          STATEV(I) = 0.D0
        END DO
      END IF

      P(1)=STATEV(1); P(2)=STATEV(2); P(3)=STATEV(3)
      H(1)=STATEV(4); H(2)=STATEV(5); H(3)=STATEV(6)
      Cm  =STATEV(7)
      Cb  =STATEV(8)

C=====================================================================
C  STEP 4: Carter-Kibler moisture — backward-Euler 2x2 ODE
C=====================================================================
      IF (NINT(moisrc) .EQ. 0) THEN
        C_total = MAX(0.D0, MIN(Csat, TEMP))
        Cm = C_total
        Cb = 0.D0
      ELSE
        IF (dtime_days .GT. 0.D0) THEN
          CKA(1,1) = 1.D0+dtime_days*(1.D0/tau_diff_days+beta_trap)
          CKA(1,2) = dtime_days*(1.D0/tau_diff_days - kappa)
          CKA(2,1) = -dtime_days*beta_trap
          CKA(2,2) = 1.D0 + dtime_days*kappa
          CKB(1)   = Cm + dtime_days*Csat/tau_diff_days
          CKB(2)   = Cb
          DET = CKA(1,1)*CKA(2,2) - CKA(1,2)*CKA(2,1)
          IF (ABS(DET) .LT. 1.D-20) THEN
            Cm_new = Cm+dtime_days*((Csat-Cm-Cb)/tau_diff_days
     &                              -beta_trap*Cm+kappa*Cb)
            Cb_new = Cb+dtime_days*(beta_trap*Cm-kappa*Cb)
          ELSE
            Cm_new=(CKB(1)*CKA(2,2)-CKA(1,2)*CKB(2))/DET
            Cb_new=(CKA(1,1)*CKB(2)-CKA(2,1)*CKB(1))/DET
          END IF
          Cm = MAX(0.D0, MIN(Csat,    Cm_new))
          Cb = MAX(0.D0, MIN(Csat-Cm, Cb_new))
        END IF
        C_total = MAX(0.D0, MIN(Csat, Cm+Cb))
      END IF

C=====================================================================
C  STEP 5: Reversible plasticisation P — driven by mobile moisture Cm
C           Exponential lag toward equilibrium P_target
C=====================================================================
      DO I=1,3
        IF (Cm .GT. 1.D-10) THEN
          P_target = 1.D0 - DEXP(-kP(I) * (Cm**mP))
        ELSE
          P_target = 0.D0
        END IF
        IF (dtime_days .GT. 0.D0 .AND. tau_P_days .GT. 0.D0) THEN
          expm = DEXP(-dtime_days / tau_P_days)
          P(I) = P_target + (P(I) - P_target)*expm
        ELSE
          P(I) = P_target
        END IF
        P(I) = MAX(0.D0, MIN(0.99D0, P(I)))
      END DO

C=====================================================================
C  STEP 6: Irreversible hydrolytic damage H — driven by bound moisture Cb
C           Ratchet: H is monotonically non-decreasing
C=====================================================================
      DO I=1,3
        IF (Cb .GT. 1.D-10 .AND. dtime_days .GT. 0.D0) THEN
          expm     = DEXP(-kH(I)*(Cb**mH)*dtime_days)
          H_target = 1.D0 - (1.D0 - H(I))*expm
          H(I)     = MAX(H(I), H_target)
        END IF
        H(I) = MAX(0.D0, MIN(0.99D0, H(I)))
      END DO

C=====================================================================
C  STEP 7: Hygroscopic swelling — isotropic, normal strains only
C=====================================================================
      DO I=1,6
        eps_h(I) = 0.D0
      END DO
      eps_h(1) = beta_h * C_total
      eps_h(2) = beta_h * C_total
      eps_h(3) = beta_h * C_total

C=====================================================================
C  STEP 8: Damage averages
C=====================================================================
      Pbar = (P(1)+P(2)+P(3)) / 3.D0
      Hbar = (H(1)+H(2)+H(3)) / 3.D0

C=====================================================================
C  STEP 9: Damage-degraded orthotropic stiffness (BOTH P and H on E)
C           Eeff(I) = E_i * (1 - bP(I)*P(I) - bH(I)*H(I))
C           This is the key correction vs CK-damage-2.for
C=====================================================================
      DO I=1,3
        Eeff(I) = PROPS(I)*(1.D0 - bP(I)*P(I) - bH(I)*H(I))
        Eeff(I) = MAX(0.05D0*PROPS(I), Eeff(I))
      END DO
      nueff12=NU12*(1.D0-0.5D0*(bP(1)*P(1)+bH(1)*H(1)
     &                         +bP(2)*P(2)+bH(2)*H(2)))
      nueff13=NU13*(1.D0-0.5D0*(bP(1)*P(1)+bH(1)*H(1)
     &                         +bP(3)*P(3)+bH(3)*H(3)))
      nueff23=NU23*(1.D0-0.5D0*(bP(2)*P(2)+bH(2)*H(2)
     &                         +bP(3)*P(3)+bH(3)*H(3)))
      nueff12 = MAX(0.01D0, MIN(0.49D0, nueff12))
      nueff13 = MAX(0.01D0, MIN(0.49D0, nueff13))
      nueff23 = MAX(0.01D0, MIN(0.49D0, nueff23))
      Geff(1)=G12*(1.D0-0.5D0*(bP(1)*P(1)+bP(2)*P(2)
     &                         +bH(1)*H(1)+bH(2)*H(2)))
      Geff(2)=G13*(1.D0-0.5D0*(bP(1)*P(1)+bP(3)*P(3)
     &                         +bH(1)*H(1)+bH(3)*H(3)))
      Geff(3)=G23*(1.D0-0.5D0*(bP(2)*P(2)+bP(3)*P(3)
     &                         +bH(2)*H(2)+bH(3)*H(3)))
      Geff(1) = MAX(0.05D0*G12, Geff(1))
      Geff(2) = MAX(0.05D0*G13, Geff(2))
      Geff(3) = MAX(0.05D0*G23, Geff(3))

      CALL ORTHO_STIFF(Eeff(1),Eeff(2),Eeff(3),
     &                 nueff12,nueff13,nueff23,
     &                 Geff(1),Geff(2),Geff(3), Cij)

C=====================================================================
C  STEP 10: Mechanical strain = total - hygroscopic
C=====================================================================
      DO I=1,6
        epsm(I) = ETOT(I) - eps_h(I)
      END DO

C=====================================================================
C  STEP 11: Elastic stress and consistent tangent
C=====================================================================
      DO I=1,6
        sig(I) = 0.D0
        DO J=1,6
          sig(I)      = sig(I) + Cij(I,J)*epsm(J)
          DDSDDE(I,J) = Cij(I,J)
        END DO
      END DO

      DO I=1,6
        STRESS(I) = sig(I)
      END DO

C=====================================================================
C  STEP 12: Tensile strength cap — fibre direction 1 only
C=====================================================================
      XT     = XT0*(1.D0 - aXT_P*P(1) - aXT_H*H(1))
      XT     = MAX(RRES*XT0, XT)
      sig1eq = STRESS(1)
      IF (sig1eq .GT. 0.D0) THEN
        fT = sig1eq / MAX(SMALL, XT)
        IF (fT .GT. 1.D0) THEN
          softf = MAX(RRES, 1.D0/fT)
          STRESS(1) = softf * STRESS(1)
          DO J=1,6
            DDSDDE(1,J) = 0.D0
            DDSDDE(J,1) = 0.D0
          END DO
        END IF
      END IF

C=====================================================================
C  STEP 13: Write STATEV back (only 8 variables — no viscous strains)
C=====================================================================
      STATEV(1)=P(1); STATEV(2)=P(2); STATEV(3)=P(3)
      STATEV(4)=H(1); STATEV(5)=H(2); STATEV(6)=H(3)
      STATEV(7)=Cm
      STATEV(8)=Cb

C=====================================================================
C  STEP 14: Strain energy density
C=====================================================================
      SSE = 0.D0
      DO I=1,NTENS
        SSE = SSE + 0.5D0*STRESS(I)*epsm(I)
      END DO
      SPD = 0.D0; SCD = 0.D0; RPL = 0.D0
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDT(I,J) = 0.D0
        END DO
        DRPLDE(I) = 0.D0
      END DO
      DRPLDT = 0.D0

C=====================================================================
C  STEP 15: Diagnostic output (element 1, GP 1 only)
C=====================================================================
      IF (NOEL.EQ.1 .AND. NPT.EQ.1) THEN
        IF (KINC.EQ.1 .OR. MOD(KINC,500).EQ.0) THEN
          WRITE(6,'(A,I2,A,I6,5(A,F10.5))')
     &      ' ST=',KSTEP,' INC=',KINC,
     &      ' t_d=',t_days,
     &      ' Cm=',Cm,' Cb=',Cb,
     &      ' P1=',P(1),' H1=',H(1)
        END IF
      END IF

      RETURN
      END

C=====================================================================
C  ORTHO_STIFF
C=====================================================================
      SUBROUTINE ORTHO_STIFF(E1,E2,E3,NU12,NU13,NU23,
     &                        G12,G13,G23, C)
      IMPLICIT NONE
      DOUBLE PRECISION E1,E2,E3,NU12,NU13,NU23,G12,G13,G23
      DOUBLE PRECISION C(6,6),S(6,6)
      DOUBLE PRECISION NU21,NU31,NU32
      INTEGER I,J
      NU21 = NU12*E2/E1
      NU31 = NU13*E3/E1
      NU32 = NU23*E3/E2
      DO I=1,6
        DO J=1,6
          S(I,J) = 0.D0
        END DO
      END DO
      S(1,1)= 1.D0/E1;  S(2,2)= 1.D0/E2;  S(3,3)= 1.D0/E3
      S(1,2)=-NU21/E2;  S(2,1)=-NU12/E1
      S(1,3)=-NU31/E3;  S(3,1)=-NU13/E1
      S(2,3)=-NU32/E3;  S(3,2)=-NU23/E2
      S(4,4)= 1.D0/G12; S(5,5)= 1.D0/G13; S(6,6)= 1.D0/G23
      CALL INV66(S,C)
      RETURN
      END

C=====================================================================
C  INV66: Gauss-Jordan with partial pivoting
C=====================================================================
      SUBROUTINE INV66(A,AINV)
      IMPLICIT NONE
      DOUBLE PRECISION A(6,6),AINV(6,6)
      DOUBLE PRECISION WORK(6,6),PIV,FACTOR
      INTEGER I,J,K,IMAX
      DO I=1,6
        DO J=1,6
          WORK(I,J) = A(I,J)
          AINV(I,J) = 0.D0
        END DO
        AINV(I,I) = 1.D0
      END DO
      DO I=1,6
        PIV  = ABS(WORK(I,I)); IMAX = I
        DO K=I+1,6
          IF (ABS(WORK(K,I)) .GT. PIV) THEN
            PIV = ABS(WORK(K,I)); IMAX = K
          END IF
        END DO
        IF (IMAX .NE. I) THEN
          DO J=1,6
            FACTOR=WORK(I,J);WORK(I,J)=WORK(IMAX,J)
            WORK(IMAX,J)=FACTOR
            FACTOR=AINV(I,J);AINV(I,J)=AINV(IMAX,J)
            AINV(IMAX,J)=FACTOR
          END DO
        END IF
        PIV = WORK(I,I)
        IF (ABS(PIV) .LT. 1.D-30) PIV = 1.D-30
        DO J=1,6
          WORK(I,J) = WORK(I,J)/PIV
          AINV(I,J) = AINV(I,J)/PIV
        END DO
        DO K=1,6
          IF (K .NE. I) THEN
            FACTOR = WORK(K,I)
            DO J=1,6
              WORK(K,J)=WORK(K,J)-FACTOR*WORK(I,J)
              AINV(K,J)=AINV(K,J)-FACTOR*AINV(I,J)
            END DO
          END IF
        END DO
      END DO
      RETURN
      END
