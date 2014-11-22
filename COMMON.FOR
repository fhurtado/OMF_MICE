C ************************************************************************
C
C     This library contains all the routines related to generating random
C     variables and condutcing numerical work used by assess.for
C         XNORM   - Generates a normal random variate
C         STUDENT - Generates a t-random variables
C         RAN?    - Generates a uniform random variate
C         GENMULT - Generates a multivariate random variate (normal / t)
C         CUMNOM  - Cumulative normal distribution
C         ERFCC   - Used by CUMNOM
C         CALC_AG - Used to calculate age-reading error given %agreement between readers
C
C ************************************************************************
C
      REAL*8 FUNCTION XNORM(INO,MEAN,SIGG,ISEED)
C
C     Uniform to Normal conversion routine
C
      INTEGER ISEED,INO
      REAL*8 Z1, Z2, MEAN, SIGG
      REAL*8 RAN1,RAN2,RAN3,RAN4
C
C     Ignore zero calls
      IF (SIGG.EQ.0.0) THEN
        XNORM = 0.0
        RETURN
      ENDIF
C
      IF (INO.EQ.1) THEN
        Z1 = RAN1(ISEED)
   1    Z2 = RAN1(ISEED)
        IF (Z2.LE.0.OR.Z2.GT.1) GOTO 1
      ENDIF
      IF (INO.EQ.2) THEN
        Z1 = RAN2(ISEED)
   2    Z2 = RAN2(ISEED)
        IF (Z2.LE.0.OR.Z2.GT.1) GOTO 2
      ENDIF
      IF (INO.EQ.3) THEN
        Z1 = RAN3(ISEED)
   3    Z2 = RAN3(ISEED)
        IF (Z2.LE.0.OR.Z2.GT.1) GOTO 3
      ENDIF
      IF (INO.EQ.4) THEN
        Z1 = RAN4(ISEED)
   4    Z2 = RAN4(ISEED)
        IF (Z2.LE.0.OR.Z2.GT.1) GOTO 4
      ENDIF
C
      XNORM = SIN(6.238319*Z1) * SQRT(-2.*LOG(Z2)) * SIGG+MEAN
C
      RETURN
      END
C
C ------------------------------------------------------------------------------
      
      SUBROUTINE GENMUL(prop,TC,C,ISEED,WithRep)

C     This subroutine generates a multinomial random variable

      IMPLICIT NONE
C
C     Global variables
      REAL*8 prop(2000),C(2000)
      INTEGER ISEED,ISEED3,TC
      LOGICAL WithRep
C
C     Local variables
      REAL*8 tot,NTOT
      INTEGER A,A1,NP,DUM(100010),I,ICOUNT
      REAL*8 RAN1,RAN3
      EXTERNAL RAN1,RAN3
C
C     Generate a random number seed
      ISEED3 = -100000*RAN1(ISEED) 

C     Check if some catches were made!
      DO 9 a = 1,2000
       C(a) =  0.0
9     CONTINUE       
      IF (TC.eq.0) RETURN

C     First total the probabilities and then normalise
      tot = 0.0
      DO 10 a = 1,2000
       tot = tot + prop(a)
10    CONTINUE
      IF (Tot.EQ.0) RETURN
      DO 11 a = 1,2000
       prop(a) = prop(a)/tot
11    CONTINUE
C
C     Adjsut to handle "without-replacement" sampling
      NTOT = 100000.0
      IF (WithRep.AND.tot.LT.14000) NTOT = tot
C
C     Do special set up for index matrix
      NP = 0
      DO 13 a = 1,2000
       ICOUNT = NINT(prop(a)*NTOT)
       IF (ICOUNT.GT.0) THEN
        DO 14 I = 1,ICOUNT
         DUM(NP+I) = a
14      CONTINUE
        NP = NP + ICOUNT
       ENDIF
13    CONTINUE

C     Do actual generation
      DO 15 I = 1,TC
C
C      Generate an age
16     a1 = INT(RAN3(ISEED3)*NP+1)
       a = DUM(a1)
       IF (a.LE.0.AND.WithRep) GOTO 16
       C(a) = C(a) + 1
       IF (WithRep) DUM(a1) = -1
C
15    CONTINUE

      RETURN
      END
C      
C ------------------------------------------------------------------------
C
      REAL*8 FUNCTION STUDENT(INO,MEAN,SIGG,d,ISEED)
C
C     This subroutine generates variables from a t-distribution.
C
      REAL*8 V1,V2,E
      REAL*8 d,MEAN,SIGG
      INTEGER ISEED,INUM
      REAL*8 RAN1,RAN2,RAN3,RAN4
      EXTERNAL RAN1,RAN2,RAN3,RAN4
C
      INUM = 0
1     IF (INO.EQ.1) THEN
        V1 = 2.0*RAN1(ISEED)-1.0
        V2 = 2.0*RAN1(ISEED)-1.0
      ENDIF
      IF (INO.EQ.2) THEN
        V1 = 2.0*RAN2(ISEED)-1.0
        V2 = 2.0*RAN2(ISEED)-1.0
      ENDIF
      IF (INO.EQ.3) THEN
        V1 = 2.0*RAN3(ISEED)-1.0
        V2 = 2.0*RAN3(ISEED)-1.0
      ENDIF
      IF (INO.EQ.4) THEN
        V1 = 2.0*RAN4(ISEED)-1.0
        V2 = 2.0*RAN4(ISEED)-1.0
      ENDIF
      IF (V1.EQ.0) GOTO 1
      IF (V1**2.0+V2**2.0.GT.1.0) GOTO 1
      INUM = INUM +1
      Y = SQRT(d)*V2/V1
      E = (-d/2.0+0.5)*LOG(1.0+Y*Y/d)
c      E = (-d/2.0-0.5)*LOG(1.0+Y*Y/d)
      IF (LOG(RAN1(ISEED)+0.0000001).GT.E) GOTO 1
      STUDENT = Y*SIGG + MEAN
C
      RETURN
      END
C ------------------------------------------------------------------------------
C
      REAL*8 FUNCTION RAN1(ISEED)
C
C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 199).
C
C     Function returns a uniform random deviate between 0.0 & 1.0.
C     Set ISEED to any negative value to reinitialize the sequence.
C     Note: ISEED is reset to 1 on exit
C
      IMPLICIT REAL*8(M)
      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      SAVE MA,INEXT,INEXTP,IFF
C
      DIMENSION MA(55)
      DATA IFF /0/
      IF (ISEED.LT.0.OR.IFF.EQ.0) THEN
        IFF = 1
        MJ = MSEED-IABS(ISEED)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO 11 I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ = MA(II)
   11      CONTINUE
        DO 13 K = 1,4
          DO 12 I = 1,55
            MA(I) = MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
   12        CONTINUE
   13      CONTINUE
        INEXT = 0
        INEXTP = 31
        ISEED = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN1 = MJ*FAC
      IF (ABS(MJ*FAC).GT.1.01) THEN
        WRITE(*,*) MJ
        WRITE(*,*) FAC
        WRITE(*,*) MA
        WRITE(*,*) MJ*FAC
        WRITE(*,*) MA(INEXT)
        WRITE(*,*) MA(INEXTP)
        STOP
      ENDIF  
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      REAL*8 FUNCTION RAN2(ISEED)
C
C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 199).
C
C     Function returns a uniform random deviate between 0.0 & 1.0.
C     Set ISEED to any negative value to reinitialize the sequence.
C     Note: ISEED is reset to 1 on exit
C
      IMPLICIT REAL*8(M)
      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      SAVE MA,INEXT,INEXTP,IFF
C
      DIMENSION MA(55)
      DATA IFF /0/
      IF (ISEED.LT.0.OR.IFF.EQ.0) THEN
        IFF = 1
        MJ = MSEED-IABS(ISEED)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO 11 I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ = MA(II)
   11      CONTINUE
        DO 13 K = 1,4
          DO 12 I = 1,55
            MA(I) = MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
   12        CONTINUE
   13      CONTINUE
        INEXT = 0
        INEXTP = 31
        ISEED = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN2 = MJ*FAC
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      REAL*8 FUNCTION RAN3(ISEED)
C
C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 199).
C
C     Function returns a uniform random deviate between 0.0 & 1.0.
C     Set ISEED to any negative value to reinitialize the sequence.
C     Note: ISEED is reset to 1 on exit
C
      IMPLICIT REAL*8(M)
      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      SAVE MA,INEXT,INEXTP,IFF
C
      DIMENSION MA(55)
      DATA IFF /0/
      IF (ISEED.LT.0.OR.IFF.EQ.0) THEN
        IFF = 1
        MJ = MSEED-IABS(ISEED)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO 11 I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ = MA(II)
   11      CONTINUE
        DO 13 K = 1,4
          DO 12 I = 1,55
            MA(I) = MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
   12        CONTINUE
   13      CONTINUE
        INEXT = 0
        INEXTP = 31
        ISEED = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN3 = MJ*FAC
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      REAL*8 FUNCTION RAN4(ISEED)
C
C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 199).
C
C     Function returns a uniform random deviate between 0.0 & 1.0.
C     Set ISEED to any negative value to reinitialize the sequence.
C     Note: ISEED is reset to 1 on exit
C
      IMPLICIT REAL*8(M)
      PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      SAVE MA,INEXT,INEXTP,IFF
C
      DIMENSION MA(55)
      DATA IFF /0/
      IF (ISEED.LT.0.OR.IFF.EQ.0) THEN
        IFF = 1
        MJ = MSEED-IABS(ISEED)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO 11 I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ = MA(II)
   11      CONTINUE
        DO 13 K = 1,4
          DO 12 I = 1,55
            MA(I) = MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
   12        CONTINUE
   13      CONTINUE
        INEXT = 0
        INEXTP = 31
        ISEED = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN4 = MJ*FAC
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE GenMult(VEC,MEANS,df,ISEED,IDIST,NPARS,TT,TVAR,SG,
     +                   RESRES,PrMean,PrrMean,NumMean,RatLike,Icnter,
     +                   MPA,MPR)
C
C     Generate from a multivariant normal (Taken from MANTST, (C) C.
C     Allison (IWC). Set IDIST to 1 for normal and to 2 for t.
C     RESRES is the probability of having selected this option.
C
      IMPLICIT NONE
C
C     Global variables
      INTEGER MPR,MPA
      REAL*8 VEC(MPA),MEANS(MPR,500),df,RESRES,UNIF,RESI,RES1,RES2
      REAL*8 TT(MPA,MPA),TVAR(MPA,MPA),SG(MPA),PrMean(500),PrrMean(500)
      REAL*8 RatLike(500)
      INTEGER ISEED,IDIST,NPARS,NumMean,IMean,JMean,Icnter(500)
      REAL*8 XNORM,STUDENT,RAN1
      EXTERNAL XNORM,STUDENT,RAN1
C
C     Local variables
      REAL*8 AA(500),X,SIGS,Sigma(500),Temp,Rweight(500)
      INTEGER*4 I,J
C
C     First Find sigma 
      DO 1000 Imean = 1,NumMean
C
C      Find Sigma (needed for weighting options)
       Rweight(Imean) = 
     +    RatLike(1)/RatLike(Imean)*PrrMean(Imean)/PrrMean(1)
       Sigma(Imean) = Rweight(Imean)**(2.0/Npars)
C
1000  CONTINUE       
C
C     Select which mean to consider
      UNIF = RAN1(Iseed)
      DO 1700 Jmean = 1,NumMean
        IF (UNIF.LE.PrMean(Jmean)) GOTO 1701
1700  CONTINUE      
1701  CONTINUE
      Icnter(Jmean) = Icnter(Jmean) + 1
C
C     Initial errors to zero
      DO 1400 I = 1,NPARS
        AA(I) = 0.0
1400  CONTINUE
C
C     Generate the multivariante normal
      DO 1500 I = 1,NPARS
        SIGS = SG(I)*SQRT(Sigma(Jmean))
        IF (IDIST.EQ.1) X = XNORM(1,0.0d0,SIGS,ISEED)
        IF (IDIST.EQ.2) X = STUDENT(1,0.0d0,SIGS,df,ISEED)
        DO 1600 J = 1,NPARS
          AA(J) = AA(J) + X*TT(J,I)
1600    CONTINUE
1500  CONTINUE
C
C     Now complete the generation
      DO 1800 I = 1,NPARS
        VEC(I) = MEANS(I,Jmean) + AA(I)
1800  CONTINUE
C
C     Find the current (log) probability for this case
      RESRES = 0.0
      DO 2000 Imean = 1,NumMean
C
C      Calculate "central" term
       RESI = 0.0
       DO 2100 I = 1,NPARS
        RES1 = VEC(I) - Means(I,Imean)
        Temp = 0.0
        DO 2110 J = 1,NPARS
         RES2 = VEC(J) - Means(J,Imean)
         Temp = Temp + TVAR(I,J)*RES2
2110    CONTINUE         
        RESI = RESI + RES1*Temp
2100   CONTINUE
       RESI = RESI / Sigma(Imean)
       IF (IDIST.EQ.1) RESI = 0.5 * RESI
       IF (IDIST.EQ.2) RESI = ((df+NPARS)/2.0d0)*LOG(1.0d0+RESI/df)
C
C      Allow for different in varco width 
       Temp = 1.0/Rweight(Imean)
       RESRES = RESRES + PrrMean(Imean) * EXP(-RESI) * Temp
c       WRITE(*,701) Imean,PrMean(Imean),PrrMean(Imean),RESI,EXP(-RESI),
c     +     Temp,Sigma(Imean),PrrMean(Imean)*Temp*EXP(-RESI)
       
2000  CONTINUE       
      IF (RESRES.GT.0) THEN
        RESRES = -LOG(RESRES)
      ELSE
        RESRES = 1000
      ENDIF    
c      WRITE(*,702) JMean,RESRES,EXP(-RESRES),UNIF
C
      RETURN
701   FORMAT(1x,I3,2(1x,F7.5),1x,F10.3,1x,2(G10.5,1x),F7.3,1x,G10.5)
702   FORMAT(1x,I3,1x,F7.3,2(1x,G10.5))
      END
C
C -------------------------------------------------------------------------
C
      REAL*8 FUNCTION ZBQLGAM(G,H,ISEED)
*
*       Returns a random number with a gamma distribution with mean
*       G/H and variance G/(H^2). (ie. shape parameter G & scale
*       parameter H)
*
      REAL*8 C,D,R,G,H,A,z1,z2,B1,B2,M
      REAL*8U1,U2,U,V,TEST,X
      REAL*8 c1,c2,c3,c4,c5,w
      INTEGER ISEED,ISEED3
      REAL*8 RAN1,RAN3
      EXTERNAL RAN1,RAN3

C
C     Generate a random number seed
      ISEED3 = -100000*RAN1(ISEED) 

      ZBQLGAM = 0.0D0

      IF ( (G.LE.0.0D0).OR.(H.LT.0.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      IF (G.LT.1.0D0) THEN
889    u=RAN3(ISEED3)
       v=RAN3(ISEED3)
       if (u.gt.exp(1.0d0)/(g+exp(1.0d0))) goto 891
       ZBQLGAM=((g+exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g)
       if (v.gt.exp(-ZBQLGAM)) then
	goto 889
       else
	goto 892
       endif
891    ZBQLGAM=-log((g+exp(1.0d0))*(1.0d0-u)/(g*exp(1.0d0)))
       if (v.gt.ZBQLGAM**(g-1.0)) goto 889
892    ZBQLGAM=ZBQLGAM/h
       RETURN
      ELSEIF (G.LT.2.0D0) THEN
       M = 0.0D0
      elseif (g.gt.10.0d0) then
       c1=g-1.0d0
       c2=(g-1.0d0/(6.0d0*g))/c1
       c3=2.0d0/c1
       c4=c3+2.0d0
       c5=1.0d0/sqrt(g)
777    u=RAN3(ISEED3)
       v=RAN3(ISEED3)
       if (g.gt.2.50d0) then
	u=v+c5*(1.0d0-1.860d0*u)
       endif 
       if (u.le.0.0d0.or.u.ge.1.0d0) goto 777 
       w=c2*v/u 
       if (c3*u+w+1.0d0/w.le.c4) goto 778 
       if (c3*log(u)-log(w)+w.ge.1.0d0) goto 777
778    ZBQLGAM=c1*w/h 
       return
      ELSE
       M = -(G-2.0D0) 
      ENDIF
      R = 0.50D0
      a = ((g-1.0d0)/exp(1.0d0))**((g-1.0d0)/(r+1.0d0))
      C = (R*(M+G)+1.0D0)/(2.0D0*R)
      D = M*(R+1.0D0)/R
      z1 = C-DSQRT(C*C-D)
*
*     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
*     is not exactly zero - this needs trapping if negative.
*
      IF ((Z1-M.LT.0.0D0).AND.(Z1-M.GT.-1.0D-12)) Z1 = M
      z2 = C+DSQRT(C*C-D)
      B1=(z1*(z1-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z1-M)/(R+1.0D0))
      B2=(z2*(z2-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z2-M)/(R+1.0D0))
50    U1=RAN3(ISEED3)
      U2=RAN3(ISEED3)
      U=A*U1
      V=B1+(B2-B1)*U2
      X=V/(U**R)
      IF (X.LE.M) GOTO 50
      TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0D0))
      IF (U.LE.TEST) THEN
       ZBQLGAM = (X-M)/H
      ELSE
       GOTO 50
      ENDIF
 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLGAM',/5X, '(both parameters must be positive)',/)

      END
C
C -------------------------------------------------------------------------
C
      REAL*8 FUNCTION CUMNOM(X)

C     Work out the cumulative normal distribution

      REAL*8 X,Z,ERFCC
      EXTERNAL ERFCC

      Z = X/SQRT(2.0d0)

      IF (X.LT.0) THEN
        CUMNOM = ERFCC(-z)/2.0
      ELSE
        CUMNOM = 1-ERFCC(Z)/2.0
      ENDIF

      RETURN
      END
C
C -------------------------------------------------------------------------
C
      REAL*8 FUNCTION ERFCC(X)

C     Returns the complementry error function ERFC(X)

      REAL*8 Z,T,X

      Z = ABS(X)
      T = 1.0/(1.0+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     +     T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     +     T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0) ERFCC = 2.0-ERFCC

      RETURN
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE Calc_ag(fractOK,SigmaOK)
C
C     The purposes of this subroutine is to discover the sigma of a normal
C     distribution which results in a %agreement of fractOK. This is acheived
C     through the use of a bisection search approach.
C
      IMPLICIT NONE
C
C     Global varaibles
      REAL*8 fractOK,SigmaOK
C
C     Local variables
      REAL*8 Sigmax,Sigmin,Term1,Term2,Term3,Fracp
      INTEGER II
      EXTERNAL CUMNOM
      REAL*8 CUMNOM
C
C     Use bisection
      Sigmax = 5
      Sigmin = 0.0
      DO 1000 II = 1,50
C
C      Bisect interval
       SigmaOK = (Sigmax+Sigmin)/2.0
C
C      Evaluate function
       Term1 = CUMNOM( 0.5/SigmaOK) - CumNOM(-0.5/SigmaOK)
       Term2 = CUMNOM( 1.5/SigmaOK) - CumNOM( 0.5/SigmaOK)
       Term3 = CUMNOM( 2.5/SigmaOK) - CumNOM( 1.5/SigmaOK)
       Fracp = Term1*Term1 + 2.0*Term2*Term2 + 2.0*Term3*Term3
C
C      Check all is OK and adjust range
       IF (ABS(Fracp-FractOK).LT.0.000001) GOTO 1001
       IF (Fracp.GT.FractOK) THEN
         Sigmin = SigmaOK
       ELSE
         Sigmax = SigmaOK
       ENDIF
C
1000  CONTINUE
1001  CONTINUE
C
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE SORT(X,M)

C     USE A QUICK-SORT TO SORT ALL THE DATA

      IMPLICIT NONE
      REAL*8 X(21000),ST1(21000),MID
      INTEGER*4 M,LEFT(21000),RIGHT(21000),STKLEN
      INTEGER*4 LEFTS,RIGHTS,LS,RS,IC

C     Check for Daft call
      IF (M.LT.2) RETURN

C     Set up initial conditions
      LEFT(1) = 1
      RIGHT(1) = M
      STKLEN = 1

99    IF (STKLEN.EQ.0) GOTO 100

C     Set up the Pointers for this run
      MID = x(LEFT(STKLEN))
      LEFTS = LEFT(STKLEN)
      RIGHTS = RIGHT(STKLEN)
      LS = LEFT(STKLEN)
      RS = RIGHT(STKLEN)
                                      
C     Do a one-level sort
      DO 10 IC = LEFT(STKLEN)+1,RIGHT(STKLEN)

C      Check whether the current is less than the middle
       IF (X(IC).GT.MID) THEN
         ST1(RIGHTS) = X(IC)
         RIGHTS = RIGHTS - 1
       ELSE
         ST1(LEFTS) = X(IC)
         LEFTS = LEFTS + 1
       ENDIF
10    CONTINUE

C     Store the middle value
      ST1(LEFTS) = x(LEFT(STKLEN))

C     Replace the data
      DO 11 IC = LEFT(STKLEN),RIGHT(STKLEN)
       x(IC) = ST1(IC)
11    CONTINUE
      STKLEN = STKLEN - 1
        
C     update right pointer
      IF ((LEFTS-LS).GT.1) THEN
        STKLEN = STKLEN + 1
        LEFT(STKLEN) = LS
        RIGHT(STKLEN) = LEFTS - 1
      ENDIF
        
C     update left pointer
      IF ((RS-RIGHTS).GT.1) THEN
        STKLEN = STKLEN + 1
        LEFT(STKLEN) = RIGHTS + 1
        RIGHT(STKLEN) = RS
      ENDIF

      GOTO 99
100   CONTINUE

      RETURN
      END

C **************************************************************************
C *   MATRIX METHODS LIBRARY
C *   ======================
C *
C *
C *************************************************************************

      SUBROUTINE prmat(UNIT,X,d1,d2,n,m)

      IMPLICIT NONE
      INTEGER d1,d2,n,m,UNIT
      REAL*8 X(d1,d2)
      INTEGER I,J

      Do 10 I = 1,n
       WRITE(UNIT,500) (X(I,J),J=1,m)
500    FORMAT(1x,100(1x,F9.3))
10    CONTINUE

      RETURN
      END

C *************************************************************************

      SUBROUTINE TRANSM(X,XT,dim1,dim2,n,m)

C     Transpose the matrix X

      IMPLICIT NONE

      INTEGER dim1,dim2,n,m,I,J
      REAL*8 X(dim1,dim2),XT(dim2,dim1)

      DO 10 I = 1,n
       DO 10 J = 1,m 
        XT(J,I) = X(I,J)
10    CONTINUE

      RETURN
      END

C **************************************************************************

      SUBROUTINE MultM(A,B,C,d1,d2,d3,d4,d5,d6,n,m,p)

C     This subroutine multiple matricies A and B to give C
 
      IMPLICIT NONE
      INTEGER d1,d2,d3,d4,d5,d6,n,m,p,i,j,k
      REAL*8 A(d1,d2),B(d3,d4),C(d5,d6)

      DO 10 I = 1,n
       DO 10 J = 1,p
        C(I,J) = 0.0
        DO 10 K = 1,m
         C(I,J) = C(I,J) + A(I,K)*B(K,J)
10    CONTINUE

      RETURN
      END

C **************************************************************************

      SUBROUTINE SOLVEM(X,Y,d1,theta,n)

C     This subroutine attempts to solve a set of linear equations

      IMPLICIT NONE

C     Global variables
      INTEGER d1,n
      REAL*8 X(d1,d1),y(d1),theta(d1),XI(1000)

      CALL INVM(X,XI,d1,n)
      CALL MULTM(XI,Y,THETA,d1,d1,d1,1,d1,1,n,n,1)
      
      RETURN
      END

C **************************************************************************

      SUBROUTINE SOLVEL(X,Y,d1,d2,theta,n,m)

C     This subroutine performs a multiple linear regression

      IMPLICIT NONE

C     Global variables
      INTEGER d1,d2,n,m
      REAL*8 X(d1,d2),y(d1),theta(d2)

C     Local variables
      REAL*8 XT(1000),C(1000),CI(1000)

      CALL TRANSM(X,XT,d1,d2,n,m)
      CALL MULTM(XT,X,C,d2,d1,d1,d2,d2,d2,m,n,m)
      CALL INVM(C,CI,d2,m)
      CALL MULTM(CI,XT,C,d2,d2,d2,d1,d2,d1,m,m,n)
      CALL MULTM(C,Y,THETA,d2,d1,d1,1,d2,1,m,n,1)

      RETURN
      END

C **************************************************************************

      SUBROUTINE INVM(A,AI,dim1,n)

C     This subroutine computes the inverse of a matrix

      IMPLICIT NONE

C     Global variables
      INTEGER dim1,n
      REAL*8 A(dim1,dim1),AI(dim1,dim1)

C     Local variables
      REAL*8 D
      INTEGER INDEX(500),I,J

C     Set up an identity matrix
      DO 12 I = 1,n
        DO 11 J = 1,n
          AI(I,J) = 0.0
11      CONTINUE
        AI(I,I) = 1.0
12    CONTINUE

C     Perform a LU-decomposition
      CALL LUDCMP(A,n,dim1,INDEX,D)

C     Find the inverse
      DO 13 J = 1,N
        CALL LUBKSB(A,n,dim1,INDEX,AI(1,J))
13    CONTINUE

      RETURN
      END

C **************************************************************************

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      IMPLICIT NONE
      INTEGER N,NP,INDX(N)
      REAL*8 A(NP,NP),B(N),SUM
      INTEGER II,I,LL,J

      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END

C ************************************************************************
C      
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      IMPLICIT NONE
      
      INTEGER NMAX
      REAL*8 TINY
      PARAMETER (NMAX=500,TINY=1.0D-50)

      INTEGER N,NP,INDX(N)
      REAL*8 A(NP,NP),VV(NMAX),SUM,AAMAX,DUM,D
      
      INTEGER I,J,K,IMAX

      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) WRITE(*,*) 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END

C ***********************************************************************

      SUBROUTINE JACOBI(A,N,NP,D,V)

C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 346).

C     Computes all eigenvalues and eigenvectors of a real symmetric matrix
C     A, which is of a size N by N, stored in a physical NP by NP array. On
C     output, elements of A above the diagonal are destroyed. D returns the
C     eigenvalues of A in its first N elements. V is a matrix with the same
C     logical and physical dimensions as A whose columns contain, on
C     output, the normalized eigenvectors of A. NROT returns the number of
C     Jacobi rotations which were required.

      IMPLICIT NONE
      INTEGER N,NP,IP,IQ,NROT,I,J,NMAX
      PARAMETER (NMAX = 500)
      REAL*8 A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX),SM,TRESH,G,H,T,C,
     +     S,TAU,THETA

      DO 12 IP = 1,N
        DO 11 IQ = 1,N
          V(IP,IQ) = 0.
   11   CONTINUE
        V(IP,IP) = 1.
   12 CONTINUE
      DO 13 IP = 1,N
        B(IP) = A(IP,IP)
        D(IP) = B(IP)
        Z(IP) = 0.
   13 CONTINUE
      NROT = 0
      DO 24 I = 1,50
        SM = 0.
        DO 14 IP = 1,N-1
        DO 14 IQ = IP+1,N
          SM = SM+ABS(A(IP,IQ))
 14     CONTINUE
        IF (SM.EQ.0.)RETURN
        IF (I.LT.4) THEN
          TRESH = 0.2*SM/N**2
        ELSE
          TRESH = 0.
        ENDIF
        DO 22 IP = 1,N-1
          DO 21 IQ = IP+1,N
            G = 100.*ABS(A(IP,IQ))
C           After 4 sweeps, skip rotation if off-diagonal element is small
            IF ((I.GT.4) .AND. (ABS(D(IP))+G .EQ. ABS(D(IP)))
     +                   .AND. (ABS(D(IQ))+G .EQ. ABS(D(IQ)))) THEN
              A(IP,IQ) = 0.0
            ELSE IF (ABS(A(IP,IQ)).GT.TRESH) THEN
              H = D(IQ)-D(IP)
              IF (ABS(H)+G.EQ.ABS(H)) THEN
                T = A(IP,IQ)/H
              ELSE
                THETA = 0.5*H/A(IP,IQ)
                T = 1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF (THETA.LT.0.)T = -T
              ENDIF
              C = 1./SQRT(1+T**2)
              S = T*C
              TAU = S/(1.+C)
              H = T*A(IP,IQ)
              Z(IP) = Z(IP)-H
              Z(IQ) = Z(IQ)+H
              D(IP) = D(IP)-H
              D(IQ) = D(IQ)+H
              A(IP,IQ) = 0.
              DO 16 J = 1,IP-1
                G = A(J,IP)
                H = A(J,IQ)
                A(J,IP) = G-S*(H+G*TAU)
                A(J,IQ) = H+S*(G-H*TAU)
   16         CONTINUE
              DO 17 J = IP+1,IQ-1
                G = A(IP,J)
                H = A(J,IQ)
                A(IP,J) = G-S*(H+G*TAU)
                A(J,IQ) = H+S*(G-H*TAU)
   17         CONTINUE
              DO 18 J = IQ+1,N
                G = A(IP,J)
                H = A(IQ,J)
                A(IP,J) = G-S*(H+G*TAU)
                A(IQ,J) = H+S*(G-H*TAU)
   18         CONTINUE
              DO 19 J = 1,N
                G = V(J,IP)
                H = V(J,IQ)
                V(J,IP) = G-S*(H+G*TAU)
                V(J,IQ) = H+S*(G-H*TAU)
   19         CONTINUE
              NROT = NROT+1
            ENDIF
   21     CONTINUE
   22   CONTINUE
        DO 23 IP = 1,N
          B(IP) = B(IP)+Z(IP)
          D(IP) = B(IP)
          Z(IP) = 0.
   23   CONTINUE
   24 CONTINUE
      WRITE (*,'(A)') '50 iterations should never happen'
      RETURN
      END

C  ***********************************************************************

      SUBROUTINE ZBRAC(FUNC,PARS,M,X1,X2,F1,F2,SUCCES)
C
C  This subroutine is based on the Numerical Recipes Software Library
C  (Refer, W.H.Press, B.P.Flannery, S.A.Teukolsky & W.T.Vetterling 1986.  
C  Numerical Recipes: The Art of Scientific Computing.  CUP, Cambridge. 818pp.
C
C  Given a function FUNC and an initial guessed range X1 to X2, the
C  routine expands the range geometrically until a root is bracketted
C  by the return values X1 and X2 (in which case SUCCES returns as
C  .TRUE.) or until the range becomes unacceptably large (in which
C  case SUCCES returns as .FALSE.).  The parameters to FUNC are passed
C  through the vector PARS.  X1 and X2 are passed to FUNC through
C  PARS(M).  The corresponding function values F1 and F2 are returned.
C  Success is guaranteed for a function which has opposite sign for
C  sufficiently large and small arguments.
C
      EXTERNAL FUNC
      REAL*8 FUNC,PARS(4),FACTOR,X1,X2,F1,F2
      INTEGER NTRY,M,IERR,J
      LOGICAL SUCCES
      PARAMETER (FACTOR=1.6D0,NTRY=50)
C
      IF(X1.EQ.X2) WRITE(*,*) 'You have to guess an initial range'
      PARS(M) = X1
      F1 = FUNC(PARS,IERR)
      IF (IERR .NE. 0) THEN
        SUCCES = .FALSE.
        RETURN
      ENDIF
      PARS(M) = X2
      F2 = FUNC(PARS,IERR)
      IF (IERR .NE. 0) THEN
        SUCCES = .FALSE.
        RETURN
      ENDIF
      SUCCES=.TRUE.
      DO 11 J=1,NTRY
        IF(F1*F2.LT.0.D0)RETURN
        IF(ABS(F1).LT.ABS(F2))THEN
          X1=X1+FACTOR*(X1-X2)
          PARS(M) = X1
          F1 = FUNC(PARS,IERR)
          IF (IERR .NE. 0) THEN
            SUCCES = .FALSE.
            RETURN
          ENDIF
        ELSE
          X2=X2+FACTOR*(X2-X1)
          PARS(M) = X2
          F2 = FUNC(PARS,IERR)
          IF (IERR .NE. 0) THEN
            SUCCES = .FALSE.
            RETURN
          ENDIF
        ENDIF
11    CONTINUE
      SUCCES=.FALSE.
      RETURN
      END
C
C  ***********************************************************************
C
      SUBROUTINE ZBRENT(FUNC,PARS,M,X1,X2,FX1,FX2,TOL,IERR)

C  This subroutine is based on the Numerical Recipes Software Library
C  (Refer, W.H.Press, B.P.Flannery, S.A.Teukolsky & W.T.Vetterling 1986.  
C  Numerical Recipes: The Art of Scientific Computing.  CUP, Cambridge. 818pp.
C
C  Using Brent's method, find the root of a function FUNC known to
C  lie between X1 and X2 (with function values FX1 and FX2).  The
C  parameters to FUNC are passed through the vector PARS.  The root
C  is returned through PARS(M)
C
      EXTERNAL FUNC
      REAL*8 FUNC,PARS(4),X1,X2,FX1,FX2,TOL,EPS,A,FA,B,FB,C,D,E,TOL1,XM,
     +       S,P,Q,R,FC
      INTEGER M,ITMAX,ITER,IERR
      PARAMETER (ITMAX=100,EPS=3.E-8)
C
      A=X1
      FA = FX1
      B=X2
      FB = FX2
C
c      IF(FB*FA.GT.0.D0) WRITE(*,*) 'Root must be bracketed for ZBRENT.'
      FC=FB
      DO 11 ITER=1,ITMAX
C
C        IF(FB/FC.GT.0.D0) THEN   *******************  AEP
C        IF(FB*FC.GT.0.D0) THEN   *******************  DLM
        IF(FB/FC.GT.0.D0) THEN
C  Rename A, B, C and adjust bounding interval D
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
C
C  Check for convergence
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.D0*EPS*ABS(B)+0.5D0*TOL
        XM=.5D0*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0)THEN
          PARS(M)=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
C  Attempt inverse quadratic interpolation
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
C
C  Check whether in bounds
          IF(P.GT.0.D0) Q=-Q
          P=ABS(P)
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
C  Accept interpolation
            E=D
            D=P/Q
          ELSE
C  Interpolation failed, use bisection
            D=XM
            E=D
          ENDIF
        ELSE
C  Bounds decreasing too slowly, use bisection
          D=XM
          E=D
        ENDIF
C  Save latest best guess
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
C  Evaluate new trial root
        PARS(M) = B
        FB = FUNC(PARS,IERR)
        IF (IERR .NE. 0) RETURN
11    CONTINUE
C
      WRITE (*,*) 'ZBRENT exceeding maximum iterations.'
      PARS(M)=B
      RETURN
      END
C
C ===========================================================================
C
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      
      IMPLICIT NONE
      
      INTEGER NMAX
      PARAMETER (NMAX=100)
      
      INTEGER N,I,K
      REAL*8 X(N),Y(N),Y2(N),U(NMAX),YP1,YPN
      REAL*8 SIG,P,QN,UN
      
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
      
      INTEGER N,K,KLO,KHI
      REAL*8 XA(N),YA(N),Y2A(N),X,Y
      REAL*8 H,A,B
      
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) WRITE (*,*) 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
C
C ===========================================================================
C
      SUBROUTINE FIT(X,SS,NDIM)

C     Set up the parameters for a fit

C     IMPLICIT NONE
      REAL*8 P(101,100),Y(101),X(100),SS,FUNK,TOL,GRD
      INTEGER*4 NDIM,I,J,ITER
      EXTERNAL FUNK

C     SET UP TOLERANCES AND GRIDDING
      TOL = 0.0000000001
      GRD = 1.2

      DO 1 I=1,NDIM+1
         DO 2 J=1,NDIM
            P(I,J)=X(J)
            IF((I-1).EQ.J) P(I,J) = GRD*P(I,J)
2        CONTINUE
1     CONTINUE

      DO 3 I=1,NDIM+1
         DO 4 J=1,NDIM
            X(J)=P(I,J)
4        CONTINUE
         Y(I)=FUNK(X)
3     CONTINUE

      CALL AMOEBA(P,Y,101,100,NDIM,TOL,ITER)

      DO 5 J = 1,NDIM
         X(J) = P(1,J)
5     CONTINUE
      SS = Y(1)

      RETURN
      END

      SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,ITER)
C
C     MULTIDIMENSIONAL MINIMISATION OF THE FUNCTION FUNK(X) WHERE X IS
C     AN NDIM-DIMENSIONAL VECTOR, BY THE DOWNHILL SIMPLEX METHOD OF
C     NELDER AND MEAD. INPUT IS A MATRIX P WHOSE NDIM+1 ROWS ARE THE
C     NDIM-DIMENSIONAL VECTORS WHICH ARE THE VERTICES OF THE STARTING
C     SIMPLEX. [LOGICAL DIMENSIONS OF P ARE P(NDIM+1,NDIM); PHYSICAL
C     DIMENSIONS ARE INPUT AS P(MP,NP).] ALSO INPUT IS THE VECTOR Y
C     OF LENGTH NDIM+1, WHOSE COMPONENTS MUST BE PRE-INITIALISED TO
C     THE VALUES OF FUNK EVALUATED AT THE NDIM+1 VERTICES (ROWS) OF P;
C     AND FTOL IS THE FRACTIONAL CONVERGENCE TOLERANCE TO BE ACHIEVED
C     IN THE FUNCTION VALUE (N.B.!). ON OUTPUT, P AND Y WILL HAVE BEEN
C     RESET TO NDIM+1 NEW POINTS ALL WITHIN FTOL OF A MINIMUM FUNCTION
C     VALUE, AND ITER GIVES THE NUMBER OF ITERATIONS TAKEN.
C
C     FROM: NUMERICAL RECIPES - THE ART OF SCIENTIFIC COMPUTING
C           BY W. H. PRESS ET AL, CAMBRIDGE UNIVERSITY PRESS
C           ISBN 0-251-30811-9
C
C     ********************************************************************
C

C     SPECIFY THE MAXIMUM NUMBER OF DIMENSIONS, THREE PARAMETERS WHICH
C     DEFINE THE EXPANSIONS AND CONTRACTIONS, AND THE MAXIMUM NUMBER OF
C     ITERATIONS ALLOWED

C     IMPLICIT NONE
      REAL*8 ALPHA,BETA,GAMMA,EPS
      INTEGER*4 ITMAX
      PARAMETER (ALPHA=1.0,BETA=0.5,GAMMA=2.0,ITMAX=5000,EPS=1.0E-10)

C     Global Data
      INTEGER*4 MP,NP,NDIM,ITER
      REAL*8 P(MP,NP),Y(MP),PR(100),PRR(100),PBAR(100),FTOL,FUNK
      EXTERNAL FUNK

C     Local Data
      REAL*8 YPR,YPRR
      INTEGER*4 I,ILO,IHI,INHI,J,MPTS

C     NOTE THAT MP IS THE PHYSICAL DIMENSION CORRESPONDING TO THE LOGICAL
C     DIMENSION MPTS, NP TO NDIM.

      MPTS = NDIM+1
      ITER=0

C     FIRST WE MUST DETERMINE WHICH POINT IS THE HIGHEST (WORST), NEXT
C     HIGHEST, AND LOWEST (BEST).

1     ILO=1
      IF(Y(1).GT.Y(2)) THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      ENDIF
      DO 11 I=1,MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI)) THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI)) THEN
            IF(I.NE.IHI) INHI=I
         ENDIF
11    CONTINUE

C     COMPUTE THE FRACTIONAL RANGE FROM THE HIGHEST TO THE LOWEST AND
C     RETURN IF SATISFACTORY

      IF(2.0*ABS(Y(IHI)-Y(ILO)).LT.FTOL*(ABS(Y(IHI))+ABS(Y(ILO))+EPS)) 
     +   RETURN
      IF(ITER.EQ.ITMAX) THEN
c        WRITE(6,200)
200      FORMAT(1H ,'AMOEBA EXCEEDING MAXIMUM ITERATIONS')
         RETURN
      ENDIF
      ITER=ITER+1

      DO 12 J=1,NDIM
         PBAR(J)=0.
12    CONTINUE

C     BEGIN A NEW ITERATION. COMPUTE THE VECTOR AVERAGE OF ALL POINTS
C     EXCEPT THE HIGHEST, I.E. THE CENTRE OF THE "FACE" OF THE SIMPLEX
C     ACROSS FROM THE HIGH POINT. WE WILL SUBSEQUENTLY EXPLORE ALONG
C     THE RAY FROM THE HIGH POINT THROUGH THE CENTRE.

      DO 14 I=1,MPTS
         IF(I.NE.IHI) THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(I,J)
13          CONTINUE
         ENDIF
14    CONTINUE

C     EXTRAPOLATE BY A FACTOR ALPHA THROUGH THE FACE, I.E. REFLECT THE
C     SIMPLEX FROM THE HIGH POINT

      DO 15 J=1,NDIM
         PBAR(J)=PBAR(J)/NDIM
         PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE

C     EVALUATE THE FUNCTION AT THE REFLECTED POINT

      YPR=FUNK(PR)

C     GIVES A BETTER RESULT THAN THE BEST POINT, SO TRY AN ADDITIONAL
C     EXTRAPOLATION BY A FACTOR GAMMA

      IF(YPR.LE.Y(ILO)) THEN
         DO 16 J=1,NDIM
            PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16       CONTINUE

C        CHECK THE FUNCTION THERE

         YPRR=FUNK(PRR)

C        THE ADDITIONAL EXTRAPOLATION SUCCEEDED, AND REPLACES THE
C        HIGHEST POINT

         IF(YPRR.LT.Y(ILO)) THEN
            DO 17 J=1,NDIM
               P(IHI,J)=PRR(J)
17          CONTINUE
            Y(IHI)=YPRR
         ELSE

C        THE ADDITIONAL EXTRAPOLATION FAILED, BUT WE CAN STILL USE THE
C        REFLECTED POINT

         DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18       CONTINUE
         Y(IHI)=YPR
      ENDIF

C     THE REFLECTED POINT IS WORSE THAN THE SECOND HIGHEST

      ELSE IF(YPR.GE.Y(INHI)) THEN

C        IF IT'S BETTER THAN THE HIGHEST, THEN REPLACE THE HIGHEST

         IF(YPR.LT.Y(IHI)) THEN
            DO 19 J=1,NDIM
               P(IHI,J)=PR(J)
19          CONTINUE
            Y(IHI)=YPR
         ENDIF

C        BUT LOOK FOR AN INTERMEDIATE LOWER POINT; IN OTHER WORDS
C        PERFORM A CONTRACTION OF THE SIMPLEX ALONG ONE DIMENSION
C        AND THEN EVALUATE THE FUNCTION

         DO 21 J=1,NDIM
            PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21       CONTINUE
         YPRR=FUNK(PRR)

C        CONTRACTION GIVES AN IMPROVEMENT, SO ACCEPT IT

         IF(YPRR.LT.Y(IHI)) THEN
            DO 22 J=1,NDIM
               P(IHI,J)=PRR(J)
22          CONTINUE
            Y(IHI)=YPRR
         ELSE

C           CAN'T SEEM TO GET RID OF THAT HIGH POINT. BETTER CONTRACT
C           AROUND THE LOWEST (BEST) POINT

            DO 24 I=1,MPTS
               IF(I.NE.ILO) THEN
                  DO 23 J=1,NDIM
                     PR(J)=0.5*(P(I,J)+P(ILO,J))
                     P(I,J)=PR(J)
23                CONTINUE
                  Y(I)=FUNK(PR)
               ENDIF
24          CONTINUE
         ENDIF

      ELSE

C        WE ARRIVE HERE IF THE ORIGINAL REFLECTION GIVES A MIDDLING
C        POINT. REPLACE THE OLD HIGH POINT AND CONTINUE

         DO 25 J=1,NDIM
            P(IHI,J)=PR(J)
25       CONTINUE
         Y(IHI)=YPR

      ENDIF

      GO TO 1

      END
C
C ===========================================================================
C
      REAL*8 FUNCTION FUNK(X)
C
C This is a dummy function to be used with FIT
C     
      REAL*8 X(100)
C
      FUNK = X(1)
C      
      RETURN
      END
C
C ===========================================================================
C
      
