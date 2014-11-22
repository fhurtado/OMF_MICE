C Pacific sardine simulation code. Written by F. Hurtado-Ferro and A.E. Punt
C Please reference the authors if you use this code.
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Local variables 
      INTEGER Year,Isim,Ispec,Istk,Jspec
C      
C     Read in the specifications
      CALL ReadInBiol()
C
C     Read in the parameters WRIof the control rule
C
      OPEN(UNIT=14,FILE="SUMMARY.OUT")
C	
C     Generate the underlying environmental signal (common for all simulations)  
      CALL CalcG()
C
      ISEED1 = ISEED1Init
      ISEED2 = ISEED2Init
      ISEED3 = ISEED3Init
      ISEED4 = ISEED4Init
      DO 10000 Isim = 1,Nsim
       WRITE(*,*) Isim
C
C      Initialize prey
       CALL InitialPrey(Isim)
C
C      Project prey
       DO 11000 Year = 1,Nyear
        CALL Proj1YrPrey(Isim,Year,.TRUE.)
11000  CONTINUE 
C
10000 CONTINUE   
      WRITE(*,*) "Unfished Simulation finished running"

      B0 = 0
      DO 30000 Isim = 1,Nsim
       DO 30001 Year = 51,Nyear
        DO 30002 Ispec = 1,Nspec
         DO 30003 Istk = 1,Nstk(Ispec)
          B0(Ispec,Istk) = B0(Ispec,Istk) + B1p(Ispec,Istk,Isim,Year)
30003 CONTINUE
30002 CONTINUE
30001 CONTINUE
30000 CONTINUE
      DO 31000 Ispec = 1,Nspec
       DO 31000 Istk = 1,Nstk(Ispec)
         B0(Ispec,Istk) = B0(Ispec,Istk) / DFLOAT(Nsim*(Nyear-50))
         WRITE(*,*) B0(Ispec,Istk)
31000 CONTINUE
      DO 32000 Jspec = 1,Npred
       B0Pred(Jspec)= 0
       DO 32100 Ispec = 1,Nspec
        DO 32100 Istk = 1,Nstk(Ispec)
         B0Pred(Jspec) = B0Pred(Jspec) + 
     +      B0(Ispec,Istk)*PredWeight(Jspec,Ispec,Istk)
32100  CONTINUE
32000 CONTINUE
C     
      ISEED1 = ISEED1Init
      ISEED2 = ISEED2Init
      ISEED3 = ISEED3Init
      ISEED4 = ISEED4Init
      DO 20000 Isim = 1,Nsim
       WRITE(*,*) Isim
C
C      Initialize prey
       CALL InitialPrey(Isim)
C
C      Project prey
       DO 21000 Year = 1,Nyear
        CALL Proj1YrPrey(Isim,Year,.FALSE.)
21000  CONTINUE 
C
C      Initialize pred
       CALL InitialPred(Isim)
C
C      Project predators
       DO 22000 Year = 1,Nyear
        CALL Proj1YrPred(Isim,Year)
22000  CONTINUE 
C
20000 CONTINUE   
      WRITE(*,*) "Simulation finished running"
C      
      WRITE(*,*) "Writing summary output"  
      CALL SUMMARY()
      CLOSE(14)
C
      STOP
      END
C
C ===============================================================================================
C
      SUBROUTINE ReadInBiol()
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Local variables
      INTEGER Age,Ispec,Istk,Jspec,TEST
C
      OPEN(UNIT=13,FILE="Sardine.Dat")
      READ(13,*)
      READ(13,*) Nspec
      READ(13,*)
      READ(13,*) (Nstk(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,*) (Nage(Ispec),Ispec=1,Nspec)
      WRITE(*,'(A5,1x,5I3)') "Nage",Nage
      READ(13,*)
      DO 10100 Ispec = 1,Nspec
       READ(13,*) (Weight(Ispec,Age),Age=0,Nage(Ispec))
10100 CONTINUE       
      READ(13,*)
      READ(13,*) (Matu(Age),Age=0,Nage(1))
      READ(13,*)
      DO 10200 Ispec=1,Nspec
       DO 10200 Istk = 1,Nstk(Ispec)
       READ(13,*) (Ninit(Ispec,Istk,Age),Age=0,Nage(Ispec))
10200 CONTINUE       
      READ(13,*)
      READ(13,*) (Selex(Age),Age=0,Nage(1))
      READ(13,*)
      READ(13,*) (SelxMC1(Age),Age=0,Nage(1))
      READ(13,*)
      READ(13,*) (SelxMC2(Age),Age=0,Nage(1))
      READ(13,*)
      READ(13,*) (SelxPNW(Age),Age=0,Nage(1))
      READ(13,*)
      READ(13,*) (Natm(Ispec),Ispec=1,Nspec)
      READ(13,*)
      READ(13,*) AlogRPS, BlogRPS, PhilogRPS
      READ(13,*)
      READ(13,*) Alpha, Beta1, Beta2
      READ(13,*)
      DO 21100 Ispec = 1,Nspec
       READ(13,*) (EpsInit(Ispec,Istk),Istk=1,Nstk(Ispec))
21100 CONTINUE
      READ(13,*)
      DO 21200 Ispec = 1,Nspec
       READ(13,*) ProwR(Ispec),SigmaR(Ispec)
21200 CONTINUE   
      READ(13,*)
      READ(13,*) Rerr,Rtype
      READ(13,*)
      READ(13,*) Vinit,Gtype,UseG,Verr
      READ(13,*)
      READ(13,*) Ybar,Period,ProwV,SigmaV,Phi,MeanSST, Amplitude
      READ(13,*)
      READ(13,*) SigmaB, Berr
      READ(13,*)
      READ(13,*) Chi
      READ(13,*)
      READ(13,*) Power
      READ(13,*)
      READ(13,*) SigmaI,Ierr
      CLOSE(13)
C
      OPEN(UNIT=13,FILE="Sardine.Spec")
      READ(13,*)
      READ(13,*) Nyear
      READ(13,*)
      READ(13,*) Nsim
      READ(13,*) 
      READ(13,*) Emsymin, Emsymax	  
      READ(13,*) 
      READ(13,*) FracMin, FracMax
      READ(13,*) 
      READ(13,*) Cutoff
      READ(13,*) 
      READ(13,*) Maxcatch
      READ(13,*) 
      READ(13,*) OFLopt, NoOFL
      READ(13,*) 
      READ(13,*) HCRopt
      READ(13,*) 
      READ(13,*) CutOpt
      READ(13,*)
      READ(13,*) Pstar,OFLmax
      READ(13,*)
      READ(13,*) ISEED1Init, ISEED2Init, ISEED3Init, ISEED4Init
      READ(13,*) 
      READ(13,*) TVselx
      READ(13,*) 
      READ(13,*) MultFleet
      READ(13,*) 
      READ(13,*) Fmex,Fcan
      READ(13,*) 
      READ(13,*) HyperSt, HPSa
      READ(13,*) 
      READ(13,*) TVM
      READ(13,*) 
      READ(13,*) SelxOut
      READ(13,*) 
      READ(13,*) UseAnn	  
      READ(13,*) 
      READ(13,*) DiffSelx	  
      CLOSE(13)
C
      OPEN(UNIT=13,FILE="Predator.DAT")
      READ(13,*)
      READ(13,*) NPred
      READ(13,*)
      READ(13,*) PredImpact
      READ(13,*)
      READ(13,*) (PNage(Ispec),Ispec=1,NPred)
      READ(13,*)
      READ(13,*) (PMat(Ispec),Ispec=1,NPred)
      READ(13,*)
      READ(13,*) (PFmax(Ispec),Ispec=1,NPred)
      READ(13,*)
      DO 41000 ISpec=1,NPred
       READ(13,*) (SP(Ispec,Age),Age=0,PNAge(Ispec))
41000 CONTINUE
      READ(13,*)
      READ(13,*) ZZ
      READ(13,*)
      DO 42000 Ispec = 1,Npred
       Read(13,*)((PredWeight(Ispec,Jspec,Istk),
     +     Istk=1,Nstk(Jspec)),Jspec=1,Nspec)
42000 CONTINUE      
      READ(13,*)
      DO 43000 Ispec = 1,Npred
       READ(13,*) PX(Ispec),PZ(Ispec),Psteep(Ispec)
43000 CONTINUE       
      READ(13,*)
      READ(13,*) TEST
      IF (TEST.NE.123456) THEN
       WRITE(*,*) "TEST IS WRONG",TEST
       STOP
      ENDIF 
      CLOSE(13)

      CALL SolveFec()
C      
      RETURN
      END
C
C ===============================================================================================
C
      SUBROUTINE SolveFec
C
C This subroutine solves for the fecundity / age-0 survival multiplier
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Local variables
      INTEGER Ispec,II
      REAL*8 SJMIN,SJMAX,SA,TARGET,MULT,SJ,PRED
C
      DO 10000 Ispec = 1,Npred
       SJmin = 0
       SJmax = 10
       SA = SP(Ispec,1)
       Target = 1.0
       DO 11000 II = 1,50
        Mult = (SJMin+SJMax)/2.0
        SJ = SP(Ispec,0)*Mult
        Pred = SA + 
     +     SJ*SA**(Pmat(Ispec)-1)
        IF (Pred.LT.Target) THEN
         SJMin = Mult
        ELSE
         SJmax = Mult
        ENDIF
11000  CONTINUE        
       SP(Ispec,0) = SP(Ispec,0)*Mult
10000 CONTINUE      
C
      DO 20000 Ispec = 1,Npred
       SJmin = 0
       SJmax = 10
       SA = SP(Ispec,1)
       Target = (1+PFMax(Ispec))**Pmat(Ispec)
       DO 21000 II = 1,50
        Mult = (SJMin+SJMax)/2.0
        SJ = SP(Ispec,0)*Mult
        Pred = (1+PFMax(Ispec))**(Pmat(Ispec)-1)*SA + 
     +     SJ*SA**(Pmat(Ispec)-1)
        IF (Pred.LT.Target) THEN
         SJMin = Mult
        ELSE
         SJmax = Mult
        ENDIF
21000  CONTINUE        
       MaxSJ(Ispec) = Mult
20000 CONTINUE      
C      
      RETURN
      END
C
C ===============================================================================================
C
      SUBROUTINE InitialPrey(Isim)
C
C This subroutine sets up variables which are common across a run
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Global variables
      INTEGER Isim
C
C     Local variables
      INTEGER Iage,Year,Ispec,Istk
      REAL*8 EpsV,EpsR,SSBV,Bio1P,Idev
      REAL*8 ANcount1, ANcount2
      REAL*8 Lselx,Jselx,MaxSelx
      REAL*8 Gslope
      REAL*8 XNORM
      EXTERNAL XNORM

C     Initial N-vector
      DO 10000 Ispec = 1,Nspec
       DO 10000 Istk = 1,Nstk(Ispec)
        DO 10100 Iage = 0,Nage(Ispec)
         N(Ispec,Istk,0,Iage) = Ninit(Ispec,Istk,Iage)
10100   CONTINUE         
        RecRep(Ispec,Istk,Isim,0) = N(Ispec,Istk,0,0)       
10000 CONTINUE
C
C     Mean Age
      ANcount1 = 0      
      ANcount2 = 0
      AgeN(Isim,0) = 0
      DO 11000 Iage = 0,Nage(1)
       ANcount1 = ANcount1 + N(1,1,0,Iage)*Iage
       ANcount2 = ANcount2 + N(1,1,0,Iage)
11000 CONTINUE
      AgeN(Isim,0) = ANcount1/ANcount2
C
C     Compute the SSB and 1+ biomass
      DO 15000 Ispec = 1,Nspec
       DO 15000 Istk = 1, Nstk(Ispec)
        SSBV = 0
        Bio1P = 0
        DO 15100 Iage = 1,Nage(Ispec)
         Bio1P = Bio1P + Weight(Ispec,Iage)*N(Ispec,Istk,0,Iage)
15100   CONTINUE
        DO 15200 Iage = 2,Nage(Ispec)
         SSBV = SSBV + Weight(Ispec,Iage)*N(Ispec,Istk,0,Iage)
15200   CONTINUE
        SSB(Ispec,Istk,Isim,0) = SSBV
        B1p(Ispec,Istk,Isim,0) = Bio1P
15000 CONTINUE        
C
C     Generate the environmental variables
      Vdev(0) = Vinit
      V(0) = Vdev(0) + MeanSST
      Vrep(Isim,0) = V(0)
C	  
      DO 31000 Year = 1,Nyear
C
C      Update
       EpsV = XNORM(1,0.0d0,1.0d0,ISEED1)*SigmaV
       EpsV = EpsV * Verr
       Vdev(Year) = ProwV*Vdev(Year-1)+(1-ProwV)*G(Year)*Amplitude
       Vdev(Year) = Vdev(Year) + SQRT(1-ProwV**2)*EpsV
       V(Year) = Vdev(Year)*UseG + MeanSST
       Vrep(Isim,Year) = V(Year)
C
31000 CONTINUE
C
C     Generate the rec_dev
      DO 32000 Ispec = 1,Nspec
       DO 32000 Istk = 1,Nstk(Ispec)
        Eps(Ispec,Istk,0) = EpsInit(Ispec,Istk)
        DO 32100 Year = 1,Nyear
         EpsR = XNORM(2,0.0d0,1.0d0,ISEED2)*SigmaR(Ispec)
         EpsR = EpsR * Rerr
         Eps(Ispec,istk,Year) = prowR(Ispec)*Eps(Ispec,Istk,Year-1)+
     +      SQRT(1-ProwR(Ispec)**2.0)*EpsR
32100  CONTINUE
32000 CONTINUE
C      
C     Add measurement error to the environmental variable
      DO 33000 Year = 0,Nyear
       Idev = XNORM(1,0.0d0,1.0d0,ISEED3)*SigmaI*Ierr
       I(Year) = V(Year)+Idev
33000 CONTINUE
C 
C     Add measurement error to the 1+ biomass (i.e. an "assessment")
      EpsB(0) = 1
      DO 34000 Year = 1,NYear   
       EpsB(Year) = XNORM(1,0.0d0,1.0d0,ISEED4)*SigmaB - 
     +                  (SigmaB**2)/2.0
       EpsB(Year) = EXP(EpsB(Year)*Berr)
34000 CONTINUE	   
C     
C     Set the default for weight and selex values
      DO 40000 Year = 0,Nyear
       DO 41000 Iage = 0,Nage(1)
        SelxSum(Year,Iage) = Selex(Iage)
41000  CONTINUE
40000 CONTINUE
C
C     SENSITIVITIES
C     2. Time varying selectivities
      IF (TVselx.EQ.1) THEN
       DO 60000 Year=0,Nyear
C       Compute L
        Lselx = -31.407868 + 1.998082*V(Year)
        Lselx = MIN(1.0,Lselx)
        Lselx = MAX(0.0,Lselx)
C       Compute J
        Jselx = -23.790353 + 1.518223*V(Year)
        Jselx = MIN(1.0,Jselx)
        Jselx = MAX(0.0,Jselx)
C       Compute the MexCal selectivity for the year
        DO 61000 Iage = 0,Nage(1)
         SelxMCy(Iage) = Lselx*SelxMC1(Iage) + 
     +                  (1-Lselx)*SelxMC2(Iage)
61000   CONTINUE
C       Normalize (set the max selx to 1)
        MaxSelx = MAXVAL(SelxMCy)
        DO 62000 Iage = 0,Nage(1)
         SelxMCy(Iage) = SelxMCy(Iage)/MaxSelx
62000   CONTINUE
C       Now, compute the total selectivity
        DO 63000 Iage = 0,Nage(1)
         Selex(Iage) = Jselx*SelxMCy(Iage) +
     +                (1-Jselx)*SelxPNW(Iage)
63000   CONTINUE
C       And normalize again
        MaxSelx = MAXVAL(Selex)
        DO 64000 Iage = 0,Nage(1)
         Selex(Iage) = Selex(Iage)/MaxSelx
64000   CONTINUE
C
C       Save the selex values
        DO 65000 Iage = 0,Nage(1)
         SelxSum(Year,Iage) = Selex(Iage)
65000   CONTINUE
60000  CONTINUE
      ENDIF
C
C     SENSITIVITIES
C     4. Time varying natural mortality
      IF(TVM.EQ.1) THEN
       DO 70000 Year=0,Nyear
        Gslope = G(Year)-G(Year-5)
        IF(Year.LT.5) Gslope = 0.0
        Myear(Year) = 0.4
        IF(Gslope.LT.0.0) Myear(Year)=0.8
70000  CONTINUE
      ENDIF
C
C     SENSITIVITIES
C     5. Different selectivities without time varying
      IF(DiffSelx.EQ.1) THEN
       DO 40001 Year = 0,Nyear
        DO 41001 Iage = 0,Nage(1)
         SelxSum(Year,Iage) = SelxMC1(Iage)
41001   CONTINUE
40001  CONTINUE
      ENDIF
      IF(DiffSelx.EQ.2) THEN
       DO 40002 Year = 0,Nyear
        DO 41002 Iage = 0,Nage(1)
         SelxSum(Year,Iage) = SelxPNW(Iage)
41002   CONTINUE
40002  CONTINUE
      ENDIF
C
C     For debugging only
C      DO 50000 Year = 1,Nyear
C      WRITE(*,600) Year,G(Year),V(Year),I(Year),Eps(Year)
C50000 CONTINUE       
C
      RETURN
600   FORMAT(1x,I4,1x,A1,1x,7(F7.3,1x))      
      END
C
C ===============================================================================================
C
      SUBROUTINE Proj1YrPrey(Isim,Iyr,ZeroHG)
C
C This subroutine updates the basic dynamics
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Global variables
      INTEGER Iyr,Isim
      REAL*8 HG, HarvGuide,Eff
      LOGICAL ZeroHG
C
C     Local variables
      INTEGER Age,II,Ispec,Istk,Mage
      REAL*8 FF,Fmin,Fmax,Cnum,Cpred
      REAL*8 SSBV,Bio1P,RecR,SSBflat
      REAL*8 ANcount1, ANcount2, Ncheck
      REAL*8 ACcount1, ACcount2
      REAL*8 qind
      REAL*8 Z(MaxSpec,MaxStk,0:MaxAge)
      REAL*8 RickerA, RickerB
C
C     Set Ricker parameters
      RickerA = -4.598
      RickerB = -3.457e-07
      IF(Rtype.EQ.4) THEN
       RickerA = -3.6978
       RickerB = -1.4680e-06	
      ENDIF
C
C     Do the "Assessment"
      Bhat = B1p(1,1,Isim,Iyr-1)*EpsB(Iyr-1)
C	  
C     SENSITIVITIES 
C     Hyperstability
      IF (HyperSt.EQ.1.0) THEN
       qind = 1
       qind = HPSa*1/SQRT(B1p(1,1,Isim,Iyr-1))
       qind = MAX(qind,1.0)
       Bhat = qind*B1p(1,1,Isim,Iyr-1)*EpsB(Iyr-1)
      ENDIF
C
C     SENSITIVITIES
C     Time varying M (for sardine)
      IF(TVM.EQ.1) NatM(1) = Myear(Iyr-1)
C
C     Calculate the harvest guideline
      HG = 0
      Cpred = 0
      Eff = 0
C
      HG = HarvGuide(Isim,Iyr-1)
C     Set the minimum catch of 2000t for bait
      HG = MAX(HG,2000.0)
C
C     SENSITIVITIES
C     3. Multiple fleets - Multiply HG by distribution
      IF(MultFleet.EQ.1.0) THEN
C       Don't forget to still set the minimum to the bait catch
        HG = MAX(HG,0.0)
C       Multiply by DISTRIBUTION (0.87)
        HG = HG * 0.87
      ENDIF
C
C     Make sure the catch does not exceed the biomass
      IF(HG.GE.0.95*B1p(1,1,Isim,Iyr-1)) HG = B1p(1,1,Isim,Iyr-1)*0.95
C
C     Set HG to zero (to calculate B0)
      IF (ZeroHG) HG = 0
C
C     Save the F values
      Eff = HG/B1p(1,1,Isim,Iyr-1)
      IF(OFLopt.EQ.1) THEN
       Eff = (Emsymax-Emsymin)/(Nsim-1)*(Isim-1) + Emsymin
      ENDIF
C
C     Default Z
      DO 10000 Ispec = 1,Nspec
       DO 10000 IStk = 1,Nstk(Ispec)
        DO 10000 Age = 0,Nage(Ispec)
         Z(Ispec,Istk,Age) = NatM(Ispec)
10000 CONTINUE         
C
C     Solve for F for species 1, stock 1
C     ==================================
C	  
C     Start calculation of mean age of the catch
      AgeC(Isim,Iyr) = 0
C     Solve for the catch
      IF (HG.LE.0) THEN
       FF = 0
      ELSE
       Fmin = 0
       Fmax = 2.996
C      Use a bisection algorithm to guarantee that the HG is taken
       DO 11000 II = 1,23
        FF = (Fmin+Fmax)/2.0
        CPred = 0
        ACcount1 = 0
        ACcount2 = 0
        DO 11100 Age = 0,Nage(1)
         Z(1,1,Age) = SelxSum(Iyr-1,Age)*FF+NatM(1)
C
C        SENSITIVITIES
C        3. Multiple fleet - Update dynamics
         IF (MultFleet.EQ.1) THEN
          Z(1,1,Age) = SelxSum(Iyr-1,Age)*FF + NatM(1) + 
     +             SelxMC2(Age)*Fmex + SelxPNW(age)*Fcan
         ENDIF
C
         Cnum = N(1,1,Iyr-1,Age)*SelxSum(Iyr-1,Age)*FF/Z(1,1,Age)*
     +          (1-exp(-Z(1,1,Age)))
         Cpred = Cpred  + Weight(1,Age)*Cnum
         ACcount1 = ACcount1 + Cnum*Age
         ACcount2 = ACcount2 + Cnum
11100   CONTINUE 
        IF (Cpred.GT.HG) THEN
         Fmax = FF
        ELSE
         Fmin = FF
        ENDIF
11000  CONTINUE        
      ENDIF 
      F(Isim,Iyr-1) = Eff
      IF(Cpred.LE.0.0001) THEN
        AgeC(Isim,Iyr-1) = 0.000
        F(Isim,Iyr-1) = 0.000
      ELSE  
       AgeC(Isim,Iyr-1) = ACcount1/ACcount2
      ENDIF
      Catch(Isim,Iyr-1) = Cpred
C
C     Update the dynamics
      DO 20000 Ispec = 1,Nspec
       Mage = Nage(Ispec)
       DO 20000 Istk = 1,Nstk(Ispec) 
        N(Ispec,Istk,Iyr,Mage) = 
     +    N(Ispec,Istk,Iyr-1,Mage)*exp(-Z(Ispec,Istk,Mage))
        N(Ispec,Istk,Iyr,Mage) = N(Ispec,Istk,Iyr,Mage) + 
     +    N(Ispec,Istk,Iyr-1,Mage-1)*exp(-Z(Ispec,Istk,Mage-1))
        DO 20000 Age = 1,Mage-1
         N(Ispec,Istk,Iyr,Age) = 
     +    N(Ispec,Istk,Iyr-1,Age-1)*exp(-Z(Ispec,Istk,Age-1))
20000 CONTINUE
C
C     Update the dynamics
      AgeN(Isim,Iyr) = 0
      ANcount1 = NAge(1)*N(1,1,Iyr,Nage(1))
      ANcount2 = N(1,1,Iyr,Nage(1))
      DO 21000 Age = 1,Nage(1)-1
       ANcount1 = ANcount1 + N(1,1,Iyr,Age)*Age
       ANcount2 = ANcount2 + N(1,1,Iyr,Age)
21000 CONTINUE
C
C     Calculate mean age of the population
      AgeN(Isim,Iyr) = ANcount1/ANcount2
      IF(ANcount2.LE.1.0) AgeN(Isim,Iyr) = 0
C	  
      Ncheck = 0
      DO 20010 Age=0,Nage(1)
       Ncheck = Ncheck + N(1,1,Iyr,Age)
20010 CONTINUE
C
      IF(Ncheck.LE.0.0001) AgeN(Isim,Iyr-1) = 0.000
C
      IF(Iyr.EQ.Nyear) THEN
       IF(Ncheck.LE.0.0001) AgeN(Isim,Iyr) = 0.000 
      ENDIF
C
C     Compute the SSB and 1+ biomass
      DO 23000 Ispec = 1,Nspec
       DO 23000 Istk = 1,Nstk(Ispec)
        Bio1P = 0
        DO 23100 Age = 1,Nage(Ispec)
         Bio1P = Bio1P + Weight(Ispec,Age)*N(Ispec,Istk,Iyr,Age)
23100   CONTINUE
C
        SSBV = 0
        DO 23200 Age = 2,Nage(Ispec)
         SSBV = SSBV + Weight(Ispec,Age)*N(Ispec,Istk,Iyr,Age)
23200   CONTINUE
C
        SSB(Ispec,Istk,Isim,Iyr) = SSBV
        B1p(Ispec,Istk,Isim,Iyr) = Bio1P
C        
23000 CONTINUE        
C
C     Get the recruitment ((EQUATION 2))
      DO 24000 Ispec = 1,Nspec
       DO 24000 Istk = 1,Nstk(Ispec)
C
C       Base case is log(R/S) model ,
        IF (Ispec.EQ.1.AND.Istk.EQ.1) THEN   
         SSBV = SSB(Ispec,Istk,Isim,Iyr)
         RecR = SSBV/1000 * EXP(AlogRPS + BlogRPS*SSBV/1000 + 
     +       PhilogRPS*V(Iyr))
C        the log(R) model, dome-shaped      
         IF(Rtype.EQ.1) THEN
          RecR = EXP(Alpha + Beta1*SSBV/1000 + Beta2*(SSBV/1000)**2 +
     +        Phi*V(Iyr))
         ENDIF
C        The log(R) model, flat-top
C        First, set the upper limit of the recruitment curve
         SSBflat = -Beta1/(2*Beta2)
         IF(Rtype.EQ.2) THEN
          IF((SSBV/1000).GT.SSBflat) THEN
           RecR = EXP(Alpha+Beta1*SSBflat+Beta2*(SSBflat)**2+
     +         Phi*V(Iyr))
          ENDIF
         ENDIF
C        A simple Ricker model with no temperature term
         IF(Rtype.EQ.3) THEN
          RecR = SSBV*EXP(RickerA + RickerB * SSBV)
         ENDIF	   
         IF(Rtype.EQ.4) THEN
          RecR = SSBV*EXP(RickerA + RickerB * SSBV)
         ENDIF	 
         Recr = Recr * 1000
        ENDIF 
C
C       Other species and stocks
        IF (Ispec.NE.1.OR.IStk.NE.1) THEN
         Recr = Ninit(Ispec,Istk,0)
        ENDIF
C
C       Store recruitment
        RecR = RecR * EXP(Eps(Ispec,Istk,Iyr)-(SigmaR(Ispec)**2)/2.0)
        N(Ispec,Istk,Iyr,0) = RecR
        RecRep(Ispec,Istk,Isim,Iyr) = RecR
24000 CONTINUE       
C      
      RETURN
600   FORMAT(I4,1x,100(F12.3,1x))      
      END
C
C ===============================================================================================
C
      SUBROUTINE InitialPred(Isim)
C
C This subroutine sets up variables which are common across a run
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Global variables
      INTEGER Isim
C
C     Local variables
      INTEGER Iage,Ispec,Nmax
      REAL*8 MatPop

C     Initial N-vector
      DO 10000 Ispec = 1,Npred
       NP(Ispec,0,0) = 1
       Nmax = PNage(Ispec)
       DO 10100 Iage = 1,Nmax
        NP(Ispec,0,Iage) = NP(Ispec,0,Iage-1)*SP(Ispec,Iage-1)
10100  CONTINUE        
       NP(Ispec,0,Nmax) = NP(Ispec,0,Nmax)/(1-SP(Ispec,Nmax))
10000 CONTINUE
C
C     Check population is in balance at K
      DO 20000 Ispec = 1,Npred
       MatPop = 0
       DO 21000 Iage = Pmat(ispec),PNage(Ispec)
        MatPop = MatPop + NP(Ispec,0,Iage)
21000  CONTINUE
       WRITE(*,*) Ispec,MatPop
       PredK(Ispec) = 0
       DO 22000 Iage = 1,PNage(Ispec)
        PredK(Ispec) = PredK(Ispec) + NP(Ispec,0,Iage)
22000  CONTINUE
       N1P(Isim,Ispec,0) = PredK(Ispec)
       MatP(Isim,Ispec,0) = MatPop
20000 CONTINUE       
C
      RETURN
      END
C
C ============================================================================
C
      SUBROUTINE Proj1YrPred(Isim,Iyr)
C
C This subroutine updates the basic dynamics
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
C     Global variables
      INTEGER Iyr,Isim
C      
C     Local variables
      INTEGER Nmax,Ispec,Age,Jspec,Jstk
      REAL*8 Depl,DensDep,PredEn,RelPreyDep
      REAL*8 TOP,BOT1,BOT2,Phi2(MaxPspec),Gamma
C
C     Compute prey base
      DO 40000 Ispec = 1,Npred
       PredEn = 0 
       DO 40100 Jspec = 1,Nspec
        DO 40100 Jstk = 1,Nstk(Jspec)
         PredEn = PredEN + 
     +    B1p(Jspec,Jstk,Isim,Iyr-1)*PredWeight(Ispec,Jspec,Jstk)
40100  CONTINUE
       RelPreyDep = PredEn/B0Pred(Ispec)-PX(Ispec)
       TOP = (PX(Ispec)+PZ(Ispec)-1)*Psteep(Ispec)*RelPreyDep
       BOT1 = (1-PX(Ispec))*PZ(Ispec)*(Psteep(Ispec)-1)
       BOT2 = PZ(Ispec)-Psteep(Ispec)*(1.0-PX(Ispec))
       Phi2(Ispec) = TOP/(BOT1+BOT2*RelPreyDep)
       IF (RelPreyDep.LT.0) Phi2(Ispec) = 0
40000 CONTINUE
C
C     Update the dynamics
      DO 10000 Ispec = 1,Npred
       Gamma = 1
       IF (PredImpact.EQ.1) Gamma=Phi2(Ispec)
       Nmax = PNage(Ispec)
       NP(Ispec,Iyr,Nmax) = (NP(Ispec,Iyr-1,Nmax)*SP(ISpec,Nmax) +
     +  NP(Ispec,Iyr-1,Nmax-1)*SP(ISpec,Nmax-1))*Gamma
       DO 11000 Age = 1,Nmax-1
         NP(Ispec,Iyr,Age) = 
     +    NP(Ispec,Iyr-1,Age-1)*SP(Ispec,Age-1)*Gamma
11000  CONTINUE
C       WRITE(*,*) (NP(Ispec,Iyr-1,Age),Age=0,Nmax)
10000 CONTINUE
C
C     Compute 1+ abundance
      DO 20000 Ispec = 1,Npred
       N1P(Isim,Ispec,Iyr) = 0
       MatP(Isim,Ispec,Iyr) = 0
       DO 21000 Age = Pmat(Ispec),PNage(Ispec)
        MatP(Isim,Ispec,Iyr)= MatP(Isim,Ispec,Iyr) + NP(Ispec,Iyr,Age)
21000  CONTINUE
       DO 22000 Age = 1,PNage(Ispec)
        N1P(Isim,Ispec,Iyr) = N1P(Isim,Ispec,Iyr) + NP(Ispec,Iyr,Age)
22000  CONTINUE
20000 CONTINUE 
C
C     Compute recruitment
      DO 30000 Ispec = 1,Npred
       Depl = N1P(Isim,Ispec,Iyr)/PredK(Ispec)
       IF (PredImpact.EQ.3) Depl = Depl*Phi2(Ispec)
       IF (Depl.LE.0) THEN
        DensDep = 1
       ELSE 
        DensDep = 1.0+(MaxSJ(Ispec)-1.0)*(1.0-Depl**ZZ)
       ENDIF 
       NP(Ispec,Iyr,0) = DensDep*MatP(Isim,Ispec,Iyr)
       IF (PredImpact.EQ.2) NP(Ispec,Iyr,0) = 
     +    NP(Ispec,Iyr,0)*Phi2(Ispec)
30000 CONTINUE       
C
      RETURN
      END
C
C ===============================================================================================
C
      SUBROUTINE Summary()
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C      
C     Local variables
      INTEGER Isim,Year,Iage,Ispec,Istk
C
      DO 10000 Year = 0,Nyear
       WRITE(14,600) Year,
     +  ((
     +    (NINT(B1p(Ispec,Istk,Isim,Year)),
     +     NINT(SSB(Ispec,Istk,Isim,Year)),
     +     NINT(RecRep(Ispec,Istk,Isim,Year)),
     +   Istk=1,Nstk(Ispec)),Ispec=1,Nspec),Catch(Isim,Year),
     +                     F(Isim,Year),ObsF(Isim,Year),Vrep(Isim,Year),
     +                     AgeN(Isim,Year),AgeC(Isim,Year),Isim=1,Nsim)
10000 CONTINUE     
C
C     Print predator stuff
      WRITE(14,800) ((Ispec,Isim,Ispec=1,Npred),Isim=1,Nsim)
      DO 20000 Year = 0,Nyear
       WRITE(14,801) Year,
     +  ((NINT(N1p(Isim,Ispec,Year)*1000.0),
     +    NINT(MatP(Isim,Ispec,Year)*1000.0),
     +    Ispec=1,Npred),Isim=1,Nsim)  
20000 CONTINUE
C
      IF(SelxOut.EQ.1.0) THEN
       DO 11000 Year = 0,Nyear
        WRITE(14,700) Year, (SelxSum(Year,Iage),Iage=0,MaxAge)
11000 CONTINUE      
      ENDIF
C
C     These are only for debugging. Comment them once the model works      
C      WRITE(14,*) F(1,2), F(1,3), F(1,4)
C      WRITE(*,*) F(1,2), F(1,3), F(1,4)
C      WRITE(14,*) Ninit
C      WRITE(14,*) N(1,0), N(1,1), N(1,2), N(1,3), N(1,4), N(1,5)
C      WRITE(14,*) SSB(1,0)
C      WRITE(14,*) V(1),V(2),V(3),V(4),V(5),V(6),V(7),V(8),V(9),V(10)
C
C      DO 11000 Year=0,Nyear
C       WRITE(14,*) Eps(Year)
C11000 CONTINUE
C	   
      RETURN
600   FORMAT(1x,I5,1x,100(15(I9,1x),1x,F12.2,1x,F6.4,1x,
     +       F6.4,1x,F7.4,1x,F6.4,1x,F6.4))
700   FORMAT(1x,I5,1x,20(F5.3,1x))    
800   FORMAT(1x," Year ",100(I6,";",I6,1x))
801   FORMAT(1x,I5,1x,100(4(I6,1x)))
      END
C
C ===============================================================================================
C
      SUBROUTINE CalcG()
C
C This subroutine sets the values for the underlying environmental signal
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
C
      INTEGER Year
C
      REAL*8 Term1,Term2,Term3,Term4
      REAL*8 TwoPi
      TwoPi = 6.283185
C
      DO 10000 Year=0,Nyear
C      Create G
C       Sine wave (c)
       Term1 = Sin(TwoPi*(Year-Ybar)/Period)
C 
C       Basic square wave (a1)
       IF (MOD(Year,NINT(Period)).LT.Period/2.0) THEN
        Term2 = 1
        ELSE
         Term2 = -1
       ENDIF
C
C       Square wave with shorter good periods (a2)
       IF (MOD(Year,NINT(Period)).LT.Period/3.0) THEN
        Term3 = 1
        ELSE
         Term3 = -1
       ENDIF
C	   
C       Square wave with slowly descending shift from god to bad (b)
       IF (MOD(Year,NINT(Period)).LT.Period/2.0) THEN
        Term4 = 1
        ELSE
         IF (MOD(Year,NINT(Period)).GE.Period/2.0+10) THEN
           Term4 = -1
         ELSE
           Term4 = -0.2*(MOD(Year,NINT(Period))-Period/2.0)+1	   
         ENDIF
       ENDIF
C
       IF (Gtype.EQ.1) G(Year) = Term1
       IF (Gtype.EQ.2) G(Year) = Term2
       IF (Gtype.EQ.3) G(Year) = Term3
       IF (Gtype.EQ.4) G(Year) = Term4
C	   
10000 CONTINUE
C
      RETURN
      END
C
C ===============================================================================================
C
      REAL*8 FUNCTION HarvGuide(Isim,Iyr)
C
      IMPLICIT NONE
      INCLUDE "SARDINE.INC"
      
      REAL*8 HGt, Emsy,IEmsy,Frac,Imovave,ABC
      INTEGER Isim, Iyr
C
C     Calculate moving average of SST
      IF(Iyr.eq.0) Imovave = I(Iyr)
      IF(Iyr.eq.1) Imovave = I(Iyr)
      IF(Iyr.eq.2) Imovave = (I(Iyr)+I(Iyr-1))/2
      IF(Iyr.ge.3) Imovave = (I(Iyr)+I(Iyr-1)+I(Iyr-2))/3
      IF(UseAnn.eq.1) Imovave = I(Iyr)
C
C     Calculate Emsy
C     By default (OFLopt=0), set Emsy=Emsymin
      Emsy = Emsymin 
      IF(OFLopt.EQ.1) THEN
       Emsy = (Emsymax-Emsymin)/(Nsim-1)*(Isim-1) + Emsymin
      ENDIF
C
C     Calculate the OFL
      OFL = Emsy*Bhat
      IF (OFLopt.EQ.2) THEN
       Emsy = - 18.46452 + 3.25209*Imovave - 0.19723*Imovave**2 + 
     +        0.0041863*Imovave**3
       IEmsy = Emsy
       Emsy = MAX(Emsy,0.0)
       Emsy = MIN(Emsy,OFLmax)
       OFL = Emsy*Bhat	  
      ENDIF	  
C	  
C     Override OFL?
      IF (NoOFL.EQ.1) OFL = 0.95*Bhat
C     
C     Calculate the ABC (p*=0.40)
      ABC = Pstar*OFL
C	  
C     Cutoff can be looped over 
      IF (CutOpt.EQ.1) Cutoff = 50000*(Isim-1)
C
C     Calculate the HG
C     ================
C     By default (HCRopt=1), set the HG = Emsy
      HGt = OFL
      ObsF(Isim,Iyr) = Emsy
C     It can also be set with a cutoff and maxcatch and pre-specified Fraction=FracMin
C
C     Constant F but constrained by Maxcatch
      IF (HCRopt.EQ.2) THEN
       Frac = FracMin
       HGt = Frac*(Bhat-Cutoff)
       HGt = MIN(HGt,Maxcatch)
       HGt = MIN(HGt,ABC)
       HGt = MAX(HGt,0.0)
       ObsF(Isim,Iyr) = Frac
      ENDIF
C
C     Or include a temperature term to determine Fraction
      IF (HCRopt.EQ.3) THEN
       Frac = Emsy
       Frac = MIN(Frac,FracMax)
       Frac = MAX(Frac,FracMin)
       HGt = Frac*(Bhat-Cutoff)
       HGt = MIN(HGt,Maxcatch)
       HGt = MIN(HGt,ABC)
       HGt = MAX(HGt,0.0)
       ObsF(Isim,Iyr) = Frac
      ENDIF
C	  
      HarvGuide = HGt
C
      RETURN
      END
C
C ===============================================================================================
C
      INCLUDE "COMMON.FOR"