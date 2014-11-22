      INTEGER MaxYr,MaxAge,MaxSim,MaxStk,MaxSpec
      INTEGER MaxPspec,MaxPage
      PARAMETER(MaxYr=20000,MaxAge=20,MaxSim=100)
      PARAMETER (MaxSpec=5,MaxStk=2)
      PARAMETER (MaxPspec=10,MaxPage=50)
C
      COMMON / General / Nsim,Nspec,Nstk(MaxSpec)
      INTEGER Nsim,Nspec,Nstk
C      Nsim                  - Number of simulations      
C      Nspec                 - Number of prey species
C      Nstk                  - Number of stocks of prey species j
      
      COMMON / BiolPar / Nage,Nyear,Weight,Matu,Selex,Ninit,NatM,
     +     R0,Steepness,Alpha,Beta1,Beta2,Rerr,Rtype,AlogRPS,
     +     BlogRPS,PhilogRPS 
      INTEGER Nyear,Nage(MaxSpec)
      REAL*8 Weight(MaxSpec,0:MaxAge),Matu(0:MaxAge)
      REAL*8 Selex(0:MaxAge)
      REAL*8 SelxPNW(0:MaxAge),Ninit(MaxSpec,MaxStk,0:MaxAge)
      REAL*8 NatM(MaxSpec),R0(MaxSpec),Steepness(MaxSpec)
      REAL*8 Alpha,Beta1,Beta2,Rerr,Rtype
      REAL*8 AlogRPS,BlogRPS,PhilogRPS	  
C      Nyear                 - Number of projection years
C      Nage(s)               - Number of age-classes year      
C      Weight(s,a)           - Weight-at-age      
C      Matu(a)               - Maturity
C      Selex(a)              - Selectivity
C      Ninit(s,a)            - Initial numbers-at-age
C      NatM(s)               - Natural mortality
C      R0(s)                 - Virgin recruitment
C      Steepness(s)          - Steepness
C      Alpha,Beta1,Beta2     - SR parameters
C      Rerr                  - Use Recruitment uncertainty?
C      Rtype                 - What rec model to use (0=dome shaped; 1=flat top)

      COMMON / Dynamics / N,B1P,SSB,F,ObsF,Catch,RecRep,AgeN,AgeC,B0
      REAL*8 N(MaxSpec,MaxStk,0:MaxYr,0:MaxAge)
      REAL*8 B1p(MaxSpec,MaxStk,MaxSim,0:MaxYr)
      REAL*8 SSB(MaxSpec,MaxStk,MaxSim,0:MaxYr)
      REAL*8 F(MaxSim,0:MaxYr),ObsF(MaxSim,0:MaxYr)
      REAL*8 Catch(MaxSim,0:MaxYr)
      REAL*8 RecRep(MaxSpec,MaxStk,MaxSim,0:MaxYr)
      REAL*8 AgeN(MaxSim,0:MaxYr),AgeC(MaxSim,0:MaxYr)
      REAL*8 B0(MaxSpec,MaxStk)
C      N(y,a)                - Numbers-at-age      
C      B1p(s,y)              - Biomass of 1+ animals for year y and simulation s
C      SSB(s,y)              - SSB for year y and simulation s
C      F(s,y)                - True fishing mortality
C      ObsF(s,y)             - Fishing mortality from HG
C      Catch(s,y)            - Catch by year
C      RecRep(s,y)           - Recruitment by year
C      AgeN(s,y)             - Mean Age of the population
C      AgeC(s,y)             - Mean age of the catch
C      B0(s,y)               - Average Biomass at F=0
	  
      COMMON / Assessment / Bhat,Emsymin,Emsymax,Cutoff,FracMin,FracMax,
     +       Maxcatch,OFL,Pstar,OFLmax,SigmaB,Berr,CutOpt,EpsB,OFLopt,
     +       HCRopt,ISEED4,NoOFL,ISEED4Init
      REAL*8 Bhat,Emsymin,Emsymax,Cutoff,FracMin,FracMax,Maxcatch,OFL
      REAL*8 OFLmax,Pstar,SigmaB,Berr,CutOpt
      REAL*8 EpsB(0:MaxYr)
      INTEGER HCRopt,OFLopt,ISEED4,NoOFL,ISEED4Init
C      Bhat                  - Estimated biomass
C      Emsy                  - Value of Emsy
C      Cuttoff               - Cutoff parameter of the HG
C      Frac                  - Fraction parameter of the HG
C      Maxcatch              - Maxcatch parameter of the HG
C      SigmaB                - Variance for biomass estimate
C      HCRopt                - Which HCR should be used? 
C      OFLopt                - Which OFL should be used? 

      COMMON / EnvVar / Gtype,ISEED1,ISEED2,ISEED3,G,V,Eps,I,Vdev,
     +       EpsInit,ProwR,SigmaR,Vinit,Ybar,Period,Chi,Power,ProwV,
     +       SigmaV,Phi,MeanSST,SigmaI,Ierr,Verr,UseG,Amplitude,Vrep,
     +       ISEED1Init,ISEED2Init,ISEED3Init
      INTEGER Gtype
      INTEGER ISEED1,ISEED2,ISEED3
      INTEGER ISEED1Init,ISEED2Init,ISEED3Init
      REAL*8 G(0:MaxYr),V(0:MaxYr),Eps(MaxSpec,MaxStk,0:MaxYr)
      REAL*8 I(0:MaxYr),Vdev(0:MaxYr)
      REAL*8 EpsInit(MaxSpec,MaxStk),ProwR(MaxSpec),SigmaR(MaxSpec)
      REAL*8 Vinit,Ybar,Period,Chi,Power,ProwV,SigmaV,Phi,MeanSST
      REAL*8 SigmaI,Ierr,Verr,UseG,Amplitude
      REAL*8 Vrep(MaxSim,0:MaxYr)
C      Gtype                 - Type of variables      
C      G(y)                  - Underlying trends in the Environment variables
C      V(y)                  - Environmental variable
C      I(y)                  - Observed environmental variable
C      Eps(y)                - Rec_dev
C      EpsInit               - last Eps
C      ProwR                 - Prow for rec_dev
C      SigmaR                - Sigma for rec_dev
C      Vinit                 - last V
C      Ybar                  - Centering parameter
C      Period                - Period of the signal
C      Chi                   - Parameter for Type 1 relationship
C      Power                 - Parameter for Type 2 relationship
C      SigmaI                - Sigma for measurement error on the environmental index
C      Verr                  - Use environment uncertainty?
C
      COMMON / Sensitivities /  TVselx,MultFleet,HyperSt,TVM,
     +       HPSa,Fmex,Fcan,SelxOut,UseAnn,DiffSelx,
     +       SelxSum,SelxMC1,SelxMC2,SelxMCy,SelxPNW,
     +       Weight90s, Myear
      REAL*8 TVselx,MultFleet,HyperSt,TVM,HPSa,Fmex,Fcan
      REAL*8 SelxOut,UseAnn,DiffSelx
      REAL*8 SelxSum(0:MaxYr,0:MaxAge)
      REAL*8 SelxMC1(0:MaxAge),SelxMC2(0:MaxAge),SelxMCy(0:MaxAge)
      REAL*8 Weight90s(MaxSpec,0:MaxAge)
      REAL*8 Myear(0:MaxYr)

      COMMON / Predators1 / Npred,PNage,Pmat,PFmax,SP,
     +  MaxSJ,ZZ,PredWeight,Psteep,PZ,PX,PredImpact
      INTEGER Npred,PNage(MaxPspec),Pmat(MaxPspec),PredImpact
      REAL*8 PFmax(MaxPspec),SP(MaxPspec,0:MaxPage)
      REAL*8 MaxSJ(MaxPspec),ZZ
      REAL*8 PredWeight(MaxPspec,MaxSpec,MaxStk)
      REAL*8 Psteep(MaxPspec),PX(MaxPspec),PZ(MaxPspec)
C      Npred                 - Number of predators
C      PNAge(s)              - Number of size-classes
C      Pmat(s)               - Age at maturity
C      Fmax(s)               - Maximum rate of increase
C      ZZ                    - Density-dependence compensation
C      PredWeight(s,t,y)     - Energy value of predator s for stock y of prey t
C      PSteep(s)             - Predation function "steepness"
C      PX(s)                 - X for predation function
C      PZ(s)                 - Z for predation function
C      PredImpact            - Functional form for predator impact

      COMMON / Predators2 / NP,PredK,N1P,MatP,B0Pred
      REAL*8 NP(MaxPspec,0:MaxYr,0:MaxPAge)
      REAL*8 PredK(MaxPspec)
      REAL*8 N1P(MaxSim,MaxPspec,0:MaxYr)
      REAL*8 MatP(MaxSim,MaxPspec,0:MaxYr)
      REAL*8 B0Pred(MaxPspec)
C      NP(s,y,a)             - Predator populations      
C      PredK(s)              - Carrying capacity
C      N1P(s,y)              - 1+ abundance
C      B0Pred(s)             - Unfished energy content


