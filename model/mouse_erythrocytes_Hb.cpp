[PROB]

Model for blood cell dynamics healthy mouse (blood cell + Hb)

original model adopted from Zheng et al., 2021. 
https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12638

Adjustments on blood cells: 
- ST-HSC amplification time: 1000 -> 15
- MPP amplification time: 1000 -> 280
- CMP amplification time: 16 -> 8
- BFU-E amplification time: 32 -> 32
- CFU-E residence time: 7 days -> 2 days
- RBC lifespan: 120 days -> 40 days
- RET maturation time: 3 days -> 2 days

Adjustment on Hb:
- VRET = VRBC = 0.05e-12
scaling based on Fukuda et al., 2019
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658564/
See more details in the readme page

[SET]
delta = 0.5

[CMT] 

// blood cells
LT0
LT1
ST01
ST02
ST11
ST12
MPP01
MPP02
MPP11
MPP12
CMP01
CMP02
CMP11
CMP12
BFUE01
BFUE02
BFUE11
BFUE12
CFUE01
CFUE02
CFUE11
CFUE12
RET0
RET1
RBC0
RBC1

// hemoglobins
alpha_RET0
beta_RET0
gamma_RET0
delta_RET0
alphabeta_RET0
alphagamma_RET0
alphadelta_RET0
HbA_RET0
HbF_RET0
HbA2_RET0
alpha_RBC0
beta_RBC0
gamma_RBC0
delta_RBC0
alphabeta_RBC0
alphagamma_RBC0
alphadelta_RBC0
HbA_RBC0
HbF_RBC0
HbA2_RBC0

alpha_RET1
beta_RET1
gamma_RET1
delta_RET1
alphabeta_RET1
alphagamma_RET1
alphadelta_RET1
HbA_RET1
HbF_RET1
HbA2_RET1
alpha_RBC1
beta_RBC1
gamma_RBC1
delta_RBC1
alphabeta_RBC1
alphagamma_RBC1
alphadelta_RBC1
HbA_RBC1
HbF_RBC1
HbA2_RBC1

[PARAM]

// total mean residence time (days);
// mean residence time in each subcompartment will be computed later;
tauLT = 100
tauST = 20
tauMPP = 2
tauCMP = 4
tauBFUE = 7
tauCFUE = 2
tauRET = 2
tauRBCendo = 40
tauRBCtrans = 40 

ssLT = 1275 // total capacity for LT-HSC

// amplification parameters (A.U.);
aST = 15
aMPP = 280
aCMP = 8
aBFUE = 32

// replication rate for LT-HSC
rLT = (1/(2.5*7));

bloodvolume = 2e-3  // unit in L


// hemoglobin related parameters
VRET = 0.05e-12 // reticulocyte volume; L 
VRBC = 0.05e-12 // erythrocyte volume; L
MWHh = 64500 // hemoglobin molecular weight, g/mol or unit in ng/nmol

// change MWHh, following suggestions in Pittman, 2011
// https://www.ncbi.nlm.nih.gov/books/NBK54103/
// MWHh = 64400

ksynalpha = 3.8e-7 // alpha globin synthesis rate; nmol/day/cell
ratiosynbeta = 0.5 // beta globin/ alpha globin synthesis rate ratio; unitless; ratio assumed based on gene copy number
ratiosyngamma = 0.03 // gamma globin/ alpha globin synthesis rate ratio
ratiosyndelta = 0.04 // delta globin/ alpha globin synthesis rate ratio
thalfmonomer = 0.25 // free monomer half life; days
kon = 1e-5*(60*60*24) // bimolecular binding on rate constant; unit nmol-1.day-1 (value adjusted for unit)
Kdalphabeta = 1e-3 // αβ dimer dissociation constant, nM; 
Kdalphagamma = 1e-5 // αγ dimer dissociation constant, nM
Kdalphadelta = 1e-2  // αδ dimer dissociation constant, nM
KdHbA = 100 // HbA (α₂β₂) tetramer dissociation constant, nM
KdHbF = 100 // HbF (α₂γ₂) tetramer dissociation constant, nM
KdHbA2 = 100 // HbA2 (α₂δ₂) tetramer dissociation rate, nM

HbA_saturation = 0.74
HbF_saturation = 0.88
HbA2_saturation = 0


[MAIN]

// LT-HSC dynamics
double kLT2ST = 1/tauLT;
double deltaLT = (rLT - kLT2ST)/ssLT; 

// ST-HSC dynamics
double tau1ST = (log2(aST) -1)/log2(aST) * tauST;
double tau2ST = tauST - tau1ST;
double krep1ST = (aST/2-1)/ tau1ST;
double krep2ST = 1/ tau2ST;
double k12ST = aST/(2*tau1ST);
double kST2MPP = 2/tau2ST;

// MPP dynamics
double tau1MPP = (log2(aMPP)-1)/log2(aMPP) * tauMPP;
double tau2MPP = tauMPP - tau1MPP;
double krep1MPP = (aMPP/2-1)/ tau1MPP;
double krep2MPP = 1/ tau2MPP;
double k12MPP = aMPP/(2*tau1MPP);
double kMPP2CMP = 2/tau2MPP;

// CMP dynamics
double tau1CMP = (log2(aCMP)-1)/log2(aCMP) * tauCMP;
double tau2CMP = tauCMP - tau1CMP;
double krep1CMP = (aCMP/2-1)/ tau1CMP;
double krep2CMP = 1/ tau2CMP;
double k12CMP = aCMP/(2*tau1CMP);
double kCMP2BFUE = 2/tau2CMP;

// BFU-E dynamics
double tau1BFUE = (log2(aBFUE)-1)/log2(aBFUE) * tauBFUE;
double tau2BFUE = tauBFUE - tau1BFUE;
double krep1BFUE = (aBFUE/2-1)/ tau1BFUE;
double krep2BFUE = 1/ tau2BFUE;
double k12BFUE = aBFUE/(2*tau1BFUE);
double kBFUE2CFUE = 2/tau2BFUE;


// RET and RBC dynamics
double kRET2RBC = 1/tauRET;
double kdeathRBCendo = 1/tauRBCendo;
double kdeathRBCtrans = 1/tauRBCtrans;


/// monomer synthesis rate
double ksynbeta = ratiosynbeta * ksynalpha;
double ksyngamma = ratiosyngamma * ksynalpha;
double ksyndelta = ratiosyndelta * ksynalpha;

// monomer degredation rate
double kdeg = log(2)/ thalfmonomer; 

// dimer dissociation rate
double koffalphagamma = Kdalphagamma * kon; 
double koffalphadelta = Kdalphadelta * kon; 
double koffalphabeta = Kdalphabeta * kon; 

// tetramer dissociation rate
double koffHbA = KdHbA * kon; 
double koffHbF = KdHbF * kon; 
double koffHbA2 = KdHbA2 * kon; 

// total RET & RBC volume 
double VRET0 = VRET * RET0 ; 
double VRBC0 = VRBC * RBC0 ;
double VRET1 = VRET * RET1 ; 
double VRBC1 = VRBC * RBC1 ;


[ODE]

// hemoglobin concentration

double HbA = (HbA_RET0 + HbA_RET1 + HbA_RBC0 + HbA_RBC1) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 
double HbF = (HbF_RET0 + HbF_RET1 + HbF_RBC0 + HbF_RBC1) * MWHh*1e-9/ (10*bloodvolume); // unit in g/dL  
double HbA2 = (HbA2_RET0 + HbA2_RET1 + HbA2_RBC0 + HbA2_RBC1) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 


// calculate O2 in the blood
double vO2 = ( HbA_saturation * HbA + HbF_saturation * HbF + HbA2_saturation * HbA2 ) * 1.34;

// calculate CFU-E amplification time
double aCFUE = 550 * exp( -0.23 * vO2 );

// CFU-E dynamics
double tau1CFUE = (log2(aCFUE)-1)/log2(aCFUE) * tauCFUE;
double tau2CFUE = tauCFUE - tau1CFUE;
double krep1CFUE = (aCFUE/2-1)/ tau1CFUE;
double krep2CFUE = 1/ tau2CFUE;
double k12CFUE = aCFUE/(2*tau1CFUE);
double kCFUE2RET = 2/tau2CFUE;


// LT-HSC compartment
dxdt_LT0 = rLT * LT0 - deltaLT * LT0 * (LT0 + LT1) - kLT2ST * LT0; // endogenous branch
dxdt_LT1 = rLT * LT1 - deltaLT * LT1 * (LT0 + LT1) - kLT2ST * LT1; // transplanted branch

// ST-HSC compartment
dxdt_ST01 = kLT2ST * LT0 + (krep1ST - k12ST) * ST01; // endogenous branch
dxdt_ST02 = k12ST * ST01 + (krep2ST - kST2MPP) * ST02; // endogenous branch

dxdt_ST11 = kLT2ST * LT1 + (krep1ST - k12ST) * ST11;  // transplanted branch 
dxdt_ST12 = k12ST * ST11 + (krep2ST - kST2MPP) * ST12; // transplanted branch 


// MPP compartment
dxdt_MPP01 = kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; // endogenous branch
dxdt_MPP02 = k12MPP * MPP01 + (krep2MPP - kMPP2CMP) * MPP02; // endogenous branch

dxdt_MPP11 = kST2MPP * ST12 + (krep1MPP - k12MPP) * MPP11; // transplanted branch 
dxdt_MPP12 = k12MPP * MPP11 + (krep2MPP - kMPP2CMP) * MPP12; // transplanted branch 


// CMP compartment
dxdt_CMP01 = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01;  // endogenous branch
dxdt_CMP02 = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE) * CMP02; // endogenous branch

dxdt_CMP11 = kMPP2CMP * MPP12 + (krep1CMP - k12CMP) * CMP11;  // transplanted branch 
dxdt_CMP12 = k12CMP * CMP11 + (krep2CMP - kCMP2BFUE) * CMP12; // transplanted branch 


// BFU-E compartment
dxdt_BFUE01 = kCMP2BFUE * CMP02 + (krep1BFUE - k12BFUE) * BFUE01;   // endogenous branch
dxdt_BFUE02 = k12BFUE * BFUE01 + (krep2BFUE - kBFUE2CFUE) * BFUE02; // endogenous branch

dxdt_BFUE11 = kCMP2BFUE * CMP12 + (krep1BFUE - k12BFUE) * BFUE11;   // transplanted branch 
dxdt_BFUE12 = k12BFUE * BFUE11 + (krep2BFUE - kBFUE2CFUE) * BFUE12; // transplanted branch 



// CFU-E compartment
dxdt_CFUE01 = kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; // endogenous branch
dxdt_CFUE02 = k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;  // endogenous branch

dxdt_CFUE11 = kBFUE2CFUE * BFUE12 + (krep1CFUE - k12CFUE) * CFUE11; // transplanted branch 
dxdt_CFUE12 = k12CFUE * CFUE11 + (krep2CFUE - kCFUE2RET) * CFUE12;  // transplanted branch 


// RET compartment
dxdt_RET0 = kCFUE2RET * CFUE02 - kRET2RBC * RET0;  // endogenous branch
dxdt_RET1 = kCFUE2RET * CFUE12 - kRET2RBC * RET1;  // transplanted branch 


// RBC
dxdt_RBC0 = kRET2RBC * RET0 - kdeathRBCendo * RBC0;  // endogenous branch
dxdt_RBC1 = kRET2RBC * RET1 - kdeathRBCtrans * RBC1; // transplanted branch

// monomer in endogenous reticulocyte; 
dxdt_alpha_RET0 = ksynalpha * RET0 - (kdeg + kRET2RBC)*alpha_RET0 - (kon*alpha_RET0*beta_RET0 + kon*alpha_RET0*gamma_RET0 + kon*alpha_RET0*delta_RET0) + (koffalphabeta*alphabeta_RET0 + koffalphagamma*alphagamma_RET0 + koffalphadelta*alphadelta_RET0)*VRET0;
dxdt_beta_RET0  = ksynbeta  * RET0 - (kdeg + kRET2RBC)*beta_RET0  - kon * alpha_RET0 * beta_RET0  + koffalphabeta *VRET0 * alphabeta_RET0  ; 
dxdt_gamma_RET0 = ksyngamma * RET0 - (kdeg + kRET2RBC)*gamma_RET0 - kon * alpha_RET0 * gamma_RET0 + koffalphagamma      *VRET0 * alphagamma_RET0 ;
dxdt_delta_RET0 = ksyndelta * RET0 - (kdeg + kRET2RBC)*delta_RET0 - kon * alpha_RET0 * delta_RET0 + koffalphadelta      *VRET0 * alphadelta_RET0 ; 


// dimers in endogenous reticulocyte; 
dxdt_alphabeta_RET0  = kon * alpha_RET0 * beta_RET0   - kRET2RBC * alphabeta_RET0  - koffalphabeta * alphabeta_RET0  * VRET0 - 2 * kon * alphabeta_RET0  * alphabeta_RET0  + 2 * koffHbA  *VRET0 * HbA_RET0 ;
dxdt_alphagamma_RET0 = kon * alpha_RET0 * gamma_RET0  - kRET2RBC * alphagamma_RET0 - koffalphagamma      * alphagamma_RET0 * VRET0 - 2 * kon * alphagamma_RET0 * alphagamma_RET0 + 2 * koffHbF  *VRET0 * HbF_RET0 ; 
dxdt_alphadelta_RET0 = kon * alpha_RET0 * delta_RET0  - kRET2RBC * alphadelta_RET0 - koffalphadelta      * alphadelta_RET0 * VRET0 - 2 * kon * alphadelta_RET0 * alphadelta_RET0 + 2 * koffHbA2 *VRET0 * HbA2_RET0; 

// tetramer in endogenous reticulocyte; 
dxdt_HbA_RET0  = kon * alphabeta_RET0  * alphabeta_RET0  - (koffHbA  *VRET0 + kRET2RBC) * HbA_RET0;
dxdt_HbF_RET0  = kon * alphagamma_RET0 * alphagamma_RET0 - (koffHbF  *VRET0 + kRET2RBC) * HbF_RET0;
dxdt_HbA2_RET0 = kon * alphadelta_RET0 * alphadelta_RET0 - (koffHbA2 *VRET0 + kRET2RBC) * HbA2_RET0;

// monomer in endogenous RBC; 
dxdt_alpha_RBC0 = kRET2RBC * alpha_RET0 - (kdeg + kdeathRBCendo) * alpha_RBC0 - (kon*alpha_RBC0*beta_RBC0 + kon*alpha_RBC0*gamma_RBC0 + kon*alpha_RBC0*delta_RBC0) + (koffalphabeta * alphabeta_RBC0 + koffalphagamma * alphagamma_RBC0 + koffalphadelta * alphadelta_RBC0) * VRBC0; 
dxdt_beta_RBC0 =  kRET2RBC * beta_RET0  - (kdeg + kdeathRBCendo) * beta_RBC0  - kon * alpha_RBC0 * beta_RBC0  + koffalphabeta *VRBC0 * alphabeta_RBC0  ;
dxdt_gamma_RBC0 = kRET2RBC * gamma_RET0 - (kdeg + kdeathRBCendo) * gamma_RBC0 - kon * alpha_RBC0 * gamma_RBC0 + koffalphagamma      *VRBC0 * alphagamma_RBC0 ;
dxdt_delta_RBC0 = kRET2RBC * delta_RET0 - (kdeg + kdeathRBCendo) * delta_RBC0 - kon * alpha_RBC0 * delta_RBC0 + koffalphadelta      *VRBC0 * alphadelta_RBC0 ;

// dimers in endogenous RBC
dxdt_alphabeta_RBC0  = kRET2RBC * alphabeta_RET0  - kdeathRBCendo * alphabeta_RBC0  + kon * alpha_RBC0 * beta_RBC0  - koffalphabeta *VRBC0 * alphabeta_RBC0  - 2*kon* alphabeta_RBC0 *alphabeta_RBC0  + 2*koffHbA *VRBC0 * HbA_RBC0;
dxdt_alphagamma_RBC0 = kRET2RBC * alphagamma_RET0 - kdeathRBCendo * alphagamma_RBC0 + kon * alpha_RBC0 * gamma_RBC0 - koffalphagamma      *VRBC0 * alphagamma_RBC0 - 2*kon* alphagamma_RBC0*alphagamma_RBC0 + 2*koffHbF *VRBC0 * HbF_RBC0;
dxdt_alphadelta_RBC0 = kRET2RBC * alphadelta_RET0 - kdeathRBCendo * alphadelta_RBC0 + kon * alpha_RBC0 * delta_RBC0 - koffalphadelta      *VRBC0 * alphadelta_RBC0 - 2*kon* alphadelta_RBC0*alphadelta_RBC0 + 2*koffHbA2*VRBC0 * HbA2_RBC0;

// tetramers in endogenous RBC
dxdt_HbA_RBC0  = kRET2RBC * HbA_RET0  - kdeathRBCendo * HbA_RBC0  + kon*alphabeta_RBC0 *alphabeta_RBC0  - koffHbA * VRBC0 * HbA_RBC0;
dxdt_HbF_RBC0  = kRET2RBC * HbF_RET0  - kdeathRBCendo * HbF_RBC0  + kon*alphagamma_RBC0*alphagamma_RBC0 - koffHbF * VRBC0 * HbF_RBC0; 
dxdt_HbA2_RBC0 = kRET2RBC * HbA2_RET0 - kdeathRBCendo * HbA2_RBC0 + kon*alphadelta_RBC0*alphadelta_RBC0 - koffHbA2* VRBC0 * HbA2_RBC0; 


// monomer in transplanted reticulocyte; 
dxdt_alpha_RET1 = ksynalpha * RET1 - (kdeg + kRET2RBC)*alpha_RET1 - (kon*alpha_RET1*beta_RET1 + kon*alpha_RET1*gamma_RET1 + kon*alpha_RET1*delta_RET1) + (koffalphabeta*alphabeta_RET1 + koffalphagamma*alphagamma_RET1 + koffalphadelta*alphadelta_RET1)*VRET1;
dxdt_beta_RET1  = ksynbeta  * RET1 - (kdeg + kRET2RBC)*beta_RET1  - kon * alpha_RET1 * beta_RET1  + koffalphabeta *VRET1 * alphabeta_RET1  ; 
dxdt_gamma_RET1 = ksyngamma * RET1 - (kdeg + kRET2RBC)*gamma_RET1 - kon * alpha_RET1 * gamma_RET1 + koffalphagamma      *VRET1 * alphagamma_RET1 ;
dxdt_delta_RET1 = ksyndelta * RET1 - (kdeg + kRET2RBC)*delta_RET1 - kon * alpha_RET1 * delta_RET1 + koffalphadelta      *VRET1 * alphadelta_RET1 ; 


// dimers in transplanted reticulocyte; 
dxdt_alphabeta_RET1  = kon * alpha_RET1 * beta_RET1   - kRET2RBC * alphabeta_RET1  - koffalphabeta       * alphabeta_RET1  * VRET1 - 2 * kon * alphabeta_RET1  * alphabeta_RET1  + 2 * koffHbA  *VRET1 * HbA_RET1 ;
dxdt_alphagamma_RET1 = kon * alpha_RET1 * gamma_RET1  - kRET2RBC * alphagamma_RET1 - koffalphagamma      * alphagamma_RET1 * VRET1 - 2 * kon * alphagamma_RET1 * alphagamma_RET1 + 2 * koffHbF  *VRET1 * HbF_RET1 ; 
dxdt_alphadelta_RET1 = kon * alpha_RET1 * delta_RET1  - kRET2RBC * alphadelta_RET1 - koffalphadelta      * alphadelta_RET1 * VRET1 - 2 * kon * alphadelta_RET1 * alphadelta_RET1 + 2 * koffHbA2 *VRET1 * HbA2_RET1; 

// tetramer in transplanted reticulocyte; 
dxdt_HbA_RET1  = kon * alphabeta_RET1  * alphabeta_RET1  - (koffHbA  *VRET1 + kRET2RBC) * HbA_RET1;
dxdt_HbF_RET1  = kon * alphagamma_RET1 * alphagamma_RET1 - (koffHbF  *VRET1 + kRET2RBC) * HbF_RET1;
dxdt_HbA2_RET1 = kon * alphadelta_RET1 * alphadelta_RET1 - (koffHbA2 *VRET1 + kRET2RBC) * HbA2_RET1;

// monomer in transplanted RBC; 
dxdt_alpha_RBC1 = kRET2RBC * alpha_RET1 - (kdeg + kdeathRBCtrans) * alpha_RBC1 - (kon*alpha_RBC1*beta_RBC1 + kon*alpha_RBC1*gamma_RBC1 + kon*alpha_RBC1*delta_RBC1) + (koffalphabeta * alphabeta_RBC1 + koffalphagamma * alphagamma_RBC1 + koffalphadelta * alphadelta_RBC1) * VRBC1; 
dxdt_beta_RBC1 =  kRET2RBC * beta_RET1  - (kdeg + kdeathRBCtrans) * beta_RBC1  - kon * alpha_RBC1 * beta_RBC1  + koffalphabeta *VRBC1 * alphabeta_RBC1  ;
dxdt_gamma_RBC1 = kRET2RBC * gamma_RET1 - (kdeg + kdeathRBCtrans) * gamma_RBC1 - kon * alpha_RBC1 * gamma_RBC1 + koffalphagamma      *VRBC1 * alphagamma_RBC1 ;
dxdt_delta_RBC1 = kRET2RBC * delta_RET1 - (kdeg + kdeathRBCtrans) * delta_RBC1 - kon * alpha_RBC1 * delta_RBC1 + koffalphadelta      *VRBC1 * alphadelta_RBC1 ;

// dimers in transplanted RBC
dxdt_alphabeta_RBC1  = kRET2RBC * alphabeta_RET1  - kdeathRBCtrans * alphabeta_RBC1  + kon * alpha_RBC1 * beta_RBC1  - koffalphabeta *VRBC1 * alphabeta_RBC1  - 2*kon* alphabeta_RBC1 *alphabeta_RBC1  + 2*koffHbA *VRBC1 * HbA_RBC1;
dxdt_alphagamma_RBC1 = kRET2RBC * alphagamma_RET1 - kdeathRBCendo * alphagamma_RBC1 + kon * alpha_RBC1 * gamma_RBC1 - koffalphagamma      *VRBC1 * alphagamma_RBC1 - 2*kon* alphagamma_RBC1*alphagamma_RBC1 + 2*koffHbF *VRBC1 * HbF_RBC1;
dxdt_alphadelta_RBC1 = kRET2RBC * alphadelta_RET1 - kdeathRBCtrans * alphadelta_RBC1 + kon * alpha_RBC1 * delta_RBC1 - koffalphadelta      *VRBC1 * alphadelta_RBC1 - 2*kon* alphadelta_RBC1*alphadelta_RBC1 + 2*koffHbA2*VRBC1 * HbA2_RBC1;

// tetramers in transplanted RBC
dxdt_HbA_RBC1  = kRET2RBC * HbA_RET1  - kdeathRBCtrans * HbA_RBC1  + kon*alphabeta_RBC1 *alphabeta_RBC1  - koffHbA * VRBC1 * HbA_RBC1;
dxdt_HbF_RBC1  = kRET2RBC * HbF_RET1  - kdeathRBCtrans * HbF_RBC1  + kon*alphagamma_RBC1*alphagamma_RBC1 - koffHbF * VRBC1 * HbF_RBC1; 
dxdt_HbA2_RBC1 = kRET2RBC * HbA2_RET1 - kdeathRBCtrans * HbA2_RBC1 + kon*alphadelta_RBC1*alphadelta_RBC1 - koffHbA2* VRBC1 * HbA2_RBC1; 

[TABLE]


// blood cell ratio
capture RBCconc = (RBC0 + RBC1)/ (bloodvolume * 1e6); // RBC per uL
capture RETconc = (RET0 + RET1)/ (bloodvolume * 1e6); // RET per uL

capture totalHb = HbA + HbF + HbA2; // unit: g/dL
capture HbinRBC = (HbA_RBC0 + HbA_RBC1 + HbF_RBC0 + HbF_RBC1 + HbA2_RBC0 + HbA2_RBC1) * MWHh*1e-9 / (VRBC0 + VRBC1); // unit: g/L

capture CFUEamp = aCFUE; 
