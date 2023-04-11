[PROB]

Implementation of the blood cell dynamics of Zheng et al., 2021. 
https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12638

[SET]
delta = 0.1

[CMT] 

LT0
LT1
ST01
ST02
ST11
ST12
ST21
ST22
MPP01
MPP02
MPP11
MPP12
MPP21
MPP22
MPP31
MPP32
CMP01
CMP02
CMP11
CMP12
CMP21
CMP22
CMP31
CMP32
CMP41
CMP42
BFUE01
BFUE02
BFUE11
BFUE12
BFUE21
BFUE22
BFUE31
BFUE32
BFUE41
BFUE42
CFUE01
CFUE02
CFUE11
CFUE12
CFUE21
CFUE22
CFUE31
CFUE32
CFUE41
CFUE42
RET0
RET1
RET2
RET3
RET4
RBC0
RBC1
RBC2
RBC3
RBC4
alpha_RET0
beta_RET0
gamma_RET0
delta_RET0
alphabeta_RET0
alphagamma_RET0
alphadelta_RET0
HbS_RET0
HbF_RET0
HbA2_RET0
alpha_RBC0
beta_RBC0
gamma_RBC0
delta_RBC0
alphabeta_RBC0
alphagamma_RBC0
alphadelta_RBC0
HbS_RBC0
HbF_RBC0
HbA2_RBC0
alpha_RET1
beta_RET1
betanew_RET1
gamma_RET1
delta_RET1
alphabeta_RET1
alphabetanew_RET1
alphagamma_RET1
alphadelta_RET1
HbA_RET1
HbS_RET1
HbF_RET1
HbA2_RET1
alpha_RBC1
beta_RBC1
betanew_RBC1
gamma_RBC1
delta_RBC1
alphabeta_RBC1
alphabetanew_RBC1
alphagamma_RBC1
alphadelta_RBC1
HbS_RBC1
HbA_RBC1
HbF_RBC1
HbA2_RBC1
alpha_RET2
beta_RET2
betanew_RET2
gamma_RET2
delta_RET2
alphabeta_RET2
alphabetanew_RET2
alphagamma_RET2
alphadelta_RET2
HbA_RET2
HbS_RET2
HbF_RET2
HbA2_RET2
alpha_RBC2
beta_RBC2
betanew_RBC2
gamma_RBC2
delta_RBC2
alphabeta_RBC2
alphabetanew_RBC2
alphagamma_RBC2
alphadelta_RBC2
HbS_RBC2
HbA_RBC2
HbF_RBC2
HbA2_RBC2
alpha_RET3
beta_RET3
betanew_RET3
gamma_RET3
delta_RET3
alphabeta_RET3
alphabetanew_RET3
alphagamma_RET3
alphadelta_RET3
HbA_RET3
HbS_RET3
HbF_RET3
HbA2_RET3
alpha_RBC3
beta_RBC3
betanew_RBC3
gamma_RBC3
delta_RBC3
alphabeta_RBC3
alphabetanew_RBC3
alphagamma_RBC3
alphadelta_RBC3
HbS_RBC3
HbA_RBC3
HbF_RBC3
HbA2_RBC3
alpha_RET4
beta_RET4
betanew_RET4
gamma_RET4
delta_RET4
alphabeta_RET4
alphabetanew_RET4
alphagamma_RET4
alphadelta_RET4
HbA_RET4
HbS_RET4
HbF_RET4
HbA2_RET4
alpha_RBC4
beta_RBC4
betanew_RBC4
gamma_RBC4
delta_RBC4
alphabeta_RBC4
alphabetanew_RBC4
alphagamma_RBC4
alphadelta_RBC4
HbS_RBC4
HbA_RBC4
HbF_RBC4
HbA2_RBC4

[PARAM]

// total mean residence time (days);
// mean residence time in each subcompartment will be computed later;
tauLT = 100
tauST = 20
tauMPP = 2
tauCMP = 4
tauBFUE = 7
tauCFUE = 7
tauRET = 3
// tauRBChealthy = 120  // this parameter will be changed from outside
tauRBCsickle = 12
tauRBCtrans = 90

ssLT = 1275 // total capacity for LT-HSC

// amplification parameters (A.U.);
aST = 1000
aMPP = 1000
aCMP = 16
aBFUE = 32
// aCFUE = 32 // see appendix for more information. This implementation has issue

// total infused CD34+ cells; 
// this parameter is introduced to calculate initial condition
totalCD34infused = 2.8 * 1e8 // dummy parameter

// replication rate for LT-HSC
rLT = (1/(2.5*7)); // this is an assumed rate. The rate is missing from the original paper; this number should be greater than 0.01
// this value is taken from Abkowitz et al., 2000
// https://pubmed.ncbi.nlm.nih.gov/11071634/

bloodvolume = 5  // unit in L

// hemoglobin related parameters
VRET = 0.09e-12 // reticulocyte volume; L 
VRBC = 0.09e-12 // erythrocyte volume; L
MWHh = 64500 // hemoglobin molecular weight, g/mol or unit in ng/nmol

// change MWHh, following suggestions in Pittman, 2011
// https://www.ncbi.nlm.nih.gov/books/NBK54103/
// MWHh = 64400

ksynalpha = 1.5e-6 // alpha globin synthesis rate; nmol/day/cell
ratiosynbeta = 0.5 // beta globin/ alpha globin synthesis rate ratio; unitless; ratio assumed based on gene copy number
ratiosynbetanew = 0.5
ratiosyngamma = 0.03 // gamma globin/ alpha globin synthesis rate ratio
ratiosyndelta = 0.04 // delta globin/ alpha globin synthesis rate ratio
thalfmonomer = 0.25 // free monomer half life; days
kon = 1e-5*(60*60*24) // bimolecular binding on rate constant; unit nmol-1.day-1 (value adjusted for unit)
Kdalphabeta_trans = 1e-3 // αβ dimer dissociation constant, nM; this is for either healthy or transgene
Kdalphabeta_sickle = 1e-2 // αβˢ dimer dissociation constant, nM; this is for sickle gene
Kdalphagamma = 1e-5 // αγ dimer dissociation constant, nM
Kdalphadelta = 1e-2  // αδ dimer dissociation constant, nM
KdHbS = 100 // HbS (α₂β₂ˢ) tetramer dissociation constant, nM
KdHbA = 100 // HbA (α₂β₂) tetramer dissociation constant, nM
KdHbF = 100 // HbF (α₂γ₂) tetramer dissociation constant, nM
KdHbA2 = 100 // HbA2 (α₂δ₂) tetramer dissociation rate, nM


// these parameters are added for further optimization
HbA_saturation = 0.74
HbS_saturation = 0.68
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
double tauRBCendo = tauRBCsickle;
double kdeathRBCendo = 1/tauRBCendo;
double kdeathRBCtrans = 1/tauRBCtrans;

// transduced cells
LT1_0 = 1e-6 * totalCD34infused; 
ST22_0 = 1e-5 * totalCD34infused; 
MPP32_0 = 1e-3 * totalCD34infused;
CMP42_0 = 1e-2 * totalCD34infused; 

// monomer synthesis rate
double ksynbeta = ratiosynbeta * ksynalpha;
double ksyngamma = ratiosyngamma * ksynalpha;
double ksyndelta = ratiosyndelta * ksynalpha;

// monomer degredation rate
double kdeg = log(2)/ thalfmonomer; // all the monomer half life time is the same, so only 1 parameter is used

// dimer dissociation rate
double koffalphabetatrans = Kdalphabeta_trans * kon;
double koffalphagamma = Kdalphagamma * kon; 
double koffalphadelta = Kdalphadelta * kon; 
double koffalphabetasickle = Kdalphabeta_sickle * kon; 

// tetramer dissociation rate
double koffHbA = KdHbA * kon; 
double koffHbS = KdHbS * kon; 
double koffHbF = KdHbF * kon; 
double koffHbA2 = KdHbA2 * kon; 

// total RET & RBC volume 
// remove the plus one indicated in the model to avoid problems
double VRET0 = VRET * RET0 ; 
double VRET1 = VRET * RET1 ; 
double VRET2 = VRET * RET2 ; 
double VRET3 = VRET * RET3 ; 
double VRET4 = VRET * RET4 ; 

double VRBC0 = VRBC * RBC0 ;
double VRBC1 = VRBC * RBC1 ;
double VRBC2 = VRBC * RBC2 ;
double VRBC3 = VRBC * RBC3 ;
double VRBC4 = VRBC * RBC4 ;

// set INIT values
// endogenous cells; all these number are guessing numbers
// LT0_0 = 1; 



[ODE]

// calculate hemoglobin concentration
double HbA = (HbA_RET1 + HbA_RET2 + HbA_RET3 + HbA_RET4 + HbA_RBC1 + HbA_RBC2 + HbA_RBC3 + HbA_RBC4) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 
double HbS = (HbS_RET0 + HbS_RET1 + HbS_RET2 + HbS_RET3 + HbS_RET4 + HbS_RBC0 + HbS_RBC1 + HbS_RBC2 + HbS_RBC3 + HbS_RBC4) * MWHh*1e-9/ (10*bloodvolume); // unit in g/dL  
double HbF = (HbF_RET0 + HbF_RET1 + HbF_RET2 + HbF_RET3 + HbF_RET4 + HbF_RBC0 + HbF_RBC1 + HbF_RBC2 + HbF_RBC3 + HbF_RBC4) * MWHh*1e-9/ (10*bloodvolume); // unit in g/dL  
double HbA2 = (HbA2_RET0 + HbA2_RET1 + HbA2_RET2 + HbA2_RET3 + HbA2_RET4 + HbA2_RBC0 + HbA2_RBC1 + HbA2_RBC2 + HbA2_RBC3 + HbA2_RBC4) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 


// calculate O2 in the blood
// double vO2 = ( 0.74 * HbA + 0.68 * HbS + 0.88 * HbF ) * 1.34; 
// parameterize how this is calculated so that this can be optimized later
// HbA_saturation = 0.74
// HbS_saturation = 0.68
// HbF_saturation = 0.88
// HbA2_saturation = 0
double vO2 = ( HbA_saturation * HbA + HbS_saturation * HbS + HbF_saturation * HbF + HbA2_saturation * HbA2 ) * 1.34;

// assuming HbA2 also carries O2, and the carrying capacity is similar to HbF
// assumption taken from Manning et al., 2020. https://link.springer.com/chapter/10.1007/978-3-030-41769-7_11
// double vO2 = ( 0.74 * HbA + 0.68 * HbS + 0.88 * HbF + 0.88 * HbA2 ) * 1.34;
// try another way, assuming all the Hb saturation in uniform
// double vO2 = (HbA + HbS + HbF + HbA2) * 0.25 * 1.34;  after this adjustment, the blood cell number works out, but the hemoglobin concentration is too high
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
dxdt_LT1 = rLT * LT1 - deltaLT * LT1 * (LT0 + LT1) - kLT2ST * LT1; // transduced branch

// ST-HSC compartment
dxdt_ST01 = kLT2ST * LT0 + (krep1ST - k12ST) * ST01; // endogenous branch
dxdt_ST02 = k12ST * ST01 + (krep2ST - kST2MPP) * ST02; // endogenous branch

dxdt_ST11 = kLT2ST * LT1 + (krep1ST - k12ST) * ST11;  // transduced branch from infused LT-HSC
dxdt_ST12 = k12ST * ST11 + (krep2ST - kST2MPP) * ST12; // transduced branch from infused LT-HSC

dxdt_ST21 = (krep1ST - k12ST) * ST21; // transduced branch from infused ST-HSC
dxdt_ST22 = k12ST * ST21 + (krep2ST - kST2MPP) * ST22; // transduced branch from infused ST-HSC

// MPP compartment
dxdt_MPP01 = kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; // endogenous branch
dxdt_MPP02 = k12MPP * MPP01 + (krep2MPP - kMPP2CMP) * MPP02; // endogenous branch

dxdt_MPP11 = kST2MPP * ST12 + (krep1MPP - k12MPP) * MPP11; // transduced branch from infused LT-HSC
dxdt_MPP12 = k12MPP * MPP11 + (krep2MPP - kMPP2CMP) * MPP12; // transduced branch from infused LT-HSC

dxdt_MPP21 = kST2MPP * ST22 + (krep1MPP - k12MPP) * MPP21;  // transduced branch from infused ST-HSC
dxdt_MPP22 = k12MPP * MPP21 + (krep2MPP - kMPP2CMP) * MPP22; // transduced branch from infused ST-HSC

dxdt_MPP31 = (krep1MPP - k12MPP) * MPP31; // transduced branch from infused MPP
dxdt_MPP32 = k12MPP * MPP31 + (krep2MPP - kMPP2CMP) * MPP32; // transduced branch from infused MPP

// CMP compartment
dxdt_CMP01 = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01;  // endogenous branch
dxdt_CMP02 = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE) * CMP02; // endogenous branch

dxdt_CMP11 = kMPP2CMP * MPP12 + (krep1CMP - k12CMP) * CMP11;  // transduced branch from infused LT-HSC
dxdt_CMP12 = k12CMP * CMP11 + (krep2CMP - kCMP2BFUE) * CMP12; // transduced branch from infused LT-HSC

dxdt_CMP21 = kMPP2CMP * MPP22 + (krep1CMP - k12CMP) * CMP21;  // transduced branch from infused ST-HSC
dxdt_CMP22 = k12CMP * CMP21 + (krep2CMP - kCMP2BFUE) * CMP22; // transduced branch from infused ST-HSC

dxdt_CMP31 = kMPP2CMP * MPP32 + (krep1CMP - k12CMP) * CMP31;  // transduced branch from infused MPP
dxdt_CMP32 = k12CMP * CMP31 + (krep2CMP - kCMP2BFUE) * CMP32; // transduced branch from infused MPP

dxdt_CMP41 = (krep1CMP - k12CMP) * CMP41; // transduced branch from infused CMP
dxdt_CMP42 = k12CMP * CMP41 + (krep2CMP - kCMP2BFUE) * CMP42; // transduced branch from infused CMP

// BFU-E compartment
dxdt_BFUE01 = kCMP2BFUE * CMP02 + (krep1BFUE - k12BFUE) * BFUE01;   // endogenous branch
dxdt_BFUE02 = k12BFUE * BFUE01 + (krep2BFUE - kBFUE2CFUE) * BFUE02; // endogenous branch

dxdt_BFUE11 = kCMP2BFUE * CMP12 + (krep1BFUE - k12BFUE) * BFUE11;   // transduced branch from infused LT-HSC
dxdt_BFUE12 = k12BFUE * BFUE11 + (krep2BFUE - kBFUE2CFUE) * BFUE12; // transduced branch from infused LT-HSC

dxdt_BFUE21 = kCMP2BFUE * CMP22 + (krep1BFUE - k12BFUE) * BFUE21;   // transduced branch from infused ST-HSC
dxdt_BFUE22 = k12BFUE * BFUE21 + (krep2BFUE - kBFUE2CFUE) * BFUE22; // transduced branch from infused ST-HSC

dxdt_BFUE31 = kCMP2BFUE * CMP32 + (krep1BFUE - k12BFUE) * BFUE31;   // transduced branch from infused MPP
dxdt_BFUE32 = k12BFUE * BFUE31 + (krep2BFUE - kBFUE2CFUE) * BFUE32; // transduced branch from infused MPP 

dxdt_BFUE41 = kCMP2BFUE * CMP42 + (krep1BFUE - k12BFUE) * BFUE41;   // transduced branch from infused CMP 
dxdt_BFUE42 = k12BFUE * BFUE41 + (krep2BFUE - kBFUE2CFUE) * BFUE42; // transduced branch from infused CMP

// CFU-E compartment
dxdt_CFUE01 = kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; // endogenous branch
dxdt_CFUE02 = k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;  // endogenous branch

// CFU-E compartment
dxdt_CFUE01 = kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; // endogenous branch
dxdt_CFUE02 = k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;  // endogenous branch

dxdt_CFUE11 = kBFUE2CFUE * BFUE12 + (krep1CFUE - k12CFUE) * CFUE11; // transduced branch from infused LT-HSC
dxdt_CFUE12 = k12CFUE * CFUE11 + (krep2CFUE - kCFUE2RET) * CFUE12;  // transduced branch from infused LT-HSC

dxdt_CFUE21 = kBFUE2CFUE * BFUE22 + (krep1CFUE - k12CFUE) * CFUE21; // transduced branch from infused ST-HSC
dxdt_CFUE22 = k12CFUE * CFUE21 + (krep2CFUE - kCFUE2RET) * CFUE22;  // transduced branch from infused ST-HSC

dxdt_CFUE31 = kBFUE2CFUE * BFUE32 + (krep1CFUE - k12CFUE) * CFUE31; // transduced branch from infused MPP
dxdt_CFUE32 = k12CFUE * CFUE31 + (krep2CFUE - kCFUE2RET) * CFUE32;  // transduced branch from infused MPP

dxdt_CFUE41 = kBFUE2CFUE * BFUE42 + (krep1CFUE - k12CFUE) * CFUE41; // transduced branch from infused CMP 
dxdt_CFUE42 = k12CFUE * CFUE41 + (krep2CFUE - kCFUE2RET) * CFUE42;  // transduced branch from infused CMP 

// RET compartment
dxdt_RET0 = kCFUE2RET * CFUE02 - kRET2RBC * RET0;  // endogenous branch
dxdt_RET1 = kCFUE2RET * CFUE12 - kRET2RBC * RET1;  // transduced branch from infused LT-HSC
dxdt_RET2 = kCFUE2RET * CFUE22 - kRET2RBC * RET2;  // transduced branch from infused ST-HSC
dxdt_RET3 = kCFUE2RET * CFUE32 - kRET2RBC * RET3;  // transduced branch from infused MPP
dxdt_RET4 = kCFUE2RET * CFUE42 - kRET2RBC * RET4;  // transduced branch from infused CMP 

// RBC
dxdt_RBC0 = kRET2RBC * RET0 - kdeathRBCendo * RBC0;  // endogenous branch
dxdt_RBC1 = kRET2RBC * RET1 - kdeathRBCtrans * RBC1; // transduced branch
dxdt_RBC2 = kRET2RBC * RET2 - kdeathRBCtrans * RBC2; // transduced branch
dxdt_RBC3 = kRET2RBC * RET3 - kdeathRBCtrans * RBC3; // transduced branch
dxdt_RBC4 = kRET2RBC * RET4 - kdeathRBCtrans * RBC4; // transduced branch

// monomer in endogenous reticulocyte; 
dxdt_alpha_RET0 = ksynalpha * RET0 - (kdeg + kRET2RBC)*alpha_RET0 - (kon*alpha_RET0*beta_RET0 + kon*alpha_RET0*gamma_RET0 + kon*alpha_RET0*delta_RET0) + (koffalphabetasickle*alphabeta_RET0 + koffalphagamma*alphagamma_RET0 + koffalphadelta*alphadelta_RET0)*VRET0;
dxdt_beta_RET0  = ksynbeta  * RET0 - (kdeg + kRET2RBC)*beta_RET0  - kon * alpha_RET0 * beta_RET0  + koffalphabetasickle *VRET0 * alphabeta_RET0  ; 
dxdt_gamma_RET0 = ksyngamma * RET0 - (kdeg + kRET2RBC)*gamma_RET0 - kon * alpha_RET0 * gamma_RET0 + koffalphagamma      *VRET0 * alphagamma_RET0 ;
dxdt_delta_RET0 = ksyndelta * RET0 - (kdeg + kRET2RBC)*delta_RET0 - kon * alpha_RET0 * delta_RET0 + koffalphadelta      *VRET0 * alphadelta_RET0 ; 


// dimers in endogenous reticulocyte; 
dxdt_alphabeta_RET0  = kon * alpha_RET0 * beta_RET0 -  koffalphabetasickle * alphabeta_RET0  * VRET0 - 2 * kon * alphabeta_RET0  * alphabeta_RET0  + 2 * koffHbS  *VRET0 * HbS_RET0  - kRET2RBC * alphabeta_RET0;
dxdt_alphagamma_RET0 = kon * alpha_RET0 * gamma_RET0 - koffalphagamma      * alphagamma_RET0 * VRET0 - 2 * kon * alphagamma_RET0 * alphagamma_RET0 + 2 * koffHbF  *VRET0 * HbF_RET0  - kRET2RBC * alphagamma_RET0; 
dxdt_alphadelta_RET0 = kon * alpha_RET0 * delta_RET0 - koffalphadelta      * alphadelta_RET0 * VRET0 - 2 * kon * alphadelta_RET0 * alphadelta_RET0 + 2 * koffHbA2 *VRET0 * HbA2_RET0 - kRET2RBC * alphadelta_RET0; 

// tetramer in endogenous reticulocyte; 
dxdt_HbS_RET0  = kon * alphabeta_RET0  * alphabeta_RET0  - (koffHbS  *VRET0 + kRET2RBC) * HbS_RET0;
dxdt_HbF_RET0  = kon * alphagamma_RET0 * alphagamma_RET0 - (koffHbF  *VRET0 + kRET2RBC) * HbF_RET0;
dxdt_HbA2_RET0 = kon * alphadelta_RET0 * alphadelta_RET0 - (koffHbA2 *VRET0 + kRET2RBC) * HbA2_RET0;

// monomer in endogenous RBC; 
dxdt_alpha_RBC0 = kRET2RBC * alpha_RET0 - (kdeg + kdeathRBCendo) * alpha_RBC0 - (kon*alpha_RBC0*beta_RBC0 + kon*alpha_RBC0*gamma_RBC0 + kon*alpha_RBC0*delta_RBC0) + (koffalphabetasickle * alphabeta_RBC0 + koffalphagamma * alphagamma_RBC0 + koffalphadelta * alphadelta_RBC0) * VRBC0; 
dxdt_beta_RBC0 =  kRET2RBC * beta_RET0  - (kdeg + kdeathRBCendo) * beta_RBC0  - kon * alpha_RBC0 * beta_RBC0  + koffalphabetasickle *VRBC0 * alphabeta_RBC0  ;
dxdt_gamma_RBC0 = kRET2RBC * gamma_RET0 - (kdeg + kdeathRBCendo) * gamma_RBC0 - kon * alpha_RBC0 * gamma_RBC0 + koffalphagamma      *VRBC0 * alphagamma_RBC0 ;
dxdt_delta_RBC0 = kRET2RBC * delta_RET0 - (kdeg + kdeathRBCendo) * delta_RBC0 - kon * alpha_RBC0 * delta_RBC0 + koffalphadelta      *VRBC0 * alphadelta_RBC0 ;

// dimers in endogenous RBC
dxdt_alphabeta_RBC0  = kRET2RBC * alphabeta_RET0  - kdeathRBCendo * alphabeta_RBC0  + kon * alpha_RBC0 * beta_RBC0  - koffalphabetasickle *VRBC0 * alphabeta_RBC0  - 2*kon* alphabeta_RBC0 *alphabeta_RBC0  + 2*koffHbS *VRBC0 * HbS_RBC0;
dxdt_alphagamma_RBC0 = kRET2RBC * alphagamma_RET0 - kdeathRBCendo * alphagamma_RBC0 + kon * alpha_RBC0 * gamma_RBC0 - koffalphagamma      *VRBC0 * alphagamma_RBC0 - 2*kon* alphagamma_RBC0*alphagamma_RBC0 + 2*koffHbF *VRBC0 * HbF_RBC0;
dxdt_alphadelta_RBC0 = kRET2RBC * alphadelta_RET0 - kdeathRBCendo * alphadelta_RBC0 + kon * alpha_RBC0 * delta_RBC0 - koffalphadelta      *VRBC0 * alphadelta_RBC0 - 2*kon* alphadelta_RBC0*alphadelta_RBC0 + 2*koffHbA2*VRBC0 * HbA2_RBC0;

// tetramers in endogenous RBC
dxdt_HbS_RBC0  = kRET2RBC * HbS_RET0  - kdeathRBCendo * HbS_RBC0  + kon*alphabeta_RBC0 *alphabeta_RBC0  - koffHbS * VRBC0 * HbS_RBC0;
dxdt_HbF_RBC0  = kRET2RBC * HbF_RET0  - kdeathRBCendo * HbF_RBC0  + kon*alphagamma_RBC0*alphagamma_RBC0 - koffHbF * VRBC0 * HbF_RBC0; 
dxdt_HbA2_RBC0 = kRET2RBC * HbA2_RET0 - kdeathRBCendo * HbA2_RBC0 + kon*alphadelta_RBC0*alphadelta_RBC0 - koffHbA2* VRBC0 * HbA2_RBC0; 

// monomer in transduced RET
dxdt_alpha_RET1   = ksynalpha                      * RET1 - (kdeg + kRET2RBC)* alpha_RET1   - (kon*alpha_RET1*beta_RET1 + kon*alpha_RET1*betanew_RET1 + kon*alpha_RET1*gamma_RET1 + kon*alpha_RET1*delta_RET1) + (koffalphabetasickle*alphabeta_RET1 + koffalphabetatrans*alphabetanew_RET1 + koffalphagamma*alphagamma_RET1 + koffalphadelta*alphadelta_RET1)*VRET1;
dxdt_beta_RET1    = ksynbeta * (1-ratiosynbetanew) * RET1 - (kdeg + kRET2RBC)* beta_RET1    - kon*alpha_RET1*beta_RET1    + koffalphabetasickle *VRET1 * alphabeta_RET1; 
dxdt_betanew_RET1 = ksynbeta * ratiosynbetanew     * RET1 - (kdeg + kRET2RBC)* betanew_RET1 - kon*alpha_RET1*betanew_RET1 + koffalphabetatrans  *VRET1 * alphabetanew_RET1;
dxdt_gamma_RET1   = ksyngamma                      * RET1 - (kdeg + kRET2RBC)* gamma_RET1   - kon*alpha_RET1*gamma_RET1   + koffalphagamma      *VRET1 * alphagamma_RET1;
dxdt_delta_RET1   = ksyndelta                      * RET1 - (kdeg + kRET2RBC)* delta_RET1   - kon*alpha_RET1*delta_RET1   + koffalphadelta      *VRET1 * alphadelta_RET1;

dxdt_alpha_RET2   = ksynalpha                      * RET2 - (kdeg + kRET2RBC)* alpha_RET2   - (kon*alpha_RET2*beta_RET2 + kon*alpha_RET2*betanew_RET2 + kon*alpha_RET2*gamma_RET2 + kon*alpha_RET2*delta_RET2) + (koffalphabetasickle*alphabeta_RET2 + koffalphabetatrans*alphabetanew_RET2 + koffalphagamma*alphagamma_RET2 + koffalphadelta*alphadelta_RET2)*VRET2;
dxdt_beta_RET2    = ksynbeta * (1-ratiosynbetanew) * RET2 - (kdeg + kRET2RBC)* beta_RET2    - kon*alpha_RET2*beta_RET2    + koffalphabetasickle *VRET2 * alphabeta_RET2 ; 
dxdt_betanew_RET2 = ksynbeta * ratiosynbetanew     * RET2 - (kdeg + kRET2RBC)* betanew_RET2 - kon*alpha_RET2*betanew_RET2 + koffalphabetatrans  *VRET2 * alphabetanew_RET2 ;
dxdt_gamma_RET2   = ksyngamma                      * RET2 - (kdeg + kRET2RBC)* gamma_RET2   - kon*alpha_RET2*gamma_RET2   + koffalphagamma      *VRET2 * alphagamma_RET2 ;
dxdt_delta_RET2   = ksyndelta                      * RET2 - (kdeg + kRET2RBC)* delta_RET2   - kon*alpha_RET2*delta_RET2   + koffalphadelta      *VRET2 * alphadelta_RET2 ;

dxdt_alpha_RET3   = ksynalpha                      * RET3 - (kdeg + kRET2RBC)* alpha_RET3   - (kon*alpha_RET3*beta_RET3 + kon*alpha_RET3*betanew_RET3 + kon*alpha_RET3*gamma_RET3 + kon*alpha_RET3*delta_RET3) + (koffalphabetasickle*alphabeta_RET3 + koffalphabetatrans*alphabetanew_RET3 + koffalphagamma*alphagamma_RET3 + koffalphadelta*alphadelta_RET3)*VRET3;
dxdt_beta_RET3    = ksynbeta * (1-ratiosynbetanew) * RET3 - (kdeg + kRET2RBC)* beta_RET3    - kon*alpha_RET3*beta_RET3    + koffalphabetasickle *VRET3* alphabeta_RET3 ; 
dxdt_betanew_RET3 = ksynbeta * ratiosynbetanew     * RET3 - (kdeg + kRET2RBC)* betanew_RET3 - kon*alpha_RET3*betanew_RET3 + koffalphabetatrans  *VRET3* alphabetanew_RET3;
dxdt_gamma_RET3   = ksyngamma                      * RET3 - (kdeg + kRET2RBC)* gamma_RET3   - kon*alpha_RET3*gamma_RET3   + koffalphagamma      *VRET3* alphagamma_RET3;
dxdt_delta_RET3  =  ksyndelta                      * RET3 - (kdeg + kRET2RBC)* delta_RET3   - kon*alpha_RET3*delta_RET3   + koffalphadelta      *VRET3* alphadelta_RET3;

dxdt_alpha_RET4   = ksynalpha                      * RET4 - (kdeg + kRET2RBC)* alpha_RET4   - (kon*alpha_RET4*beta_RET4 + kon*alpha_RET4*betanew_RET4 + kon*alpha_RET4*gamma_RET4 + kon*alpha_RET4*delta_RET4) + (koffalphabetasickle*alphabeta_RET4 + koffalphabetatrans*alphabetanew_RET4 + koffalphagamma*alphagamma_RET4 + koffalphadelta*alphadelta_RET4)*VRET4;
dxdt_beta_RET4    = ksynbeta * (1-ratiosynbetanew) * RET4 - (kdeg + kRET2RBC)* beta_RET4    - kon*alpha_RET4*beta_RET4    + koffalphabetasickle *VRET4* alphabeta_RET4; 
dxdt_betanew_RET4 = ksynbeta * ratiosynbetanew     * RET4 - (kdeg + kRET2RBC)* betanew_RET4 - kon*alpha_RET4*betanew_RET4 + koffalphabetatrans  *VRET4* alphabetanew_RET4;
dxdt_gamma_RET4   = ksyngamma                      * RET4 - (kdeg + kRET2RBC)* gamma_RET4   - kon*alpha_RET4*gamma_RET4   + koffalphagamma      *VRET4* alphagamma_RET4;
dxdt_delta_RET4   = ksyndelta                      * RET4 - (kdeg + kRET2RBC)* delta_RET4   - kon*alpha_RET4*delta_RET4   + koffalphadelta      *VRET4* alphadelta_RET4;

// dimers in transduced RET
dxdt_alphabeta_RET1    = kon * alpha_RET1 * beta_RET1    - koffalphabetasickle *VRET1* alphabeta_RET1    - 2*kon* alphabeta_RET1    * alphabeta_RET1    + 2*koffHbS *VRET1* HbS_RET1  - kRET2RBC * alphabeta_RET1;
dxdt_alphabetanew_RET1 = kon * alpha_RET1 * betanew_RET1 - koffalphabetatrans  *VRET1* alphabetanew_RET1 - 2*kon* alphabetanew_RET1 * alphabetanew_RET1 + 2*koffHbA *VRET1* HbA_RET1  - kRET2RBC * alphabetanew_RET1;
dxdt_alphagamma_RET1   = kon * alpha_RET1 * gamma_RET1   - koffalphagamma      *VRET1* alphagamma_RET1   - 2*kon* alphagamma_RET1   * alphagamma_RET1   + 2*koffHbF *VRET1* HbF_RET1  - kRET2RBC * alphagamma_RET1; 
dxdt_alphadelta_RET1   = kon * alpha_RET1 * delta_RET1   - koffalphadelta      *VRET1* alphadelta_RET1   - 2*kon* alphadelta_RET1   * alphadelta_RET1   + 2*koffHbA2*VRET1* HbA2_RET1 - kRET2RBC * alphadelta_RET1; 

dxdt_alphabeta_RET2    = kon * alpha_RET2 * beta_RET2    - koffalphabetasickle *VRET2* alphabeta_RET2    - 2*kon* alphabeta_RET2    * alphabeta_RET2    + 2*koffHbS *VRET2* HbS_RET2  - kRET2RBC * alphabeta_RET2;
dxdt_alphabetanew_RET2 = kon * alpha_RET2 * betanew_RET2 - koffalphabetatrans  *VRET2* alphabetanew_RET2 - 2*kon* alphabetanew_RET2 * alphabetanew_RET2 + 2*koffHbA *VRET2* HbA_RET2  - kRET2RBC * alphabetanew_RET2;
dxdt_alphagamma_RET2   = kon * alpha_RET2 * gamma_RET2   - koffalphagamma      *VRET2* alphagamma_RET2   - 2*kon* alphagamma_RET2   * alphagamma_RET2   + 2*koffHbF *VRET2* HbF_RET2  - kRET2RBC * alphagamma_RET2; 
dxdt_alphadelta_RET2   = kon * alpha_RET2 * delta_RET2   - koffalphadelta      *VRET2* alphadelta_RET2   - 2*kon* alphadelta_RET2   * alphadelta_RET2   + 2*koffHbA2*VRET2* HbA2_RET2 - kRET2RBC * alphadelta_RET2; 

dxdt_alphabeta_RET3    = kon * alpha_RET3 * beta_RET3    - koffalphabetasickle *VRET3* alphabeta_RET3    - 2*kon* alphabeta_RET3    * alphabeta_RET3    + 2*koffHbS *VRET3* HbS_RET3  - kRET2RBC * alphabeta_RET3;
dxdt_alphabetanew_RET3 = kon * alpha_RET3 * betanew_RET3 - koffalphabetatrans  *VRET3* alphabetanew_RET3 - 2*kon* alphabetanew_RET3 * alphabetanew_RET3 + 2*koffHbA *VRET3* HbA_RET3  - kRET2RBC * alphabetanew_RET3;
dxdt_alphagamma_RET3   = kon * alpha_RET3 * gamma_RET3   - koffalphagamma      *VRET3* alphagamma_RET3   - 2*kon* alphagamma_RET3   * alphagamma_RET3   + 2*koffHbF *VRET3* HbF_RET3  - kRET2RBC * alphagamma_RET3; 
dxdt_alphadelta_RET3   = kon * alpha_RET3 * delta_RET3   - koffalphadelta      *VRET3* alphadelta_RET3   - 2*kon* alphadelta_RET3   * alphadelta_RET3   + 2*koffHbA2*VRET3* HbA2_RET3 - kRET2RBC * alphadelta_RET3; 

dxdt_alphabeta_RET4    = kon * alpha_RET4 * beta_RET4    - koffalphabetasickle *VRET4* alphabeta_RET4    - 2*kon* alphabeta_RET4    * alphabeta_RET4    + 2*koffHbS *VRET4* HbS_RET4  - kRET2RBC * alphabeta_RET4;
dxdt_alphabetanew_RET4 = kon * alpha_RET4 * betanew_RET4 - koffalphabetatrans  *VRET4* alphabetanew_RET4 - 2*kon* alphabetanew_RET4 * alphabetanew_RET4 + 2*koffHbA *VRET4* HbA_RET4  - kRET2RBC * alphabetanew_RET4;
dxdt_alphagamma_RET4   = kon * alpha_RET4 * gamma_RET4   - koffalphagamma      *VRET4* alphagamma_RET4   - 2*kon* alphagamma_RET4   * alphagamma_RET4   + 2*koffHbF *VRET4* HbF_RET4  - kRET2RBC * alphagamma_RET4; 
dxdt_alphadelta_RET4   = kon * alpha_RET4 * delta_RET4   - koffalphadelta      *VRET4* alphadelta_RET4   - 2*kon* alphadelta_RET4   * alphadelta_RET4   + 2*koffHbA2*VRET4* HbA2_RET4 - kRET2RBC * alphadelta_RET4; 

// tetramer in transduced RET
dxdt_HbS_RET1  = kon * alphabeta_RET1    * alphabeta_RET1    - (koffHbS *VRET1 + kRET2RBC) * HbS_RET1;
dxdt_HbA_RET1  = kon * alphabetanew_RET1 * alphabetanew_RET1 - (koffHbA *VRET1 + kRET2RBC) * HbA_RET1;
dxdt_HbF_RET1  = kon * alphagamma_RET1   * alphagamma_RET1   - (koffHbF *VRET1 + kRET2RBC) * HbF_RET1;
dxdt_HbA2_RET1 = kon * alphadelta_RET1   * alphadelta_RET1   - (koffHbA2*VRET1 + kRET2RBC) * HbA2_RET1;

dxdt_HbS_RET2  = kon * alphabeta_RET2    * alphabeta_RET2    - (koffHbS *VRET2 + kRET2RBC) * HbS_RET2;
dxdt_HbA_RET2  = kon * alphabetanew_RET2 * alphabetanew_RET2 - (koffHbA *VRET2 + kRET2RBC) * HbA_RET2;
dxdt_HbF_RET2  = kon * alphagamma_RET2   * alphagamma_RET2   - (koffHbF *VRET2 + kRET2RBC) * HbF_RET2;
dxdt_HbA2_RET2 = kon * alphadelta_RET2   * alphadelta_RET2   - (koffHbA2*VRET2 + kRET2RBC) * HbA2_RET2;

dxdt_HbS_RET3  = kon * alphabeta_RET3    * alphabeta_RET3    - (koffHbS *VRET3 + kRET2RBC) * HbS_RET3;
dxdt_HbA_RET3  = kon * alphabetanew_RET3 * alphabetanew_RET3 - (koffHbA *VRET3 + kRET2RBC) * HbA_RET3;
dxdt_HbF_RET3  = kon * alphagamma_RET3   * alphagamma_RET3   - (koffHbF *VRET3 + kRET2RBC) * HbF_RET3;
dxdt_HbA2_RET3 = kon * alphadelta_RET3   * alphadelta_RET3   - (koffHbA2*VRET3 + kRET2RBC) * HbA2_RET3;

dxdt_HbS_RET4  = kon * alphabeta_RET4    * alphabeta_RET4    - (koffHbS *VRET4 + kRET2RBC) * HbS_RET4;
dxdt_HbA_RET4  = kon * alphabetanew_RET4 * alphabetanew_RET4 - (koffHbA *VRET4 + kRET2RBC) * HbA_RET4;
dxdt_HbF_RET4  = kon * alphagamma_RET4   * alphagamma_RET4   - (koffHbF *VRET4 + kRET2RBC) * HbF_RET4;
dxdt_HbA2_RET4 = kon * alphadelta_RET4   * alphadelta_RET4   - (koffHbA2*VRET4 + kRET2RBC) * HbA2_RET4;

// monomer in transduced RBC
dxdt_alpha_RBC1   = kRET2RBC * alpha_RET1   - (kdeg + kdeathRBCtrans) * alpha_RBC1 - (kon*alpha_RBC1*beta_RBC1 + kon*alpha_RBC1*betanew_RBC1 + kon*alpha_RBC1*gamma_RBC1 + kon*alpha_RBC1*delta_RBC1) + (koffalphabetasickle * alphabeta_RBC1 + koffalphabetatrans * alphabetanew_RBC1 + koffalphagamma * alphagamma_RBC1 + koffalphadelta * alphadelta_RBC1)*VRBC1; 
dxdt_beta_RBC1    = kRET2RBC * beta_RET1    - (kdeg + kdeathRBCtrans) * beta_RBC1    - kon * alpha_RBC1 * beta_RBC1    + koffalphabetasickle *VRBC1* alphabeta_RBC1;
dxdt_betanew_RBC1 = kRET2RBC * betanew_RET1 - (kdeg + kdeathRBCtrans) * betanew_RBC1 - kon * alpha_RBC1 * betanew_RBC1 + koffalphabetatrans  *VRBC1* alphabetanew_RBC1;
dxdt_gamma_RBC1   = kRET2RBC * gamma_RET1   - (kdeg + kdeathRBCtrans) * gamma_RBC1   - kon * alpha_RBC1 * gamma_RBC1   + koffalphagamma      *VRBC1* alphagamma_RBC1;
dxdt_delta_RBC1   = kRET2RBC * delta_RET1   - (kdeg + kdeathRBCtrans) * delta_RBC1   - kon * alpha_RBC1 * delta_RBC1   + koffalphadelta      *VRBC1* alphadelta_RBC1;

dxdt_alpha_RBC2   = kRET2RBC * alpha_RET2   - (kdeg + kdeathRBCtrans) * alpha_RBC2 - (kon*alpha_RBC2*beta_RBC2 + kon*alpha_RBC2*betanew_RBC2 + kon*alpha_RBC2*gamma_RBC2 + kon*alpha_RBC2*delta_RBC2) + (koffalphabetasickle * alphabeta_RBC2 + koffalphabetatrans * alphabetanew_RBC2 + koffalphagamma * alphagamma_RBC2 + koffalphadelta * alphadelta_RBC2)*VRBC2; 
dxdt_beta_RBC2    = kRET2RBC * beta_RET2    - (kdeg + kdeathRBCtrans) * beta_RBC2    - kon * alpha_RBC2 * beta_RBC2    + koffalphabetasickle *VRBC2* alphabeta_RBC2;
dxdt_betanew_RBC2 = kRET2RBC * betanew_RET2 - (kdeg + kdeathRBCtrans) * betanew_RBC2 - kon * alpha_RBC2 * betanew_RBC2 + koffalphabetatrans  *VRBC2* alphabetanew_RBC2;
dxdt_gamma_RBC2   = kRET2RBC * gamma_RET2   - (kdeg + kdeathRBCtrans) * gamma_RBC2   - kon * alpha_RBC2 * gamma_RBC2   + koffalphagamma      *VRBC2* alphagamma_RBC2;
dxdt_delta_RBC2   = kRET2RBC * delta_RET2   - (kdeg + kdeathRBCtrans) * delta_RBC2   - kon * alpha_RBC2 * delta_RBC2   + koffalphadelta      *VRBC2* alphadelta_RBC2;

dxdt_alpha_RBC3   = kRET2RBC * alpha_RET3   - (kdeg + kdeathRBCtrans) * alpha_RBC3 - (kon*alpha_RBC3*beta_RBC3 + kon*alpha_RBC3*betanew_RBC3 + kon*alpha_RBC3*gamma_RBC3 + kon*alpha_RBC3*delta_RBC3) + (koffalphabetasickle * alphabeta_RBC3 + koffalphabetatrans * alphabetanew_RBC3 + koffalphagamma * alphagamma_RBC3 + koffalphadelta * alphadelta_RBC3)*VRBC3; 
dxdt_beta_RBC3    = kRET2RBC * beta_RET3    - (kdeg + kdeathRBCtrans) * beta_RBC3    -    kon * alpha_RBC3 * beta_RBC3 + koffalphabetasickle * VRBC3 * alphabeta_RBC3;
dxdt_betanew_RBC3 = kRET2RBC * betanew_RET3 - (kdeg + kdeathRBCtrans) * betanew_RBC3 - kon * alpha_RBC3 * betanew_RBC3 + koffalphabetatrans * VRBC3 * alphabetanew_RBC3;
dxdt_gamma_RBC3   = kRET2RBC * gamma_RET3   - (kdeg + kdeathRBCtrans) * gamma_RBC3   -   kon * alpha_RBC3 * gamma_RBC3 + koffalphagamma * VRBC3 * alphagamma_RBC3;
dxdt_delta_RBC3   = kRET2RBC * delta_RET3   - (kdeg + kdeathRBCtrans) * delta_RBC3   - kon * alpha_RBC3 * delta_RBC3 + koffalphadelta * VRBC3 * alphadelta_RBC3;

dxdt_alpha_RBC4   = kRET2RBC * alpha_RET4   - (kdeg + kdeathRBCtrans) * alpha_RBC4 - (kon*alpha_RBC4*beta_RBC4 + kon*alpha_RBC4*betanew_RBC4 + kon*alpha_RBC4*gamma_RBC4 + kon*alpha_RBC4*delta_RBC4) + (koffalphabetasickle * alphabeta_RBC4 + koffalphabetatrans * alphabetanew_RBC4 + koffalphagamma * alphagamma_RBC4 + koffalphadelta * alphadelta_RBC4)*VRBC4; 
dxdt_beta_RBC4    = kRET2RBC * beta_RET4    - (kdeg + kdeathRBCtrans) * beta_RBC4    - kon * alpha_RBC4 * beta_RBC4    + koffalphabetasickle * VRBC4 * alphabeta_RBC4;
dxdt_betanew_RBC4 = kRET2RBC * betanew_RET4 - (kdeg + kdeathRBCtrans) * betanew_RBC4 - kon * alpha_RBC4 * betanew_RBC4 + koffalphabetatrans * VRBC4 * alphabetanew_RBC4;
dxdt_gamma_RBC4   = kRET2RBC * gamma_RET4   - (kdeg + kdeathRBCtrans) * gamma_RBC4   - kon * alpha_RBC4 * gamma_RBC4   + koffalphagamma * VRBC4 * alphagamma_RBC4;
dxdt_delta_RBC4   = kRET2RBC * delta_RET4   - (kdeg + kdeathRBCtrans) * delta_RBC4   - kon * alpha_RBC4 * delta_RBC4   + koffalphadelta * VRBC4 * alphadelta_RBC4;

// dimers in transduced RBC
dxdt_alphabeta_RBC1    = kRET2RBC * alphabeta_RET1    - kdeathRBCtrans * alphabeta_RBC1    + kon * alpha_RBC1 * beta_RBC1    - koffalphabetasickle *VRBC1* alphabeta_RBC1    - 2*kon* alphabeta_RBC1    *alphabeta_RBC1    + 2*koffHbS *VRBC1* HbS_RBC1;
dxdt_alphabetanew_RBC1 = kRET2RBC * alphabetanew_RET1 - kdeathRBCtrans * alphabetanew_RBC1 + kon * alpha_RBC1 * betanew_RBC1 - koffalphabetatrans  *VRBC1* alphabetanew_RBC1 - 2*kon* alphabetanew_RBC1 *alphabetanew_RBC1 + 2*koffHbA *VRBC1* HbA_RBC1;
dxdt_alphagamma_RBC1   = kRET2RBC * alphagamma_RET1   - kdeathRBCtrans * alphagamma_RBC1   + kon * alpha_RBC1 * gamma_RBC1   - koffalphagamma      *VRBC1* alphagamma_RBC1   - 2*kon* alphagamma_RBC1   *alphagamma_RBC1   + 2*koffHbF *VRBC1* HbF_RBC1;
dxdt_alphadelta_RBC1   = kRET2RBC * alphadelta_RET1   - kdeathRBCtrans * alphadelta_RBC1   + kon * alpha_RBC1 * delta_RBC1   - koffalphadelta      *VRBC1* alphadelta_RBC1   - 2*kon* alphadelta_RBC1   *alphadelta_RBC1   + 2*koffHbA2*VRBC1* HbA2_RBC1;

dxdt_alphabeta_RBC2    = kRET2RBC * alphabeta_RET2    - kdeathRBCtrans * alphabeta_RBC2    + kon * alpha_RBC2 * beta_RBC2    - koffalphabetasickle *VRBC2* alphabeta_RBC2    - 2*kon* alphabeta_RBC2    *alphabeta_RBC2    + 2*koffHbS *VRBC2* HbS_RBC2;
dxdt_alphabetanew_RBC2 = kRET2RBC * alphabetanew_RET2 - kdeathRBCtrans * alphabetanew_RBC2 + kon * alpha_RBC2 * betanew_RBC2 - koffalphabetatrans  *VRBC2* alphabetanew_RBC2 - 2*kon* alphabetanew_RBC2 *alphabetanew_RBC2 + 2*koffHbA *VRBC2* HbA_RBC2;
dxdt_alphagamma_RBC2   = kRET2RBC * alphagamma_RET2   - kdeathRBCtrans * alphagamma_RBC2   + kon * alpha_RBC2 * gamma_RBC2   - koffalphagamma      *VRBC2* alphagamma_RBC2   - 2*kon* alphagamma_RBC2   *alphagamma_RBC2   + 2*koffHbF *VRBC2* HbF_RBC2;
dxdt_alphadelta_RBC2   = kRET2RBC * alphadelta_RET2   - kdeathRBCtrans * alphadelta_RBC2   + kon * alpha_RBC2 * delta_RBC2   - koffalphadelta      *VRBC2* alphadelta_RBC2   - 2*kon* alphadelta_RBC2   *alphadelta_RBC2   + 2*koffHbA2*VRBC2* HbA2_RBC2;

dxdt_alphabeta_RBC3    = kRET2RBC * alphabeta_RET3    - kdeathRBCtrans * alphabeta_RBC3    + kon * alpha_RBC3 * beta_RBC3    - koffalphabetasickle *VRBC3* alphabeta_RBC3    - 2*kon* alphabeta_RBC3    *alphabeta_RBC3    + 2*koffHbS *VRBC3* HbS_RBC3;
dxdt_alphabetanew_RBC3 = kRET2RBC * alphabetanew_RET3 - kdeathRBCtrans * alphabetanew_RBC3 + kon * alpha_RBC3 * betanew_RBC3 - koffalphabetatrans  *VRBC3* alphabetanew_RBC3 - 2*kon* alphabetanew_RBC3 *alphabetanew_RBC3 + 2*koffHbA *VRBC3* HbA_RBC3;
dxdt_alphagamma_RBC3   = kRET2RBC * alphagamma_RET3   - kdeathRBCtrans * alphagamma_RBC3   + kon * alpha_RBC3 * gamma_RBC3   - koffalphagamma      *VRBC3* alphagamma_RBC3   - 2*kon* alphagamma_RBC3   *alphagamma_RBC3   + 2*koffHbF *VRBC3* HbF_RBC3;
dxdt_alphadelta_RBC3   = kRET2RBC * alphadelta_RET3   - kdeathRBCtrans * alphadelta_RBC3   + kon * alpha_RBC3 * delta_RBC3   - koffalphadelta      *VRBC3* alphadelta_RBC3   - 2*kon* alphadelta_RBC3   *alphadelta_RBC3   + 2*koffHbA2*VRBC3* HbA2_RBC3;

dxdt_alphabeta_RBC4    = kRET2RBC * alphabeta_RET4    - kdeathRBCtrans * alphabeta_RBC4    + kon * alpha_RBC4 * beta_RBC4    - koffalphabetasickle *VRBC4* alphabeta_RBC4    - 2*kon* alphabeta_RBC4    *alphabeta_RBC4    + 2*koffHbS* VRBC4* HbS_RBC4;
dxdt_alphabetanew_RBC4 = kRET2RBC * alphabetanew_RET4 - kdeathRBCtrans * alphabetanew_RBC4 + kon * alpha_RBC4 * betanew_RBC4 - koffalphabetatrans  *VRBC4* alphabetanew_RBC4 - 2*kon* alphabetanew_RBC4 *alphabetanew_RBC4 + 2*koffHbA* VRBC4* HbA_RBC4;
dxdt_alphagamma_RBC4   = kRET2RBC * alphagamma_RET4   - kdeathRBCtrans * alphagamma_RBC4   + kon * alpha_RBC4 * gamma_RBC4   - koffalphagamma      *VRBC4* alphagamma_RBC4   - 2*kon* alphagamma_RBC4   *alphagamma_RBC4   + 2*koffHbF* VRBC4* HbF_RBC4;
dxdt_alphadelta_RBC4   = kRET2RBC * alphadelta_RET4   - kdeathRBCtrans * alphadelta_RBC4   + kon * alpha_RBC4 * delta_RBC4   - koffalphadelta      *VRBC4* alphadelta_RBC4   - 2*kon* alphadelta_RBC4   *alphadelta_RBC4   + 2*koffHbA2*VRBC4* HbA2_RBC4;

// tetramer in transduced RBC
dxdt_HbS_RBC1  = kRET2RBC * HbS_RET1  - kdeathRBCtrans * HbS_RBC1  + kon* alphabeta_RBC1   *alphabeta_RBC1    - koffHbS *VRBC1* HbS_RBC1;
dxdt_HbA_RBC1  = kRET2RBC * HbA_RET1  - kdeathRBCtrans * HbA_RBC1  + kon* alphabetanew_RBC1*alphabetanew_RBC1 - koffHbA *VRBC1* HbA_RBC1;
dxdt_HbF_RBC1  = kRET2RBC * HbF_RET1  - kdeathRBCtrans * HbF_RBC1  + kon* alphagamma_RBC1  *alphagamma_RBC1   - koffHbF *VRBC1* HbF_RBC1; 
dxdt_HbA2_RBC1 = kRET2RBC * HbA2_RET1 - kdeathRBCtrans * HbA2_RBC1 + kon* alphadelta_RBC1  *alphadelta_RBC1   - koffHbA2*VRBC1* HbA2_RBC1; 

dxdt_HbS_RBC2  = kRET2RBC * HbS_RET2  - kdeathRBCtrans * HbS_RBC2  + kon* alphabeta_RBC2   *alphabeta_RBC2    - koffHbS *VRBC2* HbS_RBC2;
dxdt_HbA_RBC2  = kRET2RBC * HbA_RET2  - kdeathRBCtrans * HbA_RBC2  + kon* alphabetanew_RBC2*alphabetanew_RBC2 - koffHbA *VRBC2* HbA_RBC2;
dxdt_HbF_RBC2  = kRET2RBC * HbF_RET2  - kdeathRBCtrans * HbF_RBC2  + kon* alphagamma_RBC2  *alphagamma_RBC2   - koffHbF *VRBC2* HbF_RBC2; 
dxdt_HbA2_RBC2 = kRET2RBC * HbA2_RET2 - kdeathRBCtrans * HbA2_RBC2 + kon* alphadelta_RBC2  *alphadelta_RBC2   - koffHbA2*VRBC2* HbA2_RBC2; 

dxdt_HbS_RBC3  = kRET2RBC * HbS_RET3  - kdeathRBCtrans * HbS_RBC3  + kon* alphabeta_RBC3   *alphabeta_RBC3    - koffHbS *VRBC3* HbS_RBC3;
dxdt_HbA_RBC3  = kRET2RBC * HbA_RET3  - kdeathRBCtrans * HbA_RBC3  + kon* alphabetanew_RBC3*alphabetanew_RBC3 - koffHbA *VRBC3* HbA_RBC3;
dxdt_HbF_RBC3  = kRET2RBC * HbF_RET3  - kdeathRBCtrans * HbF_RBC3  + kon* alphagamma_RBC3  *alphagamma_RBC3   - koffHbF *VRBC3* HbF_RBC3; 
dxdt_HbA2_RBC3 = kRET2RBC * HbA2_RET3 - kdeathRBCtrans * HbA2_RBC3 + kon* alphadelta_RBC3  *alphadelta_RBC3   - koffHbA2*VRBC3* HbA2_RBC3; 

dxdt_HbS_RBC4  = kRET2RBC * HbS_RET4  - kdeathRBCtrans * HbS_RBC4  + kon* alphabeta_RBC4   *alphabeta_RBC4    - koffHbS *VRBC4* HbS_RBC4;
dxdt_HbA_RBC4  = kRET2RBC * HbA_RET4  - kdeathRBCtrans * HbA_RBC4  + kon* alphabetanew_RBC4*alphabetanew_RBC4 - koffHbA *VRBC4* HbA_RBC4;
dxdt_HbF_RBC4  = kRET2RBC * HbF_RET4  - kdeathRBCtrans * HbF_RBC4  + kon* alphagamma_RBC4  *alphagamma_RBC4   - koffHbF *VRBC4* HbF_RBC4; 
dxdt_HbA2_RBC4 = kRET2RBC * HbA2_RET4 - kdeathRBCtrans * HbA2_RBC4 + kon* alphadelta_RBC4  *alphadelta_RBC4   - koffHbA2*VRBC4* HbA2_RBC4; 




[TABLE]

// blood cell numbers
capture LTHSC = LT0 + LT1;
capture STHSC = ST01 + ST02 + ST11 + ST12 + ST21 + ST22; 
capture MPP = MPP01 + MPP02 + MPP11 + MPP12 + MPP21 + MPP22 + MPP31 + MPP32;
capture CMP = CMP01 + CMP02 + CMP11 + CMP12 + CMP21 + CMP22 + CMP32 + CMP32 + CMP41 + CMP42; 
capture BFUE = BFUE01 + BFUE02 + BFUE11 + BFUE12 + BFUE21 + BFUE22 + BFUE31 + BFUE32 + BFUE41 + BFUE42; 
capture CFUE = CFUE01 + CFUE02 + CFUE11 + CFUE12 + CFUE21 + CFUE22 + CFUE31 + CFUE32 + CFUE41 + CFUE42; 
capture RET = RET0 + RET1 + RET2 + RET3 + RET4; 
capture RBC = (RBC0 + RBC1 + RBC2 + RBC3 + RBC4); 
capture CD34 = LTHSC + STHSC + MPP + CMP; 

// blood cell ratio
capture RBCconc = RBC/ (bloodvolume * 1e6); // RBC per uL
capture RETconc = RET/ (bloodvolume * 1e6); // RET per uL

capture RET2RBCperc = RET/ RBC* 100; // reticulocytes, %

// hemoglobin concentration
capture totalHbA = HbA;  // unit in g/dL 
capture totalHbS = HbS; // unit in g/dL  
capture totalHbF = HbF; // unit in g/dL  
capture totalHbA2 = HbA2;  // unit in g/dL 
capture totalHb = HbA + HbS + HbF + HbA2;

// blood O2 level
capture bloodO2 = vO2;

// hemoglobin concentration in RBC
// this is to troubleshoot the hemologlobin synthesis 
capture HbinRBC = (HbA + HbS + HbF + HbA2) * (10 * bloodvolume) / (VRBC0 + VRBC1 + VRBC2 + VRBC3 + VRBC4);
// concentration of Hb inside RBC should be ~ 330g/L (from https://www.ncbi.nlm.nih.gov/books/NBK54103/ )
