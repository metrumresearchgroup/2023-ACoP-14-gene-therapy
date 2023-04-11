[PROB]

Model for blood cell dynamics healthy mouse (blood cell + Hb + lymphocyte + granulocyte)

original model adopted from Zheng et al., 2021. 
https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12638

Adjustments on blood cells: 
- ST-HSC amplification time: 1000 -> 90
- MPP amplification time: 1000 -> 450
- CMP amplification time: 16 -> 8
- BFU-E amplification time: 32 -> 16
- CFU-E residence time: 7 days -> 2 days
- RBC lifespan: 120 days -> 40 days
- RET maturation time: 3 days -> 2 days

Adjustment on Hb:
- VRET = VRBC = 0.05e-12
scaling based on Fukuda et al., 2019
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5658564/
See more details in the readme page

Lymphocytes and granulocytes dynamics is taken from Busch et al., 2015
https://www.nature.com/articles/nature14242

additional adjustment: introduction of dummy parameters to match fluxes towards lymphocyte progenitors; see comments that follow parameters

[SET]
delta = 0.1, end = 70

[GLOBAL]  
#define min(a,b)  ( a<b ? a : b) 

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

// Lymphocytes progenitors
CLP0
CLP1

// B cells in bone marrow
Boe0 // propreB cells
Bi0  // immature B cells
Boe1 // propreB cells
Bi1  // immature B cells

// spleen B cells
Bt0  // transitional B cells
BMspl0 // mature B cells in spleen
Bt1  // transitional B cells
BMspl1 // mature B cells in spleen

// mature recirculating cells
BMrec0 
BMrec1

// double negative thymocyte (4 stages)
N00
N10
N20
N30
N40
N01
N11
N21
N31
N41

// double positive thymocytes (7 stages)
P00
P10
P20
P30
P40
P50
P60
P70
P01
P11
P21
P31
P41
P51
P61
P71

// single positive CD4 thymocyte
S400
S410
S420
S401
S411
S421

// single positive CD8 thymocyte
S800
S810
S820
S801
S811
S821

// T cells in blood
cd4rec0 
cd8rec0
cd4rec1 
cd8rec1
  
// T cells in the secondary lymphoid organs
cd4lym0
cd8lym0
cd4lym1
cd8lym1
  

// myeloid 
GMP01
GMP02
GM0
GMP11
GMP12
GM1

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
aST = 90
aMPP = 450
aCMP = 8
aBFUE = 16

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

// CLP parameters
alphaMPP2CLP = 0.03 // differentiation rate, MPP -> CLP; unit day-1
betaCLP = 3 // CLP proliferation rate; unit day-1
alphaCLP2proB = 3.1 // differentiation rate, CLP -> proB; unit day-1
kCLP = 0.015 // CLP death rate, day-1 

CLP2proB_amp = 4 // dummy parameter; introduced to match the flow
CLP2DN_amp = 256 // dummy parameter; introduced to match the cell flow
alphaCLP2DN = 2.5e-4 // assumed, based on  Zlotoff and Bhandoola, 2012; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076003/

//--- B cell parameter ---//

K = 6e6  // propreB and immature B cell maximum capacity in bone marrow

// preproB cell rates (per day)
gamma = 0.3 * 4 // propreB cell self-renew
delta_oe = 0.5 * 4 // death rate of propreB cells
delta_r = 0.05 * 4 // immature B cells -> propreB cells

// immature B cell rates (per day)
mu_i = 0.1 * 4 // immature B cell death rate
delta_i_t = 0.6 * 4 // immature B cells -> transitional B cells
delta_i_re = 0.19 * 4  // immature B cells -> mature recirculating B cells 

// mature recirculating B cells (per day)
phi_s = 0.03 * 4 // mature B cell in spleen -> mature recirculting B cells
mu_re = 0.008 * 4 // mature B cell death
phi_BM = 0.94 * 4 // mature recirculating B cells -> mature B cells in spleen

// transitional B cells in spleen (per day)
mu_t = 0.03 * 4 // death rate of transitional B cells
delta_t = 0.03 * 4 // transitional B cells -> mature B cells in spleen

// mature B cells death in spleen (per day)
epsilon_spl = 0.008 * 4

//--- T cell parameters ---//

// proliferation 
pN = 0.23   // DN cells
pP = 4.5    // DP cells
pS = 0.23   // SP cells

// selection DP into SP4/ SP8 (fraction)
alpha4 = 0.06    // DP -> SP4
alpha8 = 0.01    // DP -> SP8

// differentiation parameters
alpha_muN = 0.29 // DN 
alpha_muP = 0.2029  // DP
alpha_e = 0.994   // SP4
alpha_r = 0.48   // SP8

// exponent of differentiation function
n = 127 

// death rate for DN, DP, SP
delta_dn = 0
delta_dp = 0
delta_sp = 0

// removal rate last stage DP (LP)
muLP = 0.37

// T cell entering and leaving secondary lymphoid organs
enter_lymph_cd4 = 10
enter_lymph_cd8 = 5
exit_lymph = 1.2

// death rate of T cells in blood
death_cd4b = 1/31
death_cd8b = 1/72

// death rate of T cells in other organs
death_cd4l = 0.002
death_cd8l = 0.001

// GMP and GM related dynamics
kCMP2GMP = 2.5 // differentiation rate, CMP -> GMP; unit day-1
betaGM = 3 // GM death rate; unit day-1
aGMP = 32 // GMP proliferation date, day-1
tauGMP = 0.12 //mean residence time, day-1

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

// GMP dynamics
double tau1GMP = ( log2(aGMP)-1 )/ log2(aGMP) * tauGMP;
double tau2GMP = tauGMP - tau1GMP;
double krep1GMP = ( aGMP/2-1 )/ tau1GMP;
double krep2GMP = 1/tau2GMP; 
double k12GMP = aGMP/(2*tau1GMP);
double kGMP2GM = 2/tau2GMP;

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
dxdt_MPP02 = k12MPP * MPP01 + (krep2MPP - kMPP2CMP - alphaMPP2CLP) * MPP02; // endogenous branch

dxdt_MPP11 = kST2MPP * ST12 + (krep1MPP - k12MPP) * MPP11; // transplanted branch 
dxdt_MPP12 = k12MPP * MPP11 + (krep2MPP - kMPP2CMP - alphaMPP2CLP) * MPP12; // transplanted branch 


// CMP compartment
dxdt_CMP01 = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01;  // endogenous branch
dxdt_CMP02 = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE - kCMP2GMP) * CMP02; // endogenous branch

dxdt_CMP11 = kMPP2CMP * MPP12 + (krep1CMP - k12CMP) * CMP11;  // transplanted branch 
dxdt_CMP12 = k12CMP * CMP11 + (krep2CMP - kCMP2BFUE - kCMP2GMP) * CMP12; // transplanted branch 

// GMP compartment
dxdt_GMP01 = kCMP2GMP * CMP02 + (krep1GMP-k12GMP) * GMP01;
dxdt_GMP02 = k12GMP * GMP01 + (krep2GMP - kGMP2GM) * GMP02;
dxdt_GMP11 = kCMP2GMP * CMP12 + (krep1GMP-k12GMP) * GMP11;
dxdt_GMP12 = k12GMP * GMP11 + (krep2GMP - kGMP2GM) * GMP12;

// GM compartment
dxdt_GM0 = kGMP2GM * GMP02 - betaGM * GM0; 
dxdt_GM1 = kGMP2GM * GMP12 - betaGM * GM1; 


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


// CLP compartment
dxdt_CLP0 = alphaMPP2CLP * MPP02 - kCLP * CLP0 + betaCLP * CLP0 - alphaCLP2proB * CLP0 - alphaCLP2DN * CLP0; 
dxdt_CLP1 = alphaMPP2CLP * MPP12 - kCLP * CLP1 + betaCLP * CLP1 - alphaCLP2proB * CLP1 - alphaCLP2DN * CLP1; 


// B cell models; assume CLP -> propreB has another round of division
dxdt_Boe0 = CLP2proB_amp * alphaCLP2proB * CLP0 + ( gamma * ( 1 -  (Boe0 + BMrec0 + Boe1 + BMrec1)/K ) - delta_oe ) * Boe0 + delta_r * Bi0; 
dxdt_Boe1 = CLP2proB_amp * alphaCLP2proB * CLP1 + ( gamma * ( 1 -  (Boe0 + BMrec0 + Boe1 + BMrec1)/K ) - delta_oe ) * Boe1 + delta_r * Bi1; 

dxdt_Bi0 = delta_oe * Boe0 - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi0;
dxdt_Bi1 = delta_oe * Boe1 - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi1;

dxdt_BMrec0 = delta_i_re * Bi0 + phi_s * BMspl0 - (mu_re + phi_BM) * BMrec0; 
dxdt_BMrec1 = delta_i_re * Bi1 + phi_s * BMspl1 - (mu_re + phi_BM) * BMrec1; 

dxdt_Bt0 = delta_i_t * Bi0 - (mu_t + delta_t) * Bt0;
dxdt_Bt1 = delta_i_t * Bi1 - (mu_t + delta_t) * Bt1;

dxdt_BMspl0 = delta_t * Bt0 + phi_BM * BMrec0 - (phi_s + epsilon_spl) * BMspl0; 
dxdt_BMspl1 = delta_t * Bt1 + phi_BM * BMrec1 - (phi_s + epsilon_spl) * BMspl1; 



// dynamics of DN cells
dxdt_N00 = CLP2DN_amp * alphaCLP2DN * CLP0 - (pN + delta_dn) * N00; 
dxdt_N10 =  2 * pN * N00 - ( pN + delta_dn +  min( pow(alpha_muN * 1, n), 100) ) * N10; 
dxdt_N20 =  2 * pN * N10 - ( pN + delta_dn +  min( pow(alpha_muN * 2, n), 100) ) * N20; 
dxdt_N30 =  2 * pN * N20 - ( pN + delta_dn +  min( pow(alpha_muN * 3, n), 100) ) * N30; 
dxdt_N40 =  2 * pN * N30 - ( pN + delta_dn +  min( pow(alpha_muN * 4, n), 100) ) * N40; 

dxdt_N01 = CLP2DN_amp * alphaCLP2DN * CLP1 - (pN + delta_dn) * N01; 
dxdt_N11 =  2 * pN * N01 - ( pN + delta_dn +  min( pow(alpha_muN * 1, n), 100) ) * N11; 
dxdt_N21 =  2 * pN * N11 - ( pN + delta_dn +  min( pow(alpha_muN * 2, n), 100) ) * N21; 
dxdt_N31 =  2 * pN * N21 - ( pN + delta_dn +  min( pow(alpha_muN * 3, n), 100) ) * N31; 
dxdt_N41 =  2 * pN * N31 - ( pN + delta_dn +  min( pow(alpha_muN * 4, n), 100) ) * N41; 


// dynamics of DP cells
double sum_mu_N0 = min( pow(alpha_muN * 1, n), 100) * N10 + min( pow(alpha_muN * 2, n), 100) * N20 + min( pow(alpha_muN * 3, n), 100) * N30 + min( pow(alpha_muN * 4, n), 100) * N40; 
double sum_mu_N1 = min( pow(alpha_muN * 1, n), 100) * N11 + min( pow(alpha_muN * 2, n), 100) * N21 + min( pow(alpha_muN * 3, n), 100) * N31 + min( pow(alpha_muN * 4, n), 100) * N41; 

dxdt_P00 = sum_mu_N0 + 2 * pN * N40 - (pP + delta_dp) * P00; 
dxdt_P10 = 2 * pP * P00 - ( pP + delta_dp + min( pow(alpha_muP * 1, n), 100) ) * P10; 
dxdt_P20 = 2 * pP * P10 - ( pP + delta_dp + min( pow(alpha_muP * 2, n), 100) ) * P20; 
dxdt_P30 = 2 * pP * P20 - ( pP + delta_dp + min( pow(alpha_muP * 3, n), 100) ) * P30; 
dxdt_P40 = 2 * pP * P30 - ( pP + delta_dp + min( pow(alpha_muP * 4, n), 100) ) * P40; 
dxdt_P50 = 2 * pP * P40 - ( pP + delta_dp + min( pow(alpha_muP * 5, n), 100) ) * P50; 
dxdt_P60 = 2 * pP * P50 - ( pP + delta_dp + min( pow(alpha_muP * 6, n), 100) ) * P60; 

dxdt_P01 = sum_mu_N1 + 2 * pN * N41 - (pP + delta_dp) * P01; 
dxdt_P11 = 2 * pP * P01 - ( pP + delta_dp + min( pow(alpha_muP * 1, n), 100) ) * P11; 
dxdt_P21 = 2 * pP * P11 - ( pP + delta_dp + min( pow(alpha_muP * 2, n), 100) ) * P21; 
dxdt_P31 = 2 * pP * P21 - ( pP + delta_dp + min( pow(alpha_muP * 3, n), 100) ) * P31; 
dxdt_P41 = 2 * pP * P31 - ( pP + delta_dp + min( pow(alpha_muP * 4, n), 100) ) * P41; 
dxdt_P51 = 2 * pP * P41 - ( pP + delta_dp + min( pow(alpha_muP * 5, n), 100) ) * P51; 
dxdt_P61 = 2 * pP * P51 - ( pP + delta_dp + min( pow(alpha_muP * 6, n), 100) ) * P61; 

double sum_mu_P0 = min( pow(alpha_muP * 1, n), 100) * P10 +  min( pow(alpha_muP * 2, n), 100) * P20 + min( pow(alpha_muP * 3, n), 100) * P30 +
  min( pow(alpha_muP * 4, n), 100) * P40 +  min( pow(alpha_muP * 5, n), 100) * P50 + min( pow(alpha_muP * 6, n), 100) * P60; 

double sum_mu_P1 = min( pow(alpha_muP * 1, n), 100) * P11 +  min( pow(alpha_muP * 2, n), 100) * P21 + min( pow(alpha_muP * 3, n), 100) * P31 +
  min( pow(alpha_muP * 4, n), 100) * P41 +  min( pow(alpha_muP * 5, n), 100) * P51 + min( pow(alpha_muP * 6, n), 100) * P61;

dxdt_P70 = sum_mu_P0 + 2 * pP * P60 - muLP * P70;
dxdt_P71 = sum_mu_P1 + 2 * pP * P61 - muLP * P71;  

// dynamics of CD4+ cells
dxdt_S400 = alpha4 * muLP * P70 - (pS + delta_sp) * S400; 
dxdt_S410 = 2 * pS * S400 - (pS + delta_sp + min( pow(alpha_e * 1, n), 100) ) * S410; 
dxdt_S420 = 2 * pS * S410 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S420; 

dxdt_S401 = alpha4 * muLP * P71 - (pS + delta_sp) * S401; 
dxdt_S411 = 2 * pS * S401 - (pS + delta_sp + min( pow(alpha_e * 1, n), 100) ) * S411; 
dxdt_S421 = 2 * pS * S411 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S421; 

double sum_SP40 =  min( pow(alpha_e * 1, n), 100) * S410 +  min( pow(alpha_e * 2, n), 100) * S420; 
double sum_SP41 =  min( pow(alpha_e * 1, n), 100) * S411 +  min( pow(alpha_e * 2, n), 100) * S421; 


// dynamics of CD8+ cells

dxdt_S800 = alpha8 * muLP * P70 - (pS + delta_sp) * S800; 
dxdt_S810 = 2 * pS * S800 - (pS + delta_sp +  min( pow(alpha_e * 1, n), 100) ) * S810; 
dxdt_S820 = 2 * pS * S810 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S820; 

dxdt_S801 = alpha8 * muLP * P71 - (pS + delta_sp) * S801; 
dxdt_S811 = 2 * pS * S801 - (pS + delta_sp +  min( pow(alpha_e * 1, n), 100) ) * S811; 
dxdt_S821 = 2 * pS * S811 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S821; 

double sum_SP80 =  min( pow(alpha_e * 1, n), 100) * S810 +  min( pow(alpha_e * 2, n), 100) * S820; 
double sum_SP81 =  min( pow(alpha_e * 1, n), 100) * S811 +  min( pow(alpha_e * 2, n), 100) * S821; 

// dynamics of T cells in blood
dxdt_cd4rec0 = sum_SP40 - death_cd4b * cd4rec0 + exit_lymph * cd4lym0 - enter_lymph_cd4 * cd4rec0;
dxdt_cd8rec0 = sum_SP80 - death_cd8b * cd8rec0 + exit_lymph * cd8lym0 - enter_lymph_cd8 * cd8rec0;

dxdt_cd4rec1 = sum_SP41 - death_cd4b * cd4rec1 + exit_lymph * cd4lym1 - enter_lymph_cd4 * cd4rec1;
dxdt_cd8rec1 = sum_SP81 - death_cd8b * cd8rec1 + exit_lymph * cd8lym1 - enter_lymph_cd8 * cd8rec1;

// dynamics of T cells in secondary lymphoid system
dxdt_cd4lym0 = enter_lymph_cd4 * cd4rec0 - (exit_lymph + death_cd4l) * cd4lym0;
dxdt_cd8lym0 = enter_lymph_cd8 * cd8rec0 - (exit_lymph + death_cd8l) * cd8lym0; 

dxdt_cd4lym1 = enter_lymph_cd4 * cd4rec1 - (exit_lymph + death_cd4l) * cd4lym1;
dxdt_cd8lym1 = enter_lymph_cd8 * cd8rec1 - (exit_lymph + death_cd8l) * cd8lym1; 


[TABLE]
// hemoglobin 
capture totalHb = HbA + HbF + HbA2; // unit: g/dL
capture HbinRBC = (HbA_RBC0 + HbA_RBC1 + HbF_RBC0 + HbF_RBC1 + HbA2_RBC0 + HbA2_RBC1) * MWHh*1e-9 / (VRBC0 + VRBC1); // unit: g/L

// fluxes
capture CLP2proB0 = CLP2proB_amp * alphaCLP2proB * CLP0 + CLP2proB_amp * alphaCLP2proB * CLP1; // CLP -> propreB cells; endogenous branch
capture CLP2DN = CLP2DN_amp * alphaCLP2DN * CLP0; 
capture CLPexport2thymus = alphaCLP2DN * CLP0; 

// intermediate progenitors
capture ST = ST01 + ST02 + ST11 + ST12;
capture MPP = MPP01 + MPP02 + MPP11 + MPP12; 
capture CMP = CMP01 + CMP02 + CMP11 + CMP12; 
capture GMP = GMP01 + GMP02 + GMP11 + GMP12;
capture DN = N00 + N10 + N20 + N30 + N40 + N10 + N11 + N21 + N31 + N41; 
capture DP = P00 + P10 + P20 + P30 + P40 + P50 + P60 + P70 + P01 + P11 + P21 + P31 + P41 + P51 + P61 + P71; 
capture Bcell_spleen = Bt0 + BMspl0 + Bt1 + BMspl1;
capture BcellBM = Boe0 + Bi0 + Boe1 + Bi1;
capture SP4 = S400 + S410 + S420 + S401 + S411 + S421;
capture SP8 = S800 + S810 + S820 + S801 + S811 + S821;
capture thymocyte = DN + DP + SP4 + SP8;
capture Bspleen = Bt0 + Bt1 + BMspl0 + BMspl1; 
capture BM = LT0 + LT1 + ST + MPP + CMP + GMP + CLP0 + CLP1 + BFUE01 + BFUE02 + BFUE11 + BFUE12 + CFUE01 + CFUE02 + CFUE11 + CFUE12 + Boe0 + Bi0 + Boe1 + Bi1;


// end point cell count
capture GM_conc = (GM0 + GM1)/(bloodvolume * 1e6) ; // granulocyte count per uL blood
capture RBCconc = (RBC0 + RBC1)/ (bloodvolume * 1e6); // RBC per uL
capture RETconc = (RET0 + RET1)/ (bloodvolume * 1e6); // RET per uL
capture B = (BMrec0 + BMrec1)/(bloodvolume * 1e6); // B cells in periphery blood
capture T = (cd4rec0 + cd4rec1 + cd8rec0 + cd8rec1)/(bloodvolume * 1e6); // T cells in periphery blood
capture lymphocyte_conc = (cd4rec0 + cd4rec1 + cd8rec0 + cd8rec1 + BMrec0 + BMrec1)/ (bloodvolume * 1e6); // lymphocyte count per uL periphery blood
capture blood_CD4_CD8_ratio = (cd4rec0 + cd4rec1)/ (cd8rec0 + cd8rec1); 


// transduced branch
capture BMtrans = LT1 + ST11 + ST12 + MPP11 + MPP12 + CMP11 + CMP12 + BFUE11 + BFUE12 + CFUE11 + CFUE12 + CLP1 + Boe1 + Bi1;
capture thytrans = N11 + N21 + N31 + N41 + P01 + P11 + P21 + P31 + P41 + P51 + P61 + P71 + S401 + S411 + S421 + S801 + S811 + S821; 
capture DNtrans = N11 + N21 + N31 + N41;
capture DPtrans = P01 + P11 + P21 + P31 + P41 + P51 + P61 + P71;
capture SP4trans = S401 + S411 + S421;
capture SP8trans = S801 + S811 + S821;

// ratios
capture transB = BMrec1/(BMrec0 + BMrec1); 
capture transT = (cd4rec1 + cd8rec1)/(cd4rec0 + cd8rec0 + cd4rec1 + cd8rec1);


[capture]
aCFUE, sum_SP40, sum_SP80
