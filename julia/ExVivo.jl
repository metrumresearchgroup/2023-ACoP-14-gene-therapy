using DifferentialEquations, LinearAlgebra, ComponentArrays
using Parameters: @unpack
using DataFrames, CSV, DataFramesMeta, Plots

function pow(a, n)
    return a^n;
end

function HSChomo!(du, u, p,t)

    @unpack cis_erythrocyte, trans_erythrocyte, cis_hb_ret, cis_hb_rbc, trans_hb_ret, trans_hb_rbc, myeloid, cisB, transB, cisT, transT = u;

    @unpack LT0, ST01, ST02, MPP01, MPP02, CMP01, CMP02, BFUE01, BFUE02, CFUE01, CFUE02, RET0, RBC0 = cis_erythrocyte
    @unpack LT1, ST11, ST12, MPP11, MPP12, CMP11, CMP12, BFUE11, BFUE12, CFUE11, CFUE12, RET1, RBC1 = trans_erythrocyte

    @unpack alpha_RET0, beta_RET0, gamma_RET0, delta_RET0, alphabeta_RET0, alphagamma_RET0, alphadelta_RET0, HbA_RET0, HbF_RET0, HbA2_RET0 = cis_hb_ret
    @unpack alpha_RBC0, beta_RBC0, gamma_RBC0, delta_RBC0, alphabeta_RBC0, alphagamma_RBC0, alphadelta_RBC0, HbA_RBC0, HbF_RBC0, HbA2_RBC0 = cis_hb_rbc
    @unpack alpha_RET1, beta_RET1, gamma_RET1, delta_RET1, alphabeta_RET1, alphagamma_RET1, alphadelta_RET1, HbA_RET1, HbF_RET1, HbA2_RET1 = trans_hb_ret
    @unpack alpha_RBC1, beta_RBC1, gamma_RBC1, delta_RBC1, alphabeta_RBC1, alphagamma_RBC1, alphadelta_RBC1, HbA_RBC1, HbF_RBC1, HbA2_RBC1 = trans_hb_rbc

    @unpack GMP01, GMP02, GMP11, GMP12, GM0, GM1 = myeloid

    @unpack CLP0, Boe0, Bi0, BMrec0, Bt0, BMspl0 = cisB 
    @unpack CLP1, Boe1, Bi1, BMrec1, Bt1, BMspl1 = transB

    @unpack N00, N10, N20, N30, N40, P00, P10, P20, P30, P40, P50, P60, P70, S400, S410, S420, S800, S810, S820, cd4rec0, cd8rec0, cd4peripheral0, cd8peripheral0, cd4lym0, cd8lym0 = cisT
    @unpack N01, N11, N21, N31, N41, P01, P11, P21, P31, P41, P51, P61, P71, S401, S411, S421, S801, S811, S821, cd4rec1, cd8rec1, cd4peripheral1, cd8peripheral1, cd4lym1, cd8lym1 = transT

    @unpack tauLT, tauST, tauMPP, tauCMP, tauBFUE, tauCFUE, tauRET, tauRBCendo, tauRBCtrans, ssLT, aST, aMPP, aCMP, aBFUE, rLT, VRET, VRBC, MWHh, 
    ksynalpha, ratiosynbeta, ratiosyngamma, ratiosyndelta, thalfmonomer, kon, Kdalphabeta, Kdalphagamma, Kdalphadelta, KdHbA, KdHbF, KdHbA2, 
    HbA_saturation, HbF_saturation, HbA2_saturation, bloodvolume, 
    kCMP2GMP, betaGM, aGMP, tauGMP, 
    alphaMPP2CLP, betaCLP, alphaCLP2proB, kCLP, CLP2proB_amp, K, gamma, delta_oe, delta_r, mu_i, delta_i_t, delta_i_re, phi_s, mu_re, phi_BM, mu_t, delta_t, epsilon_spl, epsilon_spl_r, 
    alphaCLP2DN, CLP2DN_amp, pN, pP, pS, alpha4, alpha8, alpha_muN, alpha_muP, alpha_e, alpha_r, delta_dn, delta_dp, delta_dp_r, delta_sp, muLP, death_naiveT_cd8, death_naiveT_cd4, 
    enter_lymph, enter_peripheral, exit_lymph, exit_peripheral_cd4, exit_peripheral_cd8, env_naiveT, k_nT_pro, peripheralvolume = p;

    # LT-HSC dynamics
    kLT2ST = 1/tauLT;
    deltaLT = (rLT - kLT2ST)/ssLT; 
    # ST-HSC dynamics
    tau1ST = (log2(aST) -1)/log2(aST) * tauST;
    tau2ST = tauST - tau1ST;
    krep1ST = (aST/2-1)/ tau1ST;
    krep2ST = 1/ tau2ST;
    k12ST = aST/(2*tau1ST);
    kST2MPP = 2/tau2ST;
    # MPP dynamics
    tau1MPP = (log2(aMPP)-1)/log2(aMPP) * tauMPP;
    tau2MPP = tauMPP - tau1MPP;
    krep1MPP = (aMPP/2-1)/ tau1MPP;
    krep2MPP = 1/ tau2MPP;
    k12MPP = aMPP/(2*tau1MPP);
    kMPP2CMP = 2/tau2MPP;
    # CMP dynamics
    tau1CMP = (log2(aCMP)-1)/log2(aCMP) * tauCMP;
    tau2CMP = tauCMP - tau1CMP;
    krep1CMP = (aCMP/2-1)/ tau1CMP;
    krep2CMP = 1/ tau2CMP;
    k12CMP = aCMP/(2*tau1CMP);
    kCMP2BFUE = 2/tau2CMP;
    # BFU-E dynamics
    tau1BFUE = (log2(aBFUE)-1)/log2(aBFUE) * tauBFUE;
    tau2BFUE = tauBFUE - tau1BFUE;
    krep1BFUE = (aBFUE/2-1)/ tau1BFUE;
    krep2BFUE = 1/ tau2BFUE;
    k12BFUE = aBFUE/(2*tau1BFUE);
    kBFUE2CFUE = 2/tau2BFUE;
    # RET and RBC dynamics
    kRET2RBC = 1/tauRET;
    kdeathRBCendo = 1/tauRBCendo;
    kdeathRBCtrans = 1/tauRBCtrans;
    # monomer synthesis rate
    ksynbeta = ratiosynbeta * ksynalpha;
    ksyngamma = ratiosyngamma * ksynalpha;
    ksyndelta = ratiosyndelta * ksynalpha;
    # monomer degredation rate
    kdeg = log(2)/ thalfmonomer; 
    # dimer dissociation rate
    koffalphagamma = Kdalphagamma * kon; 
    koffalphadelta = Kdalphadelta * kon; 
    koffalphabeta = Kdalphabeta * kon; 
    # tetramer dissociation rate
    koffHbA = KdHbA * kon; 
    koffHbF = KdHbF * kon; 
    koffHbA2 = KdHbA2 * kon; 
    # total RET & RBC volume 
    VRET0 = VRET * RET0 ; 
    VRBC0 = VRBC * RBC0 ;
    VRET1 = VRET * RET1 ; 
    VRBC1 = VRBC * RBC1 ;

    # Hb concentration & oxygen level
    HbA = (HbA_RET0 + HbA_RET1 + HbA_RBC0 + HbA_RBC1) * MWHh*1e-9/ (10*bloodvolume);  # unit in g/dL 
    HbF = (HbF_RET0 + HbF_RET1 + HbF_RBC0 + HbF_RBC1) * MWHh*1e-9/ (10*bloodvolume); # unit in g/dL  
    HbA2 = (HbA2_RET0 + HbA2_RET1 + HbA2_RBC0 + HbA2_RBC1) * MWHh*1e-9/ (10*bloodvolume);  # unit in g/dL 
    vO2 = ( HbA_saturation * HbA + HbF_saturation * HbF + HbA2_saturation * HbA2 ) * 1.34;
    #  CFU-E amplification time
    aCFUE = 550 * exp( -0.23 * vO2 );
    tau1CFUE = (log2(aCFUE)-1)/log2(aCFUE) * tauCFUE;
    tau2CFUE = tauCFUE - tau1CFUE;
    krep1CFUE = (aCFUE/2-1)/ tau1CFUE;
    krep2CFUE = 1/ tau2CFUE;
    k12CFUE = aCFUE/(2*tau1CFUE);
    kCFUE2RET = 2/tau2CFUE;

    # granulocyte and progenitor dynamics
    tau1GMP = ( log2(aGMP)-1 )/ log2(aGMP) * tauGMP;
    tau2GMP = tauGMP - tau1GMP;
    krep1GMP = ( aGMP/2-1 )/ tau1GMP;
    krep2GMP = 1/tau2GMP; 
    k12GMP = aGMP/(2*tau1GMP);
    kGMP2GM = 2/tau2GMP;

    # endogenous erythrocyte dynamics
    du.cis_erythrocyte.LT0    = rLT * LT0 - deltaLT * LT0 * (LT0 + LT1) - kLT2ST * LT0;
    du.cis_erythrocyte.ST01   = kLT2ST * LT0 + (krep1ST - k12ST) * ST01;
    du.cis_erythrocyte.ST02   = k12ST * ST01 + (krep2ST - kST2MPP) * ST02; 
    du.cis_erythrocyte.MPP01  = kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; 
    du.cis_erythrocyte.MPP02  = k12MPP * MPP01 + (krep2MPP - kMPP2CMP - alphaMPP2CLP) * MPP02; 
    du.cis_erythrocyte.CMP01  = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01; 
    du.cis_erythrocyte.CMP02  = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE - kCMP2GMP) * CMP02;
    du.cis_erythrocyte.BFUE01 = kCMP2BFUE * CMP02 + (krep1BFUE - k12BFUE) * BFUE01;  
    du.cis_erythrocyte.BFUE02 = k12BFUE * BFUE01 + (krep2BFUE - kBFUE2CFUE) * BFUE02; 
    du.cis_erythrocyte.CFUE01 = kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; 
    du.cis_erythrocyte.CFUE02 = k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;
    du.cis_erythrocyte.RET0   = kCFUE2RET * CFUE02 - kRET2RBC * RET0;  
    du.cis_erythrocyte.RBC0   = kRET2RBC * RET0 - kdeathRBCendo * RBC0; 

    # Hb dynamics in endogenous RET/ RBC
    flux_alpha_beta_RET0  = kon * alpha_RET0 * beta_RET0  - koffalphabeta  *VRET0 * alphabeta_RET0;
    flux_alpha_gamma_RET0 = kon * alpha_RET0 * gamma_RET0 - koffalphagamma *VRET0 * alphagamma_RET0;
    flux_alpha_delta_RET0 = kon * alpha_RET0 * delta_RET0 - koffalphadelta *VRET0 * alphadelta_RET0; 
    
    du.cis_hb_ret.alpha_RET0 = ksynalpha * RET0 - (kdeg + kRET2RBC)*alpha_RET0 - flux_alpha_beta_RET0 - flux_alpha_gamma_RET0 - flux_alpha_delta_RET0;
    du.cis_hb_ret.beta_RET0  = ksynbeta  * RET0 - (kdeg + kRET2RBC)*beta_RET0  - flux_alpha_beta_RET0; 
    du.cis_hb_ret.gamma_RET0 = ksyngamma * RET0 - (kdeg + kRET2RBC)*gamma_RET0 - flux_alpha_gamma_RET0;
    du.cis_hb_ret.delta_RET0 = ksyndelta * RET0 - (kdeg + kRET2RBC)*delta_RET0 - flux_alpha_delta_RET0; 

    flux_alphabeta_RET0  = kon * alphabeta_RET0  * alphabeta_RET0  - koffHbA *VRET0 * HbA_RET0; 
    flux_alphagamma_RET0 = kon * alphagamma_RET0 * alphagamma_RET0 - koffHbF *VRET0 * HbF_RET0;
    flux_alphadelta_RET0 = kon * alphadelta_RET0 * alphadelta_RET0 - koffHbA2*VRET0 * HbA2_RET0; 

    du.cis_hb_ret.alphabeta_RET0  = flux_alpha_beta_RET0   - kRET2RBC * alphabeta_RET0  - 2 * flux_alphabeta_RET0;
    du.cis_hb_ret.alphagamma_RET0 = flux_alpha_gamma_RET0  - kRET2RBC * alphagamma_RET0 - 2 * flux_alphagamma_RET0; 
    du.cis_hb_ret.alphadelta_RET0 = flux_alpha_delta_RET0  - kRET2RBC * alphadelta_RET0 - 2 * flux_alphadelta_RET0; 

    du.cis_hb_ret.HbA_RET0  = flux_alphabeta_RET0  - kRET2RBC * HbA_RET0;
    du.cis_hb_ret.HbF_RET0  = flux_alphagamma_RET0 - kRET2RBC * HbF_RET0;
    du.cis_hb_ret.HbA2_RET0 = flux_alphadelta_RET0 - kRET2RBC * HbA2_RET0;

    flux_alpha_beta_RBC0 = kon * alpha_RBC0 * beta_RBC0  - koffalphabeta *VRBC0 * alphabeta_RBC0;
    flux_alpha_gamma_RBC0= kon * alpha_RBC0 * gamma_RBC0 - koffalphagamma*VRBC0 * alphagamma_RBC0;
    flux_alpha_delta_RBC0= kon * alpha_RBC0 * delta_RBC0 - koffalphadelta*VRBC0 * alphadelta_RBC0;

    du.cis_hb_rbc.alpha_RBC0 = kRET2RBC * alpha_RET0 - (kdeg + kdeathRBCendo) * alpha_RBC0 - flux_alpha_beta_RBC0 - flux_alpha_gamma_RBC0 - flux_alpha_delta_RBC0;
    du.cis_hb_rbc.beta_RBC0 =  kRET2RBC * beta_RET0  - (kdeg + kdeathRBCendo) * beta_RBC0  - flux_alpha_beta_RBC0;
    du.cis_hb_rbc.gamma_RBC0 = kRET2RBC * gamma_RET0 - (kdeg + kdeathRBCendo) * gamma_RBC0 - flux_alpha_gamma_RBC0;
    du.cis_hb_rbc.delta_RBC0 = kRET2RBC * delta_RET0 - (kdeg + kdeathRBCendo) * delta_RBC0 - flux_alpha_delta_RBC0;

    du.cis_hb_rbc.alphabeta_RBC0  = kRET2RBC * alphabeta_RET0  - kdeathRBCendo * alphabeta_RBC0  + flux_alpha_beta_RBC0  - 2*kon* alphabeta_RBC0 *alphabeta_RBC0  + 2*koffHbA *VRBC0 * HbA_RBC0;
    du.cis_hb_rbc.alphagamma_RBC0 = kRET2RBC * alphagamma_RET0 - kdeathRBCendo * alphagamma_RBC0 + flux_alpha_gamma_RBC0 - 2*kon* alphagamma_RBC0*alphagamma_RBC0 + 2*koffHbF *VRBC0 * HbF_RBC0;
    du.cis_hb_rbc.alphadelta_RBC0 = kRET2RBC * alphadelta_RET0 - kdeathRBCendo * alphadelta_RBC0 + flux_alpha_delta_RBC0 - 2*kon* alphadelta_RBC0*alphadelta_RBC0 + 2*koffHbA2*VRBC0 * HbA2_RBC0;

    du.cis_hb_rbc.HbA_RBC0  = kRET2RBC * HbA_RET0  - kdeathRBCendo * HbA_RBC0  + kon*alphabeta_RBC0 *alphabeta_RBC0  - koffHbA * VRBC0 * HbA_RBC0;
    du.cis_hb_rbc.HbF_RBC0  = kRET2RBC * HbF_RET0  - kdeathRBCendo * HbF_RBC0  + kon*alphagamma_RBC0*alphagamma_RBC0 - koffHbF * VRBC0 * HbF_RBC0; 
    du.cis_hb_rbc.HbA2_RBC0 = kRET2RBC * HbA2_RET0 - kdeathRBCendo * HbA2_RBC0 + kon*alphadelta_RBC0*alphadelta_RBC0 - koffHbA2* VRBC0 * HbA2_RBC0; 

    # transduced erythrocyte dynamics
    du.trans_erythrocyte.LT1    = rLT * LT1 - deltaLT * LT1 * (LT0 + LT1) - kLT2ST * LT1;
    du.trans_erythrocyte.ST11   = kLT2ST * LT1 + (krep1ST - k12ST) * ST11; 
    du.trans_erythrocyte.ST12   = k12ST * ST11 + (krep2ST - kST2MPP) * ST12;
    du.trans_erythrocyte.MPP11  = kST2MPP * ST12 + (krep1MPP - k12MPP) * MPP11; 
    du.trans_erythrocyte.MPP12  = k12MPP * MPP11 + (krep2MPP - kMPP2CMP - alphaMPP2CLP) * MPP12;
    du.trans_erythrocyte.CMP11  = kMPP2CMP * MPP12 + (krep1CMP - k12CMP) * CMP11; 
    du.trans_erythrocyte.CMP12  = k12CMP * CMP11 + (krep2CMP - kCMP2BFUE - kCMP2GMP) * CMP12;
    du.trans_erythrocyte.BFUE11 = kCMP2BFUE * CMP12 + (krep1BFUE - k12BFUE) * BFUE11; 
    du.trans_erythrocyte.BFUE12 = k12BFUE * BFUE11 + (krep2BFUE - kBFUE2CFUE) * BFUE12; 
    du.trans_erythrocyte.CFUE11 = kBFUE2CFUE * BFUE12 + (krep1CFUE - k12CFUE) * CFUE11;
    du.trans_erythrocyte.CFUE12 = k12CFUE * CFUE11 + (krep2CFUE - kCFUE2RET) * CFUE12;
    du.trans_erythrocyte.RET1   = kCFUE2RET * CFUE12 - kRET2RBC * RET1;  
    du.trans_erythrocyte.RBC1   = kRET2RBC * RET1 - kdeathRBCtrans * RBC1;

    # Hb dynamics in transduced RET/ RBC
    flux_alpha_beta_RET1  = kon * alpha_RET1 * beta_RET1  - koffalphabeta  *VRET1 * alphabeta_RET1;
    flux_alpha_gamma_RET1 = kon * alpha_RET1 * gamma_RET1 - koffalphagamma *VRET1 * alphagamma_RET1;
    flux_alpha_delta_RET1 = kon * alpha_RET1 * delta_RET1 - koffalphadelta *VRET1 * alphadelta_RET1; 
    
    du.trans_hb_ret.alpha_RET1 = ksynalpha * RET1 - (kdeg + kRET2RBC)*alpha_RET1 - flux_alpha_beta_RET1 - flux_alpha_gamma_RET1 - flux_alpha_delta_RET1;
    du.trans_hb_ret.beta_RET1  = ksynbeta  * RET1 - (kdeg + kRET2RBC)*beta_RET1  - flux_alpha_beta_RET1; 
    du.trans_hb_ret.gamma_RET1 = ksyngamma * RET1 - (kdeg + kRET2RBC)*gamma_RET1 - flux_alpha_gamma_RET1;
    du.trans_hb_ret.delta_RET1 = ksyndelta * RET1 - (kdeg + kRET2RBC)*delta_RET1 - flux_alpha_delta_RET1; 

    flux_alphabeta_RET1  = kon * alphabeta_RET1  * alphabeta_RET1  - koffHbA *VRET1 * HbA_RET1; 
    flux_alphagamma_RET1 = kon * alphagamma_RET1 * alphagamma_RET1 - koffHbF *VRET1 * HbF_RET1;
    flux_alphadelta_RET1 = kon * alphadelta_RET1 * alphadelta_RET1 - koffHbA2*VRET1 * HbA2_RET1; 

    du.trans_hb_ret.alphabeta_RET1  = flux_alpha_beta_RET1   - kRET2RBC * alphabeta_RET1  - 2 * flux_alphabeta_RET1;
    du.trans_hb_ret.alphagamma_RET1 = flux_alpha_gamma_RET1  - kRET2RBC * alphagamma_RET1 - 2 * flux_alphagamma_RET1; 
    du.trans_hb_ret.alphadelta_RET1 = flux_alpha_delta_RET1  - kRET2RBC * alphadelta_RET1 - 2 * flux_alphadelta_RET1; 

    du.trans_hb_ret.HbA_RET1  = flux_alphabeta_RET1  - kRET2RBC * HbA_RET1;
    du.trans_hb_ret.HbF_RET1  = flux_alphagamma_RET1 - kRET2RBC * HbF_RET1;
    du.trans_hb_ret.HbA2_RET1 = flux_alphadelta_RET1 - kRET2RBC * HbA2_RET1;

    flux_alpha_beta_RBC1 = kon * alpha_RBC1 * beta_RBC1  - koffalphabeta *VRBC1 * alphabeta_RBC1;
    flux_alpha_gamma_RBC1= kon * alpha_RBC1 * gamma_RBC1 - koffalphagamma*VRBC1 * alphagamma_RBC1;
    flux_alpha_delta_RBC1= kon * alpha_RBC1 * delta_RBC1 - koffalphadelta*VRBC1 * alphadelta_RBC1;

    du.trans_hb_rbc.alpha_RBC1 = kRET2RBC * alpha_RET1 - (kdeg + kdeathRBCendo) * alpha_RBC1 - flux_alpha_beta_RBC1 - flux_alpha_gamma_RBC1 - flux_alpha_delta_RBC1;
    du.trans_hb_rbc.beta_RBC1 =  kRET2RBC * beta_RET1  - (kdeg + kdeathRBCendo) * beta_RBC1  - flux_alpha_beta_RBC1;
    du.trans_hb_rbc.gamma_RBC1 = kRET2RBC * gamma_RET1 - (kdeg + kdeathRBCendo) * gamma_RBC1 - flux_alpha_gamma_RBC1;
    du.trans_hb_rbc.delta_RBC1 = kRET2RBC * delta_RET1 - (kdeg + kdeathRBCendo) * delta_RBC1 - flux_alpha_delta_RBC1;

    du.trans_hb_rbc.alphabeta_RBC1  = kRET2RBC * alphabeta_RET1  - kdeathRBCendo * alphabeta_RBC1  + flux_alpha_beta_RBC1  - 2*kon* alphabeta_RBC1 *alphabeta_RBC1  + 2*koffHbA *VRBC1 * HbA_RBC1;
    du.trans_hb_rbc.alphagamma_RBC1 = kRET2RBC * alphagamma_RET1 - kdeathRBCendo * alphagamma_RBC1 + flux_alpha_gamma_RBC1 - 2*kon* alphagamma_RBC1*alphagamma_RBC1 + 2*koffHbF *VRBC1 * HbF_RBC1;
    du.trans_hb_rbc.alphadelta_RBC1 = kRET2RBC * alphadelta_RET1 - kdeathRBCendo * alphadelta_RBC1 + flux_alpha_delta_RBC1 - 2*kon* alphadelta_RBC1*alphadelta_RBC1 + 2*koffHbA2*VRBC1 * HbA2_RBC1;

    du.trans_hb_rbc.HbA_RBC1  = kRET2RBC * HbA_RET1  - kdeathRBCendo * HbA_RBC1  + kon*alphabeta_RBC1 *alphabeta_RBC1  - koffHbA * VRBC1 * HbA_RBC1;
    du.trans_hb_rbc.HbF_RBC1  = kRET2RBC * HbF_RET1  - kdeathRBCendo * HbF_RBC1  + kon*alphagamma_RBC1*alphagamma_RBC1 - koffHbF * VRBC1 * HbF_RBC1; 
    du.trans_hb_rbc.HbA2_RBC1 = kRET2RBC * HbA2_RET1 - kdeathRBCendo * HbA2_RBC1 + kon*alphadelta_RBC1*alphadelta_RBC1 - koffHbA2* VRBC1 * HbA2_RBC1; 

    # GMP compartment
    du.myeloid.GMP01 = kCMP2GMP * CMP02 + (krep1GMP-k12GMP) * GMP01;
    du.myeloid.GMP02 = k12GMP * GMP01 + (krep2GMP - kGMP2GM) * GMP02;
    du.myeloid.GMP11 = kCMP2GMP * CMP12 + (krep1GMP-k12GMP) * GMP11;
    du.myeloid.GMP12 = k12GMP * GMP11 + (krep2GMP - kGMP2GM) * GMP12;

    # GM compartment
    du.myeloid.GM0 = kGMP2GM * GMP02 - betaGM * GM0; 
    du.myeloid.GM1 = kGMP2GM * GMP12 - betaGM * GM1; 

    # CLP compartment
    du.cisB.CLP0 = alphaMPP2CLP * MPP02 - kCLP * CLP0 + betaCLP * CLP0 - alphaCLP2proB * CLP0 - alphaCLP2DN * CLP0; 
    du.transB.CLP1 = alphaMPP2CLP * MPP12 - kCLP * CLP1 + betaCLP * CLP1 - alphaCLP2proB * CLP1 - alphaCLP2DN * CLP1; 

    # B cell models
    du.cisB.Boe0 = CLP2proB_amp * alphaCLP2proB * CLP0 + ( gamma * ( 1 -  (Boe0 + BMrec0 + Boe1 + BMrec1)/K ) - delta_oe ) * Boe0 + delta_r * Bi0; 
    du.transB.Boe1 = CLP2proB_amp * alphaCLP2proB * CLP1 + ( gamma * ( 1 -  (Boe0 + BMrec0 + Boe1 + BMrec1)/K ) - delta_oe ) * Boe1 + delta_r * Bi1; 

    du.cisB.Bi0 = delta_oe * Boe0 - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi0;
    du.transB.Bi1 = delta_oe * Boe1 - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi1;

    du.cisB.BMrec0 = delta_i_re * Bi0 + phi_s * BMspl0 - (mu_re + phi_BM) * BMrec0; 
    du.transB.BMrec1 = delta_i_re * Bi1 + phi_s * BMspl1 - (mu_re + phi_BM) * BMrec1; 

    du.cisB.Bt0 = delta_i_t * Bi0 - (mu_t + delta_t) * Bt0;
    du.transB.Bt1 = delta_i_t * Bi1 - (mu_t + delta_t) * Bt1;

    du.cisB.BMspl0 = delta_t * Bt0 + phi_BM * BMrec0 - (phi_s + epsilon_spl) * BMspl0; 
    du.transB.BMspl1 = delta_t * Bt1 + phi_BM * BMrec1 - (phi_s + epsilon_spl_r) * BMspl1; 

    n = 127.0

    # dynamics of DN cells
    du.cisT.N00 = CLP2DN_amp * alphaCLP2DN * CLP0 - (pN + delta_dn) * N00; 
    du.cisT.N10 =  2 * pN * N00 - ( pN + delta_dn +  min( pow(alpha_muN * 1, n), 100) ) * N10; 
    du.cisT.N20 =  2 * pN * N10 - ( pN + delta_dn +  min( pow(alpha_muN * 2, n), 100) ) * N20; 
    du.cisT.N30 =  2 * pN * N20 - ( pN + delta_dn +  min( pow(alpha_muN * 3, n), 100) ) * N30; 
    du.cisT.N40 =  2 * pN * N30 - ( pN + delta_dn +  min( pow(alpha_muN * 4, n), 100) ) * N40; 

    du.transT.N01 = CLP2DN_amp * alphaCLP2DN * CLP1 - (pN + delta_dn) * N01; 
    du.transT.N11 =  2 * pN * N01 - ( pN + delta_dn +  min( pow(alpha_muN * 1, n), 100) ) * N11; 
    du.transT.N21 =  2 * pN * N11 - ( pN + delta_dn +  min( pow(alpha_muN * 2, n), 100) ) * N21; 
    du.transT.N31 =  2 * pN * N21 - ( pN + delta_dn +  min( pow(alpha_muN * 3, n), 100) ) * N31; 
    du.transT.N41 =  2 * pN * N31 - ( pN + delta_dn +  min( pow(alpha_muN * 4, n), 100) ) * N41; 

    # dynamics of DP cells
    sum_mu_N0 = min( pow(alpha_muN * 1, n), 100) * N10 + min( pow(alpha_muN * 2, n), 100) * N20 + min( pow(alpha_muN * 3, n), 100) * N30 + min( pow(alpha_muN * 4, n), 100) * N40; 
    sum_mu_N1 = min( pow(alpha_muN * 1, n), 100) * N11 + min( pow(alpha_muN * 2, n), 100) * N21 + min( pow(alpha_muN * 3, n), 100) * N31 + min( pow(alpha_muN * 4, n), 100) * N41; 

    du.cisT.P00 = sum_mu_N0 + 2 * pN * N40 - (pP + delta_dp) * P00; 
    du.cisT.P10 = 2 * pP * P00 - ( pP + delta_dp + min( pow(alpha_muP * 1, n), 100) ) * P10; 
    du.cisT.P20 = 2 * pP * P10 - ( pP + delta_dp + min( pow(alpha_muP * 2, n), 100) ) * P20; 
    du.cisT.P30 = 2 * pP * P20 - ( pP + delta_dp + min( pow(alpha_muP * 3, n), 100) ) * P30; 
    du.cisT.P40 = 2 * pP * P30 - ( pP + delta_dp + min( pow(alpha_muP * 4, n), 100) ) * P40; 
    du.cisT.P50 = 2 * pP * P40 - ( pP + delta_dp + min( pow(alpha_muP * 5, n), 100) ) * P50; 
    du.cisT.P60 = 2 * pP * P50 - ( pP + delta_dp + min( pow(alpha_muP * 6, n), 100) ) * P60; 

    du.transT.P01 = sum_mu_N1 + 2 * pN * N41 - (pP + delta_dp_r) * P01; 
    du.transT.P11 = 2 * pP * P01 - ( pP + delta_dp_r + min( pow(alpha_muP * 1, n), 100) ) * P11; 
    du.transT.P21 = 2 * pP * P11 - ( pP + delta_dp_r + min( pow(alpha_muP * 2, n), 100) ) * P21; 
    du.transT.P31 = 2 * pP * P21 - ( pP + delta_dp_r + min( pow(alpha_muP * 3, n), 100) ) * P31; 
    du.transT.P41 = 2 * pP * P31 - ( pP + delta_dp_r + min( pow(alpha_muP * 4, n), 100) ) * P41; 
    du.transT.P51 = 2 * pP * P41 - ( pP + delta_dp_r + min( pow(alpha_muP * 5, n), 100) ) * P51; 
    du.transT.P61 = 2 * pP * P51 - ( pP + delta_dp_r + min( pow(alpha_muP * 6, n), 100) ) * P61; 

    sum_mu_P0 = min( pow(alpha_muP * 1, n), 100) * P10 +  min( pow(alpha_muP * 2, n), 100) * P20 + min( pow(alpha_muP * 3, n), 100) * P30 +
                   min( pow(alpha_muP * 4, n), 100) * P40 +  min( pow(alpha_muP * 5, n), 100) * P50 + min( pow(alpha_muP * 6, n), 100) * P60; 

    sum_mu_P1 = min( pow(alpha_muP * 1, n), 100) * P11 +  min( pow(alpha_muP * 2, n), 100) * P21 + min( pow(alpha_muP * 3, n), 100) * P31 +
                   min( pow(alpha_muP * 4, n), 100) * P41 +  min( pow(alpha_muP * 5, n), 100) * P51 + min( pow(alpha_muP * 6, n), 100) * P61;

    du.cisT.P70 = sum_mu_P0 + 2 * pP * P60 - muLP * P70;
    du.transT.P71 = sum_mu_P1 + 2 * pP * P61 - muLP * P71;  

    # dynamics of CD4+ cells
    du.cisT.S400 = alpha4 * muLP * P70 - (pS + delta_sp) * S400; 
    du.cisT.S410 = 2 * pS * S400 - (pS + delta_sp + min( pow(alpha_e * 1, n), 100) ) * S410; 
    du.cisT.S420 = 2 * pS * S410 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S420; 

    du.transT.S401 = alpha4 * muLP * P71 - (pS + delta_sp) * S401; 
    du.transT.S411 = 2 * pS * S401 - (pS + delta_sp + min( pow(alpha_e * 1, n), 100) ) * S411; 
    du.transT.S421 = 2 * pS * S411 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S421; 

    sum_SP40 =  min( pow(alpha_e * 1, n), 100) * S410 +  min( pow(alpha_e * 2, n), 100) * S420; 
    sum_SP41 =  min( pow(alpha_e * 1, n), 100) * S411 +  min( pow(alpha_e * 2, n), 100) * S421; 

    # dynamics of CD8+ cells
    du.cisT.S800 = alpha8 * muLP * P70 - (pS + delta_sp) * S800; 
    du.cisT.S810 = 2 * pS * S800 - (pS + delta_sp +  min( pow(alpha_e * 1, n), 100) ) * S810; 
    du.cisT.S820 = 2 * pS * S810 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S820; 

    du.transT.S801 = alpha8 * muLP * P71 - (pS + delta_sp) * S801; 
    du.transT.S811 = 2 * pS * S801 - (pS + delta_sp +  min( pow(alpha_e * 1, n), 100) ) * S811; 
    du.transT.S821 = 2 * pS * S811 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S821; 

    sum_SP80 =  min( pow(alpha_e * 1, n), 100) * S810 +  min( pow(alpha_e * 2, n), 100) * S820; 
    sum_SP81 =  min( pow(alpha_e * 1, n), 100) * S811 +  min( pow(alpha_e * 2, n), 100) * S821; 

    # dynamics of T cells in blood
    du.cisT.cd4rec0 = sum_SP40 - death_naiveT_cd4 * cd4rec0 + exit_lymph * cd4lym0 - enter_lymph * cd4rec0 + exit_peripheral_cd4 * cd4peripheral0 - enter_peripheral * cd4rec0 + k_nT_pro/bloodvolume * (cd4rec0/bloodvolume)/(cd4rec0/bloodvolume + env_naiveT);
    du.cisT.cd8rec0 = sum_SP80 - death_naiveT_cd8 * cd8rec0 + exit_lymph * cd8lym0 - enter_lymph * cd8rec0 + exit_peripheral_cd8 * cd8peripheral0 - enter_peripheral * cd8rec0 + k_nT_pro/bloodvolume * (cd8rec0/bloodvolume)/(cd8rec0/bloodvolume + env_naiveT);

    du.transT.cd4rec1 = sum_SP41 - death_naiveT_cd4 * cd4rec1 + exit_lymph * cd4lym1 - enter_lymph * cd4rec1 + exit_peripheral_cd4 * cd4peripheral1 - enter_peripheral * cd4rec1 + k_nT_pro/bloodvolume * (cd4rec1/bloodvolume)/(cd4rec1/bloodvolume + env_naiveT);
    du.transT.cd8rec1 = sum_SP81 - death_naiveT_cd8 * cd8rec1 + exit_lymph * cd8lym1 - enter_lymph * cd8rec1 + exit_peripheral_cd8 * cd8peripheral1 - enter_peripheral * cd8rec1 + k_nT_pro/bloodvolume * (cd8rec1/bloodvolume)/(cd8rec1/bloodvolume + env_naiveT);

    # dynamics of T cells in peripheral tissues
    du.cisT.cd4peripheral0 = k_nT_pro/peripheralvolume * (cd4peripheral0/peripheralvolume)/(cd4peripheral0/peripheralvolume + env_naiveT) - death_naiveT_cd4 * cd4peripheral0 - exit_peripheral_cd4 * cd4peripheral0 + enter_peripheral * cd4rec0;
    du.cisT.cd8peripheral0 = k_nT_pro/peripheralvolume * (cd8peripheral0/peripheralvolume)/(cd8peripheral0/peripheralvolume + env_naiveT) - death_naiveT_cd8 * cd8peripheral0 - exit_peripheral_cd8 * cd8peripheral0 + enter_peripheral * cd8rec0;

    du.transT.cd4peripheral1 = k_nT_pro/peripheralvolume * (cd4peripheral1/peripheralvolume)/(cd4peripheral1/peripheralvolume + env_naiveT) - death_naiveT_cd4 * cd4peripheral1 - exit_peripheral_cd4 * cd4peripheral1 + enter_peripheral * cd4rec1;
    du.transT.cd8peripheral1 = k_nT_pro/peripheralvolume * (cd8peripheral1/peripheralvolume)/(cd8peripheral1/peripheralvolume + env_naiveT) - death_naiveT_cd8 * cd8peripheral1 - exit_peripheral_cd8 * cd8peripheral1 + enter_peripheral * cd8rec1;

    # dynamics of T cells in secondary lymphoid system
    du.cisT.cd4lym0 = enter_lymph * cd4rec0 - exit_lymph * cd4lym0;
    du.cisT.cd8lym0 = enter_lymph * cd8rec0 - exit_lymph * cd8lym0; 

    du.transT.cd4lym1 = enter_lymph * cd4rec1 - exit_lymph * cd4lym1;
    du.transT.cd8lym1 = enter_lymph * cd8rec1 - exit_lymph * cd8lym1;
end


p1 = ComponentArray(
tauLT = 100.0, tauST = 20.0, tauMPP = 2.0, tauCMP = 4.0, tauBFUE = 7.0, tauCFUE = 7.0, tauRET = 3.0, tauRBCendo = 120.0, tauRBCtrans = 120.0, 
ssLT = 1275.0, aST = 1000.0, aMPP = 1000.0, aCMP = 36.0, aBFUE = 32.0, rLT = (1/(2.5*7)), 
VRET = 0.09e-12 , VRBC = 0.09e-12, MWHh = 64500.0, 
ksynalpha = 6e-7, ratiosynbeta = 0.5, ratiosyngamma = 0.03, ratiosyndelta = 0.04, thalfmonomer = 0.25, kon = 1.0e-5*(60*60*24), 
Kdalphabeta = 1.0e-3, Kdalphagamma = 1.0e-5, Kdalphadelta = 1.0e-2, KdHbA = 100.0, KdHbF = 100.0, KdHbA2 = 100.0, 
HbA_saturation = 0.74, HbF_saturation = 0.88, HbA2_saturation = 0.0, bloodvolume = 5.0, 
kCMP2GMP=2.5, betaGM=3.0, aGMP=128.0, tauGMP=0.12, 
alphaMPP2CLP = 0.022, betaCLP = 3.0, alphaCLP2proB = 3.0, kCLP = 0.015, CLP2proB_amp = 4.0, K = 8.4e9, gamma = 2.2, delta_oe = 2.0, delta_r = 0.2, mu_i = 0.4, delta_i_t = 2.4, delta_i_re = 0.76, phi_s = 0.12, mu_re = 0.032, phi_BM = 3.76, mu_t = 0.12, delta_t = 0.12, epsilon_spl = 0.032, epsilon_spl_r = 0.032, 
alphaCLP2DN = 2.5e-4, CLP2DN_amp = 512.0, pN = 0.23, pP = 4.5, pS = 0.23, alpha4 = 0.06, alpha8 = 0.01, alpha_muN = 0.29, alpha_muP = 0.2029, alpha_e = 0.994, alpha_r = 0.48, 
delta_dn = 0.0, delta_dp = 1.0e-6, delta_dp_r = 1.0e-6, delta_sp = 0.0, muLP = 0.37, death_naiveT_cd8 = 0.002, death_naiveT_cd4 = 0.002, 
enter_lymph = 0.07, enter_peripheral = 4.8, exit_lymph = 1.13, exit_peripheral_cd4= 0.55, exit_peripheral_cd8= 1.0, env_naiveT = 1.0e9, k_nT_pro = 3.2e8 , peripheralvolume = 60.0
);

init_cis_erythrocyte = ComponentArray(LT0=1000.0, ST01=0.0, ST02=0.0, MPP01=0.0, MPP02=0.0, CMP01=0.0, CMP02=0.0, BFUE01=0.0, BFUE02=0.0, CFUE01=0.0, CFUE02=0.0, RET0=0.0, RBC0=0.0);
init_trans_erythrocyte = ComponentArray(LT1=100.0, ST11=0.0, ST12=0.0, MPP11=0.0, MPP12=0.0, CMP11=0.0, CMP12=0.0, BFUE11=0.0, BFUE12=0.0, CFUE11=0.0, CFUE12=0.0, RET1=0.0, RBC1=0.0);
init_cis_hb_ret = ComponentArray(alpha_RET0=0.0, beta_RET0=0.0, gamma_RET0=0.0, delta_RET0=0.0, alphabeta_RET0=0.0, alphagamma_RET0=0.0, alphadelta_RET0=0.0, HbA_RET0=0.0, HbF_RET0=0.0, HbA2_RET0=0.0);
init_cis_hb_rbc = ComponentArray(alpha_RBC0=0.0, beta_RBC0=0.0, gamma_RBC0=0.0, delta_RBC0=0.0, alphabeta_RBC0=0.0, alphagamma_RBC0=0.0, alphadelta_RBC0=0.0, HbA_RBC0=0.0, HbF_RBC0=0.0, HbA2_RBC0=0.0);
init_trans_hb_ret = ComponentArray(alpha_RET1=0.0, beta_RET1=0.0, gamma_RET1=0.0, delta_RET1=0.0, alphabeta_RET1=0.0, alphagamma_RET1=0.0, alphadelta_RET1=0.0, HbA_RET1=0.0, HbF_RET1=0.0, HbA2_RET1=0.0);
init_trans_hb_rbc = ComponentArray(alpha_RBC1=0.0, beta_RBC1=0.0, gamma_RBC1=0.0, delta_RBC1=0.0, alphabeta_RBC1=0.0, alphagamma_RBC1=0.0, alphadelta_RBC1=0.0, HbA_RBC1=0.0, HbF_RBC1=0.0, HbA2_RBC1=0.0);
init_myeloid = ComponentArray(GMP01=0.0, GMP02=0.0, GMP11=0.0, GMP12=0.0, GM0=0.0, GM1=0.0); 
init_cisB = ComponentArray(CLP0=0.0, Boe0=0.0, Bi0=0.0, BMrec0=0.0, Bt0=0.0, BMspl0=0.0);
init_transB = ComponentArray(CLP1=0.0, Boe1=0.0, Bi1=0.0, BMrec1=0.0, Bt1=0.0, BMspl1=0.0);
init_cisT = ComponentArray(N00=0.0, N10=0.0, N20=0.0, N30=0.0, N40=0.0, P00=0.0, P10=0.0, P20=0.0, P30=0.0, P40=0.0, P50=0.0, P60=0.0, P70=0.0, S400=0.0, S410=0.0, S420=0.0, S800=0.0, S810=0.0, S820=0.0, 
                           cd4rec0=0.0, cd8rec0=0.0, cd4peripheral0=0.0, cd8peripheral0=0.0, cd4lym0=0.0, cd8lym0 =0.0);
init_transT = ComponentArray(N01=0.0, N11=0.0, N21=0.0, N31=0.0, N41=0.0, P01=0.0, P11=0.0, P21=0.0, P31=0.0, P41=0.0, P51=0.0, P61=0.0, P71=0.0, S401=0.0, S411=0.0, S421=0.0, S801=0.0, S811=0.0, S821=0.0, 
                             cd4rec1=0.0, cd8rec1=0.0, cd4peripheral1=0.0, cd8peripheral1=0.0, cd4lym1=0.0, cd8lym1 =0.0);

u0 = ComponentArray(cis_erythrocyte = init_cis_erythrocyte, trans_erythrocyte = init_trans_erythrocyte, 
    cis_hb_ret = init_cis_hb_ret, cis_hb_rbc = init_cis_hb_rbc, trans_hb_ret = init_trans_hb_ret, trans_hb_rbc = init_trans_hb_rbc, 
    myeloid = init_myeloid, cisB = init_cisB, transB = init_transB, cisT = init_cisT, transT = init_transT);

prob = ODEProblem(HSChomo!, u0, (0.0, 600.0), p1);
sol = solve(prob, saveat=1.0, reltol = 1e-8);

LT0 = [sol.u[i].cis_erythrocyte.LT0 for i in 1:length(sol.t)]
LT1 = [sol.u[i].trans_erythrocyte.LT1 for i in 1:length(sol.t)]

RBC0 = [sol.u[i].cis_erythrocyte.RBC0 for i in 1:length(sol.t)]
RBC1 = [sol.u[i].trans_erythrocyte.RBC1 for i in 1:length(sol.t)]

GM0 = [sol.u[i].myeloid.GM0 for i in 1:length(sol.t)]
GM1 = [sol.u[i].myeloid.GM1 for i in 1:length(sol.t)]

BMrec0 = [sol.u[i].cisB.BMrec0 for i in 1:length(sol.t)]
BMrec1 = [sol.u[i].transB.BMrec1 for i in 1:length(sol.t)]

cd4rec0 = [sol.u[i].cisT.cd4rec0 for i in 1:length(sol.t)]
cd8rec0 = [sol.u[i].cisT.cd8rec0 for i in 1:length(sol.t)]

cd4rec1 = [sol.u[i].transT.cd4rec1 for i in 1:length(sol.t)]
cd8rec1 = [sol.u[i].transT.cd8rec1 for i in 1:length(sol.t)]

obs = CSV.read("rbc.csv",DataFrame);
@rsubset!(obs, :time in [1, 10, 30, 50, 70, 100, 200, 300, 400, 500, 600]);

pRBC = plot(title = "RBC", palette=:PiYG_4, legendtitle = "Source, solver", legend = :outerright);
plot!(sol.t, RBC0, alpha = 0.9, label = "endogenous, julia", linewidth = 2);
plot!(sol.t, RBC1, alpha = 0.9, label = "transduced, julia", linewidth = 2);
scatter!(obs.time, obs.RBC0, alpha = 0.9, label = "endogenous, mrgsolve");
scatter!(obs.time, obs.RBC1, alpha = 0.9, label = "transduced, mrgsolve");
xlabel!("time (day)");

pGM = plot(title = "Granulocytes", palette=:PiYG_4, legendtitle = "Source, solver", legend = :outerright);
plot!(sol.t, GM0, alpha = 0.9, label = "endogenous, julia", linewidth = 2);
plot!(sol.t, GM1, alpha = 0.9, label = "transduced, julia", linewidth = 2);
scatter!(obs.time, obs.GM0, alpha = 0.9, label = "endogenous, mrgsolve");
scatter!(obs.time, obs.GM1, alpha = 0.9, label = "transduced, mrgsolve");
xlabel!("time (day)");

pB = plot(title = "circulating B", palette=:PiYG_4, legendtitle = "Source, solver", legend = :outerright);
plot!(sol.t, BMrec0, alpha = 0.9, label = "endogenous, julia", linewidth = 2);
plot!(sol.t, BMrec1, alpha = 0.9, label = "transduced, julia", linewidth = 2);
scatter!(obs.time, obs.BMrec0, alpha = 0.9, label = "endogenous, mrgsolve");
scatter!(obs.time, obs.BMrec1, alpha = 0.9, label = "transduced, mrgsolve");
xlabel!("time (day)");

pT = plot(title = "circulating T", palette=:PiYG_8, legendtitle = "Source, solver", legend = :outerright);
plot!(sol.t, cd4rec0, alpha = 0.9, label = "endogenous CD4+, julia", linewidth = 2);
plot!(sol.t, cd4rec1, alpha = 0.9, label = "transduced CD4+, julia", linewidth = 2);
plot!(sol.t, cd8rec0, alpha = 0.9, label = "endougenous CD8+, julia", linewidth = 2);
plot!(sol.t, cd8rec1, alpha = 0.9, label = "transduced CD8+, julia", linewidth = 2);
scatter!(obs.time, obs.cd4rec0, alpha = 0.9, label = "endogenous CD4+, mrgsolve");
scatter!(obs.time, obs.cd4rec1, alpha = 0.9, label = "transduced CD4+, mrgsolve");
scatter!(obs.time, obs.cd8rec0, alpha = 0.9, label = "endougenous CD8+, mrgsolve");
scatter!(obs.time, obs.cd8rec1, alpha = 0.9, label = "transduced CD8+, mrgsolve");
xlabel!("time (day)");

plot(pRBC, pGM, pB, pT, layout=(2,2), size = (1200,800))
savefig("JuliaVerification.png")
