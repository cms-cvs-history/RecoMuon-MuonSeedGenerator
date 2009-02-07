import FWCore.ParameterSet.Config as cms

dphiScale = cms.PSet(
    # =============================================================== 
    #              idx              p0       p1      p2      r1      r2 
    # --------------------------------------------------------------- 
    CSC_01_1_scale = cms.vdouble(0.409773, -1.35896, 0.0, -3.316374, 0.0),
    CSC_12_1_scale = cms.vdouble(0.082276, -7.24509, 0.0, -88.058295, 0.0),
    CSC_12_2_scale = cms.vdouble(0.211236, -1.838651, 0.0, -8.704266, 0.0),
    CSC_12_3_scale = cms.vdouble(0.773482, -1.050894, 0.0, -1.358653, 0.0),
    CSC_13_2_scale = cms.vdouble(0.26387, -3.923963, 0.0, -14.87081, 0.0),
    CSC_13_3_scale = cms.vdouble(0.90177, -1.75464, 0.0, -1.945774, 0.0),
    CSC_14_3_scale = cms.vdouble(0.957255, -2.821265, 0.0, -2.947245, 0.0),
    CSC_23_1_scale = cms.vdouble(0.063805, -7.109899, 0.0, -111.431855, 0.0),
    CSC_23_2_scale = cms.vdouble(0.114126, -2.568468, 0.0, -22.505524, 0.0),
    CSC_24_1_scale = cms.vdouble(0.216754, -0.652681, 0.0, -3.011161, 0.0),
    CSC_34_1_scale = cms.vdouble(0.060732, -2.801265, 0.0, -46.125126, 0.0),
    DT_12_1_scale = cms.vdouble(0.212302, -2.194977, 0.0, -10.338952, 0.0),
    DT_12_2_scale = cms.vdouble(0.181672, -1.956069, 0.0, -10.767036, 0.0),
    DT_13_1_scale = cms.vdouble(0.364866, -3.123576, 0.0, -8.560875, 0.0),
    DT_13_2_scale = cms.vdouble(0.341619, -2.537996, 0.0, -7.429325, 0.0),
    DT_14_1_scale = cms.vdouble(0.417764, -4.474353, 0.0, -10.71025, 0.0),
    DT_14_2_scale = cms.vdouble(0.395697, -3.801596, 0.0, -9.607341, 0.0),
    DT_23_1_scale = cms.vdouble(0.150846, -3.325731, 0.0, -22.047226, 0.0),
    DT_23_2_scale = cms.vdouble(0.130436, -3.082495, 0.0, -23.632156, 0.0),
    DT_24_1_scale = cms.vdouble(0.208306, -6.083214, 0.0, -29.203249, 0.0),
    DT_24_2_scale = cms.vdouble(0.193048, -5.179807, 0.0, -26.831762, 0.0),
    DT_34_1_scale = cms.vdouble(0.053875, -14.200023, 0.0, -263.573067, 0.0),
    DT_34_2_scale = cms.vdouble(0.048121, -10.68666, 0.0, -222.080722, 0.0),
    OL_1213_0_scale = cms.vdouble(0.305767, -2.727901, 0.0, -8.921499, 0.0),
    OL_1222_0_scale = cms.vdouble(0.251284, -4.058942, 0.0, -16.152822, 0.0),
    OL_1232_0_scale = cms.vdouble(0.208174, -4.753092, 0.0, -22.832256, 0.0),
    OL_2213_0_scale = cms.vdouble(0.131266, -4.700717, 0.0, -35.810724, 0.0),
    OL_2222_0_scale = cms.vdouble(0.125065, -5.775142, 0.0, -46.177302, 0.0),
    SMB_10_0_scale = cms.vdouble(1.412012, 2.608279, 0.0, 1.847207, 0.0),
    SMB_11_0_scale = cms.vdouble(1.422246, 2.746887, 0.0, 1.931372, 0.0),
    SMB_12_0_scale = cms.vdouble(1.285671, 2.252121, 0.0, 1.751709, 0.0),
    SMB_20_0_scale = cms.vdouble(1.044853, 1.367932, 0.0, 1.30921, 0.0),
    SMB_21_0_scale = cms.vdouble(1.01585, 1.254693, 0.0, 1.235117, 0.0),
    SMB_22_0_scale = cms.vdouble(0.902138, 1.399135, 0.0, 1.550911, 0.0),
    SMB_30_0_scale = cms.vdouble(0.578687, -3.820392, 0.0, -6.601832, 0.0),
    SMB_31_0_scale = cms.vdouble(0.544333, -3.404327, 0.0, -6.254127, 0.0),
    SMB_32_0_scale = cms.vdouble(0.501362, -3.683311, 0.0, -7.346607, 0.0),
    SME_11_0_scale = cms.vdouble(1.055973, 2.059006, 0.0, 1.949866, 0.0),
    SME_12_0_scale = cms.vdouble(1.001529, 1.820654, 0.0, 1.817874, 0.0),
    SME_13_0_scale = cms.vdouble(0.623258, -2.74588, 0.0, -4.40569, 0.0),
    SME_21_0_scale = cms.vdouble(0.49473, -1.3401, 0.0, -2.708752, 0.0),
    SME_22_0_scale = cms.vdouble(0.414984, -0.424861, 0.0, -1.023801, 0.0)
)


