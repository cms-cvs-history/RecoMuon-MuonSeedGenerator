#! /bin/sh

eval `scramv1 runtime -sh`
mkdir h_mix210pre6
echo 'block ptSeedParameterization = {' > h_mix210pre6/ptSeedParameterization.cfi
echo '# DT Parameters ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# Det_Stations          p1        p2        p3       ep1      ep2      ep3' >> h_mix210pre6/ptSeedParameterization.cfi
echo '#------------------------------------------------------------------------#'>> h_mix210pre6/ptSeedParameterization.cfi
echo '.x SeedParaFit.C(3,12,2,100,3)' | root.exe -b
echo '.x SeedParaFit.C(3,13,2,90,3)' | root.exe -b
echo '.x SeedParaFit.C(3,14,2,80,3)' | root.exe -b
echo '.x SeedParaFit.C(3,23,2,90,3)' | root.exe -b
echo '.x SeedParaFit.C(3,24,2,80,3)' | root.exe -b
echo '.x SeedParaFit.C(3,34,2,80,3)' | root.exe -b

echo '  ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# CSC Parameters ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# Det_Stations      p1        p2        p3       ep1      ep2      ep3' >> h_mix210pre6/ptSeedParameterization.cfi
echo '#------------------------------------------------------------------------#'>> h_mix210pre6/ptSeedParameterization.cfi
echo '.x SeedParaFit.C(1,11,55,68,2)' | root.exe -b
echo '.x SeedParaFit.C(1,12,2,60,3)' | root.exe -b
echo '.x SeedParaFit.C(1,12,65,140,3)' | root.exe -b
echo '.x SeedParaFit.C(1,13,12,52,3)'  | root.exe -b
echo '.x SeedParaFit.C(1,13,70,140,3)' | root.exe -b
echo '.x SeedParaFit.C(1,14,80,140,3)' | root.exe -b
echo '.x SeedParaFit.C(1,23,2,140,3)' | root.exe -b
echo '.x SeedParaFit.C(1,24,80,140,3)' | root.exe -b
echo '.x SeedParaFit.C(1,34,80,140,3)' | root.exe -b

echo '  ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# OL Parameters ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# Det_Layers        p1        p2        p3       ep1      ep2      ep3' >> h_mix210pre6/ptSeedParameterization.cfi
echo '#------------------------------------------------------------------------#'>> h_mix210pre6/ptSeedParameterization.cfi
echo '.x SeedParaFit.C(2,1213,5,35,2)' | root.exe -b
echo '.x SeedParaFit.C(2,1222,18,27,2)' | root.exe -b
echo '.x SeedParaFit.C(2,1232,25,40,2)' | root.exe -b
echo '.x SeedParaFit.C(2,2213,5,35,2)' | root.exe -b
echo '.x SeedParaFit.C(2,2222,18,27,2)' | root.exe -b

echo '# Sigle CSC Segments ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# CSCDet_Id         p1        p2        p3       ep1      ep2      ep3' >> h_mix210pre6/ptSeedParameterization.cfi
echo '#------------------------------------------------------------------------#'>> h_mix210pre6/ptSeedParameterization.cfi
echo '.x SeedParaFit.C(5,11,62,138,2)' | root.exe -b
echo '.x SeedParaFit.C(5,12,22,60,2)' | root.exe -b
echo '.x SeedParaFit.C(5,13,1,10,2)' | root.exe -b
echo '.x SeedParaFit.C(5,21,65,138,2)' | root.exe -b
echo '.x SeedParaFit.C(5,22,2,52,2)' | root.exe -b
echo '.x SeedParaFit.C(5,31,78,138,2)' | root.exe -b
echo '.x SeedParaFit.C(5,32,12,65,2)' | root.exe -b
echo '.x SeedParaFit.C(5,41,85,138,2)' | root.exe -b

echo '# Sigle DT Segments ' >> h_mix210pre6/ptSeedParameterization.cfi
echo '# DTDet_Id          p1        p2        p3       ep1      ep2      ep3' >> h_mix210pre6/ptSeedParameterization.cfi
echo '#------------------------------------------------------------------------#'>> h_mix210pre6/ptSeedParameterization.cfi
echo '.x SeedParaFit.C(4,10,1,25,2)' | root.exe -b
echo '.x SeedParaFit.C(4,11,35,75,2)' | root.exe -b
echo '.x SeedParaFit.C(4,12,88,110,2)' | root.exe -b
echo '.x SeedParaFit.C(4,20,1,21,2)' | root.exe -b
echo '.x SeedParaFit.C(4,21,30,66,2)' | root.exe -b
echo '.x SeedParaFit.C(4,22,75,102,2)' | root.exe -b
echo '.x SeedParaFit.C(4,30,1,17,2)' | root.exe -b
echo '.x SeedParaFit.C(4,31,25,57,2)' | root.exe -b
echo '.x SeedParaFit.C(4,32,62,88,2)' | root.exe -b
echo ' '>> h_mix210pre6/ptSeedParameterization.cfi
echo '} '>> h_mix210pre6/ptSeedParameterization.cfi

