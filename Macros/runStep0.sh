#!/bin/bash

. setupJSONs.sh
. setupPaths.sh

if [[ ! -d ${NTUPLE_DIR} ]] ; then
  mkdir -p ${NTUPLE_DIR}
fi

cd Step0Scripts

# Submit job for all MC
qsub Step0_GetPUWeights_TTbarJob.sh
qsub Step0_GetPUWeights_TTbarLOJob.sh
qsub Step0_GetPUWeights_TTLepJob.sh
qsub Step0_GetPUWeights_TTpowJob.sh
qsub Step0_GetPUWeights_Wjets_70to100Job.sh
qsub Step0_GetPUWeights_Wjets_100to200Job.sh
qsub Step0_GetPUWeights_Wjets_200to400Job.sh
qsub Step0_GetPUWeights_Wjets_400to600Job.sh
qsub Step0_GetPUWeights_Wjets_600to800Job.sh
qsub Step0_GetPUWeights_Wjets_800to1200Job.sh
qsub Step0_GetPUWeights_Wjets_1200to2500Job.sh
qsub Step0_GetPUWeights_Wjets_2500toInfJob.sh
qsub Step0_GetPUWeights_otherMCJob1.sh
qsub Step0_GetPUWeights_otherMCJob2.sh
qsub Step0_GetPUWeights_otherMCJob3.sh
qsub Step0_GetPUWeights_otherMCJob4.sh

# Submit job for all Data
qsub Step0_GetPUWeights_DataJob.sh

# Submit job for each group of signal points with a given stop mass
qsub Step0_GetPUWeights_Stop250Job.sh
qsub Step0_GetPUWeights_Stop275Job.sh
qsub Step0_GetPUWeights_Stop300Job.sh
qsub Step0_GetPUWeights_Stop325Job.sh
qsub Step0_GetPUWeights_Stop350Job.sh
qsub Step0_GetPUWeights_Stop375Job.sh
qsub Step0_GetPUWeights_Stop400Job.sh
qsub Step0_GetPUWeights_Stop425Job.sh
qsub Step0_GetPUWeights_Stop450Job.sh
qsub Step0_GetPUWeights_Stop475Job.sh
qsub Step0_GetPUWeights_Stop500Job.sh
qsub Step0_GetPUWeights_Stop525Job.sh
qsub Step0_GetPUWeights_Stop550Job.sh
qsub Step0_GetPUWeights_Stop575Job.sh
qsub Step0_GetPUWeights_Stop600Job.sh
qsub Step0_GetPUWeights_Stop625Job.sh
qsub Step0_GetPUWeights_Stop650Job.sh
qsub Step0_GetPUWeights_Stop675Job.sh
qsub Step0_GetPUWeights_Stop700Job.sh
qsub Step0_GetPUWeights_Stop725Job.sh
qsub Step0_GetPUWeights_Stop750Job.sh
qsub Step0_GetPUWeights_Stop775Job.sh
qsub Step0_GetPUWeights_Stop800Job.sh

cd -
