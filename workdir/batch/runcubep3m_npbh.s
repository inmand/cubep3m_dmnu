#!/bin/bash
#SBATCH --job-name=cubep3m_npbh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=8GB

source /home/dbi208/.bashrc
iload

SUBDIR=$HOME/src/cubep3m_pbh/workdir/batch
RUNDIR=$SCRATCH/PBH/TestPBH_0

echo ${SLURM_JOB_NODELIST} > $RUNDIR/nodelist

export OMP_STACKSIZE="100M"
export OMP_NUM_THREADS="4"

cd /home/dbi208/src/p3dfft.2.5.1_fix
source dobuild_nyu >& ${RUNDIR}/p3dfft_log

cd ${SUBDIR}
cp ../parameters ${RUNDIR}
cp ./run_cpbh.s ${RUNDIR}

source ./COMPILE_CPBH.csh >& ${RUNDIR}/compile_log

mpirun ../utils/dist_init/dist_init_cpbh >& ${RUNDIR}/dist_init_dmu_dm_log
mpirun ../source_threads/cubep3m_dm >& ${RUNDIR}/cubep3m_dm_log
mpirun ../utils/cic_power/ngp_power_dm_init >& ${RUNDIR}/ngp_power_dm_log
