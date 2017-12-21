#!/bin/bash
#SBATCH --job-name=cubep3m_dmbh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=4GB

export OMP_STACKSIZE="100M"
export OMP_NUM_THREADS="8"

source /home/dbi208/.bashrc
iload

ORGDIR=$HOME/src/cubep3m_dmnu
RUNDIR=$SCRATCH/PBH/Test11_0_0

cp -r $ORGDIR $RUNDIR/source
SUBDIR=$RUNDIR/source/workdir/batch

echo ${SLURM_JOB_NODELIST} > $RUNDIR/list_of_nodes
cat /proc/cpuinfo | grep "model name" | head -n1 >> $RUNDIR/list_of_nodes

LOG=_log

cd /home/dbi208/src/p3dfft.2.5.1_fix
source dobuild_nyu >& ${RUNDIR}/p3dfft${LOG}

cd ${SUBDIR}
cd ..
ln -s parameters_dmbh parameters
cd  ./batch/
source ./COMPILE_dmbh.csh >& ${RUNDIR}/compile${LOG}

srun ../utils/dist_init_dmbh/dist_init_dmbh >& ${RUNDIR}/dist_init_dmbh${LOG}
#srun ../utils/dist_init_dmbh/ic_pbh >& ${RUNDIR}/ic_pbh${LOG}
srun ../source_threads/cubep3m >& ${RUNDIR}/cubep3m${LOG}
srun ../utils/dist_init_dmbh/ngp_power >& ${RUNDIR}/ngp_power${LOG}
srun ../utils/dist_init_dmbh/ngp_veldivg >& ${RUNDIR}/ngp_veldivg${LOG}

