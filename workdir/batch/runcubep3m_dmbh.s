#!/bin/bash
#SBATCH --job-name=cubep3m_dmbh
#SBATCH --nodes=8
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=16GB

export OMP_STACKSIZE="100M"
export OMP_NUM_THREADS="8"

source /home/dbi208/.bashrc
iload

cd $SLURM_SUBMIT_DIR
SUBDIR=$SLURM_SUBMIT_DIR
RUNDIR=$SUBDIR/../../

echo ${SLURM_JOB_NODELIST} > $RUNDIR/list_of_nodes
cat /proc/cpuinfo | grep "model name" | head -n1 >> $RUNDIR/list_of_nodes

LOG=_log

#cd /home/dbi208/src/p3dfft.2.5.1_fix
#source dobuild_nyu >& ${RUNDIR}/p3dfft${LOG}

cd ${SUBDIR}
cd ..
ln -s parameters_dmbh parameters
cd  ./batch/
source ./COMPILE_dmbh.csh >& ${RUNDIR}/compile${LOG}

srun ../utils/dist_init_dmbh/dist_init_dmbh >& ${RUNDIR}/dist_init_dmbh${LOG}
srun ../source_threads/cubep3m >& ${RUNDIR}/cubep3m${LOG}
srun ../utils/dist_init_dmbh/ngp_power >& ${RUNDIR}/ngp_power${LOG}
#srun ../utils/dist_init_dmbh/ngp_veldivg >& ${RUNDIR}/ngp_veldivg${LOG}

