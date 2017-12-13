#!/bin/bash
#SBATCH --job-name=cubep3m_dmnu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem=8GB

export OMP_STACKSIZE="100M"
export OMP_NUM_THREADS="4"

source /home/dbi208/.bashrc
iload

ORGDIR=$HOME/src/cubep3m_dmnu
RUNDIR=$SCRATCH/DMNU/Test14/cubep3m_dm

cp -r $ORGDIR $RUNDIR/source
SUBDIR=$RUNDIR/source/workdir/batch

echo ${SLURM_JOB_NODELIST} > $RUNDIR/list_of_nodes
cat /proc/cpuinfo | grep "model name" | head -n1 >> $RUNDIR/list_of_nodes

LOG=_log

cd /home/dbi208/src/p3dfft.2.5.1_fix
source dobuild_nyu >& ${RUNDIR}/p3dfft${LOG}

cd ${SUBDIR}
cd ..
ln -s parameters_dm parameters
cd  ./batch/
source ./COMPILE_dmnu.csh >& ${RUNDIR}/compile${LOG}

mpirun ../utils/dist_init_dmnu/dist_init_dmnu_dm >& ${RUNDIR}/dist_init_dmnu_dm${LOG}
mpirun ../source_threads/cubep3m >& ${RUNDIR}/cubep3m${LOG}
mpirun ../utils/cic_power/ngp_power_dm >& ${RUNDIR}/ngp_power_dm${LOG}
mpirun ../utils/cic_power/ngp_power_dm_init >& ${RUNDIR}/ngp_power_dm_init${LOG}
#mpirun ../utils/cic_velpower/ngp_veldivg >& ${RUNDIR}/ngp_veldivg${LOG}

