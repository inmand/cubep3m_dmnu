#!/bin/bash
#SBATCH --job-name=cubep3m_npbh
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=8GB

export OMP_STACKSIZE="100M"
export OMP_NUM_THREADS="4"

source /home/dbi208/.bashrc
iload

ORGDIR=$HOME/src/cubep3m_dmnu
RUNDIR=$SCRATCH/PBH/Test6

cp -r $ORGDIR $RUNDIR/source
SUBDIR=$RUNDIR/source/workdir/batch

echo ${SLURM_JOB_NODELIST} > $RUNDIR/list_of_nodes
cat /proc/cpuinfo | grep "model name" | head -n1 >> $RUNDIR/list_of_nodes

LOG=_log

cd /home/dbi208/src/p3dfft.2.5.1_fix
source dobuild_nyu >& ${RUNDIR}/p3dfft${LOG}

cd ${SUBDIR}
cd ..
ln -s parameters_npbh parameters
cd  ./batch/
source ./COMPILE_npbh.csh >& ${RUNDIR}/compile${LOG}

mpirun ../utils/npbh/ic_cdm >& ${RUNDIR}/ic_cdm${LOG}
mpirun ../utils/npbh/ic_pbh >& ${RUNDIR}/ic_pbh${LOG}
mpirun ../source_threads/cubep3m >& ${RUNDIR}/cubep3m${LOG}
mpirun ../utils/npbh/ngp_power >& ${RUNDIR}/ngp_power${LOG}
mpirun ../utils/npbh/ngp_veldivg >& ${RUNDIR}/ngp_veldivg${LOG}

