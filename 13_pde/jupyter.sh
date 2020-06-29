#!/bin/sh
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=24:00:00
#$ -o log/$JOB_ID.out
#$ -j y
#$ -N linsho

. /etc/profile.d/modules.sh
module load cuda/10.0.130
module load cudnn/7.4 
module load python/3.8.3
source /home/2/15B03631/dev/hpc_lecture/13_pde/venv/bin/activate
/home/2/15B03631/dev/kaggle/venv/bin/jupyter notebook --no-browser --ip=0.0.0.0 --port=8888
