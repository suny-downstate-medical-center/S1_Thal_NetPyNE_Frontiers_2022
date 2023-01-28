#!/bin/bash
#SBATCH --nodes=35                    # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=24:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_bprod
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_0_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_1_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_2_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_3_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_4_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_5_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_6_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_7_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_8_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_9_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v103_batch2/v103_batch2_0_cfg.json netParams=../data/v103_batch2/v103_batch2_netParams.py

