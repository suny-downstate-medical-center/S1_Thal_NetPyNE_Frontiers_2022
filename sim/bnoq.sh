#!/bin/bash
#SBATCH --nodes=35                    # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=12:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_bprod
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_0_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_1_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_2_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_3_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_4_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v104_batch3/v104_batch3_0_cfg.json netParams=../data/v104_batch3/v104_batch3_netParams.py

