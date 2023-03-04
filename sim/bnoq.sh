#!/bin/bash
#SBATCH --nodes=35                    # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=22:00:00               # time limits: 1 hour
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_bprod
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_0_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_0_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_1_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v108_batch1/v108_batch1_0_cfg.json netParams=../data/v108_batch1/v108_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v108_batch1/v108_batch1_1_cfg.json netParams=../data/v108_batch1/v108_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v108_batch1/v108_batch1_2_cfg.json netParams=../data/v108_batch1/v108_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v108_batch1/v108_batch1_3_cfg.json netParams=../data/v108_batch1/v108_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v108_batch1/v108_batch1_4_cfg.json netParams=../data/v108_batch1/v108_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_2_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_3_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py
srun ./x86_64/special -mpi -python init.py simConfig=../data/v107_batch1/v107_batch1_4_cfg.json netParams=../data/v107_batch1/v107_batch1_netParams.py

