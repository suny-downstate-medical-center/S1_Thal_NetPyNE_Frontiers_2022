#!/bin/bash
#SBATCH --nodes=16            # node
#SBATCH --ntasks-per-node=48   # tasks per node
#SBATCH --time=24:00:00               # time limit
#SBATCH --partition=g100_usr_bmem
#SBATCH --account=icei_H_King

srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py
srun ./x86_64/special -mpi -python init.py

