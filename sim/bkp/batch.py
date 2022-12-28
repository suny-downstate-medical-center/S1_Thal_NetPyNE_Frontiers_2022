"""
batch.py 

Batch simulation for S1-thalamus model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom():
    params = specs.ODict()
    
    # params[('seeds', 'conn')] =  [1234]

    params[('rateStimI')] = [9.0]
    params[('rateStimE')] = [9.0]

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py')

    return b

# ----------------------------------------------------------------------------------------------
# Inhibitory connections
# ----------------------------------------------------------------------------------------------
def inhib():
    params = specs.ODict()
    
    params[('IEGain')] = [0.5, 0.75, 1.0, 1.25, 1.5]
    params[('IIGain')] = [0.5, 0.75, 1.0, 1.25, 1.5]
    params[('seeds', 'conn')] =  [0, 1, 2, 3, 4]
    params[('seeds', 'conn')] =  [0, 1, 2, 3, 4]

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py')

    return b
    
# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin' or type=='mpi':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 64,
            'script': 'init.py',
            'mpiCommand': 'mpiexec', # --use-hwthread-cpus
            'skip': True}

    elif type=='mpi_direct2':
        b.runCfg = {'type': 'mpi_direct',
            'mpiCommand': 'mpirun -n 80 ./x86_64/special -mpi -python init.py', 
            'skip': True}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default',
            'walltime': '72:00:00', 
            'nodes': 1,
            'coresPerNode': 40,
            'email': 'fernandodasilvaborges@gmail.com',
            'folder': '/home/ext_fernandodasilvaborges_gmail_/S1_netpyne/sim/', 
            'script': 'init.py', 
            'mpiCommand': 'mpirun',
            'skipCustom': '_raster_gid.png'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    b = custom() #

    b.batchLabel = 'v7_batch1'  
    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'mpi_direct2')
    b.run() # run batch
