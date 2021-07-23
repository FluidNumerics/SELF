#!/usr/bin/python


import json
from datetime import datetime
import os

WORKSPACE=os.getenv('WORKSPACE')
GIT_SHA=os.getenv('GIT_SHA')
BUILD_ID=os.getenv('BUILD_ID')
PLATFORM=os.getenv('PLATFORM')
NODE_COUNT=os.getenv('NODE_COUNT')
MACHINE_TYPE=os.getenv('MACHINE_TYPE')
GPU_COUNT=os.getenv('GPU_COUNT')
GPU_TYPE=os.getenv('GPU_TYPE')
BUILD_CONTAINER_PLATFORM=os.getenv('BUILD_CONTAINER_PLATFORM')
RUN_CONTAINER_PLATFORM=os.getenv('RUN_CONTAINER_PLATFORM')
BUILD_TYPE=os.getenv('BUILD_TYPE')
TARGET_PLATFORM=os.getenv('TARGET_PLATFORM')
COMPILER=os.getenv('COMPILER')
MPI_PROVIDER=os.getenv('MPI_PROVIDER')



def main():

    with open(WORKSPACE+'/tests.json','r')as f:          
      tests = json.load(f)

    utc = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

    results = []
    for test in tests :
        results.append({"cli_command": test['cli_command'],
                        "execution_command": test['execution_command'],
                        "exit_code": test['exit_code'],
                        "stdout": test['stdout'],
                        "stderr": test['stderr'],
                        "build_id": BUILD_ID,
                        "build_type": BUILD_TYPE,
                        "compiler": COMPILER,
                        "container_platform": BUILD_CONTAINER_PLATFORM,
                        "container_platform_runtime": RUN_CONTAINER_PLATFORM,
                        "mpi_provider": MPI_PROVIDER,
                        "target_platform": TARGET_PLATFORM,
                        "machine_type": MACHINE_TYPE,
                        "gpu_type":GPU_TYPE,
                        "gpu_count":int(GPU_COUNT),
                        "node_count":int(NODE_COUNT),
                        "git_sha": GIT_SHA,
                        "datetime": utc})

    with open(WORKSPACE+'/flat_results.json','w')as f:          
        for res in results:
            f.write(json.dumps(res))
            f.write('\n')
      
if __name__ == '__main__':
    main()
