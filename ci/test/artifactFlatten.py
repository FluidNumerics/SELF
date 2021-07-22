#!/usr/bin/python


import json
from datetime import datetime
import os

WORKSPACE=os.getenv('WORKSPACE')
GIT_SHA=os.getenv('SHORT_SHA')
BUILD_ID=os.getenv('BUILD_ID')


def main():

    with open(WORKSPACE+'/tests.json','r')as f:          
      tests = json.load(f)

    utc = datetime.utcnow().strftime("%Y-%m-%dT-%H:%M:%S")

    results = []
    for test in tests :
        benchmark_info = test['benchmark_info']
        funcOpts = benchmark_info['function_opts']
        addlOpts = benchmark_info['additional_opts']
        results.append({"benchmark_info.cli_command": benchmark_info['cli_command'],
                         "benchmark_info.control_quadrature" : benchmark_info['control_quadrature'],
                         "benchmark_info.control_degree" : benchmark_info['control_degree'],
                         "benchmark_info.nelements": benchmark_info['nelements'],
                         "benchmark_info.nvar": benchmark_info['nvar'],
                         "benchmark_info.target_quadrature" : benchmark_info['target_quadrature'], 
                         "benchmark_info.target_degree" : benchmark_info['target_degree'], 
                         "benchmark_info.additional_opts.name" : addlOpts['name'],
                         "benchmark_info.additional_opts.value" : addlOpts['value'],
                         "benchmark_info.function_opts.name" : funcOpts['name'],
                         "benchmark_info.function_opts.value" : funcOpts['value'],
                         "git_sha": GIT_SHA,
                         "build_id": BUILD_ID,
                         "datetime": utc,
                         "stdout": test['stdout'],
                         "stderr": test['stderr'],
                         "exit_code": test['exit_code']})

    with open(WORKSPACE+'/flat_results.json','w')as f:          
        for res in results:
            f.write(json.dumps(res))
            f.write('\n')
      
if __name__ == '__main__':
    main()
