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

    flatFile = []
    for test in ci_conf['tests'] :
      for nel in test['nelems'] :
        for nvar in test['nvars'] :
          for cQuad in test['control_quadrature'] :
            for cDeg in test['control_degree'] :
              for tQuad in test['target_quadrature'] :
                for tDeg in test['target_degree'] :
                  for addlOpts in test['additional_opts'] :
                    for funcOpts in test['function_opts'] :

                      flatFile.append({"cli_command": test['cli_command'],
                                       "control_quadrature" : cQuad,
                                       "control_degree" : cDeg,
                                       "nelems": nel,
                                       "nvars": nvar,
                                       "target_quadrature" : tQuad, 
                                       "target_degree" : tDeg, 
                                       "additional_opts" : addlOpts,
                                       "function_opts" : funcOpts,
                                       "git_sha": GIT_SHA,
                                       "build_id": BUILD_ID,
                                       "datetime": utc,
                                       "stdout": test['stdout'],
                                       "stderr": test['stderr'],
                                       "exit_code": test['exit_code']})

    with open(WORKSPACE+'/flat_restuls.json','w')as f:          
      f.write(json.dumps(tests,indent=2))
      
    print(json.dumps(tests,indent=2))      
