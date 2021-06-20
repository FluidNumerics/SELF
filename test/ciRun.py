#!/usr/bin/python3

import json
import os
import subprocess
import shutil

WORKSPACE=os.getenv('WORKSPACE')
INSTALL_ROOT=os.getenv('INSTALL_ROOT')
GPU_ACCEL=os.getenv('GPU_ACCEL')


def main():

    with open(INSTALL_ROOT+'/test/ci.json','r') as f:
      ci_conf = json.load(f)

    with open(INSTALL_ROOT+'/tests.json','r')as f:          
      tests = json.load(f)

    k = 0
    # Create commands to test with
    for test in ci_conf['tests'] :
      for nel in test['nelems'] :
        for nvar in test['nvars'] :
          for cQuad in test['control_quadrature'] :
            for cDeg in test['control_degree'] :
              for tQuad in test['target_quadrature'] :
                for tDeg in test['target_degree'] :
                  for addlOpts in test['additional_opts'] :
                    for funcOpts in test['function_opts'] :

                      workdir = INSTALL_ROOT+'/test/'
                      workdir += test['cli_command']+'/'
                      workdir += 'nel_{}'.format(nel)+'/'
                      workdir += 'nvar_{}'.format(nvar)+'/'
                      workdir += 'cQuad_{}'.format(cQuad)+'/'
                      workdir += 'cDeg_{}'.format(cDeg)+'/'
                      workdir += 'tQuad_{}'.format(tQuad)+'/'
                      workdir += 'tDeg_{}'.format(tDeg)+'/'
                      workdir += '{}'.format(addlOpts['name'])+'/'
                      workdir += '{}'.format(funcOpts['name'])+'/'

                      print('Running {}'.format(workdir+'/test.sh \n'))
                      proc = subprocess.Popen([workdir+'/test.sh'],
                                               shell=True,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
                      stdout, stderr = proc.communicate()
                      print(stdout.decode("utf-8"))
                      print(stderr.decode("utf-8"))
                      tests[k]['stdout'] = stdout.decode("utf-8")
                      tests[k]['stderr'] = stderr.decode("utf-8")
                      tests[k]['exit_code'] = proc.returncode

                      k+=1
                                        
    with open(WORKSPACE+'/tests.json','w')as f:          
      f.write(json.dumps(tests,indent=2))
      
    print(json.dumps(tests,indent=2))      

if __name__=='__main__':
    main()
