#!/usr/bin/python3

import json
import os

WORKSPACE=os.getenv('WORKSPACE')
GPU_ACCEL=os.getenv('GPU_ACCEL')
SCHEDULER=os.getenv('SCHEDULER')


def main():

    with open(WORKSPACE+'/test/ci.json','r') as f:
        ci_conf = json.load(f)

    tests = []
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

                        workdir = WORKSPACE+'/test/'
                        workdir += test['cli_command']+'/'
                        workdir += 'nel_{}'.format(nel)+'/'
                        workdir += 'nvar_{}'.format(nvar)+'/'
                        workdir += 'cQuad_{}'.format(cQuad)+'/'
                        workdir += 'cDeg_{}'.format(cDeg)+'/'
                        workdir += 'tQuad_{}'.format(tQuad)+'/'
                        workdir += 'tDeg_{}'.format(tDeg)+'/'
                        workdir += '{}'.format(addlOpts['name'])+'/'
                        workdir += '{}'.format(funcOpts['name'])+'/'

                        if SCHEDULER=='none'
                            proc = subprocess.Popen([workdir+'/test.sh'],shell=True,)
                            with open(workdir+'test.sh','w') as f:
                            f.write(cmd)
                        os.chmod(workdir+'test.sh',0o755)

                        tests.append({'test_script':workdir+'test.sh',
                                       'status':'ready',
                                       'exit_code':999,
                                       'profile':{},
                                       'system': {},
                                       'benchmark_info':{
                                         'gpu_accelerated':GPU_ACCEL,
                                         'runtime_seconds':0.0,
                                         'control_degree':cDeg,
                                         'control_quadrature':cQuad,
                                         'target_degree':tDeg,
                                         'target_quadrature':tQuad,
                                         'nelements':nel,
                                         'nvar':nvar,
                                         'cli_command':test['cli_command'],
                                         'additional_opts':addlOpts,
                                         'function_opts':funcOpts
                                       }
                                    })

    with open(WORKSPACE+'/tests.json','w')as f:          
      f.write(json.dumps(tests,indent=2))




                  



if __name__=='__main__':
    main()
