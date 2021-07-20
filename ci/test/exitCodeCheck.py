#!/usr/bin/python


import json
import sys

WORKSPACE=os.getenv('WORKSPACE')


def main():

    with open(WORKSPACE+'/tests.json','r')as f:          
      tests = json.load(f)

    sysExitCode = 0
    for test in ci_conf['tests'] :
      npass = 0
      nfail = 0
      for nel in test['nelems'] :
        for nvar in test['nvars'] :
          for cQuad in test['control_quadrature'] :
            for cDeg in test['control_degree'] :
              for tQuad in test['target_quadrature'] :
                for tDeg in test['target_degree'] :
                  for addlOpts in test['additional_opts'] :
                    for funcOpts in test['function_opts'] :

                      if int(test['exit_code']) == 0:
                          npass += 1
                      else:
                          nfail += 1
                          sysExitCode=1

      print('========================')
      print('                        ')
      print('  {}'.format(test['cli_command']))
      print('  > PASS : {}/{}'.format(str(npass),str(npass+nfail)))
      print('  > FAIL : {}/{}'.format(str(fail),str(npass+nfail)))
      print('  > FAILURE RATE : {}%'.format(str(fail/(npass+nfail)*100)))
      print('                        ')

    sys.exit(sysExitCode) 
