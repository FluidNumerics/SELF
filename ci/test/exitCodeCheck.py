#!/usr/bin/python


import json
import sys
import os

WORKSPACE=os.getenv('WORKSPACE')


def main():

    with open(WORKSPACE+'/tests.json','r')as f:          
      tests = json.load(f)

    results = {}
    sysExitCode = 0
    for test in tests :
        cli_command = test['cli_command']
        if cli_command in results.keys():
            if test['exit_code'] == 0:
                results[cli_command]['npass'] += 1
            else:
                results[cli_command]['nfail'] += 1
                sysExitCode = 1
        else:
            results[cli_command] = {'npass':0,'nfail':0}
            if test['exit_code'] == 0:
                results[cli_command]['npass'] += 1
            else:
                results[cli_command]['nfail'] += 1
                sysExitCode = 1

    print('============================')
    print('')
    print('Program Exit Code Check')
    print('')
    for cli in results.keys():
      npass = results[cli]['npass']
      nfail = results[cli]['nfail']
      print('============================')
      print('                        ')
      print('  {}'.format(cli))
      print('  > PASS : {}/{}'.format(str(npass),str(npass+nfail)))
      print('  > FAIL : {}/{}'.format(str(nfail),str(npass+nfail)))
      print('  > PASS RATE : {}%'.format(str(npass/(npass+nfail)*100)))
      print('                        ')

    print('============================')
    sys.exit(sysExitCode) 

if __name__ == '__main__':
    main()
