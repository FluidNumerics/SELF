#!/usr/bin/env python3

import errno
import subprocess
import os
import time
from datetime import datetime
import json
import sys


# Get the path to this script's directory
cipath = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/'
cifile = cipath+'fluid-run.json'
os.environ['WORKSPACE']=cipath

def localRun():
    """Executes all execution_commands sequentially on your local system"""

    try:
        with open(cifile,'r')as f: 
            tests = json.load(f)
    except:
        print('Error opening CI file {}'.format(cifile))
        sys.exit(-1)
    utc = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
    k=0
    for test in tests['tests'] :

        workdir=cipath
        os.chdir(workdir)

        cmd = test['execution_command']

        print('Running {}\n'.format(cmd),flush=True)
        t0 = time.time()
        proc = subprocess.Popen(cmd,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        t1 = time.time()
        print(stdout.decode("utf-8"),flush=True)
        print(stderr.decode("utf-8"),flush=True)
        tests['tests'][k]['stdout'] = stdout.decode("utf-8")
        tests['tests'][k]['stderr'] = stderr.decode("utf-8")
        tests['tests'][k]['exit_code'] = proc.returncode
        tests['tests'][k]['build_id'] = ""
        tests['tests'][k]['machine_type'] = ""
        tests['tests'][k]['node_count'] = 1
        tests['tests'][k]['gpu_type'] = ""
        tests['tests'][k]['gpu_count'] = ""
        tests['tests'][k]['git_sha'] = ""
        tests['tests'][k]['datetime'] = utc
        tests['tests'][k]['runtime'] = float(t1-t0)
        tests['tests'][k]['allocated_cpus'] = ""
        tests['tests'][k]['compiler'] = ""
        tests['tests'][k]['target_arch'] = ""

        k+=1
                                        
    return tests

#END localRun

def checkExitCodes(tests):
    """Parses the results.json output and reports exit code statistics"""

    results = {}
    sysExitCode = 0
    for test in tests['tests'] :
        cli_command = test['command_group']
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

    print('============================',flush=True)
    print('',flush=True)
    print('Exit Code Check',flush=True)
    print('')
    for cli in results.keys():
      npass = results[cli]['npass']
      nfail = results[cli]['nfail']
      print('============================',flush=True)
      print('                        ',flush=True)
      print('  {}'.format(cli),flush=True)
      print('  > PASS : {}/{}'.format(str(npass),str(npass+nfail)),flush=True)
      print('  > FAIL : {}/{}'.format(str(nfail),str(npass+nfail)),flush=True)
      print('  > PASS RATE : {}%'.format(str(npass/(npass+nfail)*100)),flush=True)
      print('                        ',flush=True)

    print('============================',flush=True)

    sys.exit(sysExitCode)

#END checkExitCodes

def main():


  res = localRun()
  checkExitCodes(res)


if __name__=='__main__':
    main()
