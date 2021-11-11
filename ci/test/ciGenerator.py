#!/usr/bin/python3

import json
import os

# Get the path to this script's directory
cipath = '/'.join(os.path.abspath(__file__).split('/')[:-1])+'/'

PARTITIONS=['c2-standard-4','v100']


def main():

  with open(cipath+'ci.json','r') as f:
    ci_conf = json.load(f)
        
  with open(cipath+'cmd.tmpl') as f:
    cmd_tmpl = f.read()

  ntests = 0
  tests = {"tests":[]}
  GPU_ACCEL = 'false'

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

                    workdir = 'test/'
                    workdir += test['cli_command']+'/'
                    workdir += 'nel_{}'.format(nel)+'/'
                    workdir += 'nvar_{}'.format(nvar)+'/'
                    workdir += 'cQuad_{}'.format(cQuad)+'/'
                    workdir += 'cDeg_{}'.format(cDeg)+'/'
                    workdir += 'tQuad_{}'.format(tQuad)+'/'
                    workdir += 'tDeg_{}'.format(tDeg)+'/'
                    workdir += '{}'.format(addlOpts['name'])+'/'
                    workdir += '{}'.format(funcOpts['name'])+'/'

                    cmd = cmd_tmpl

                    cmd = cmd.replace('@PROFILER@','to-do')
                    cmd = cmd.replace('@GPU_ACCEL@',GPU_ACCEL)
                    cmd = cmd.replace('@CONTROL_QUADRATURE@',cQuad)
                    cmd = cmd.replace('@CONTROL_DEGREE@',str(cDeg))
                    cmd = cmd.replace('@TARGET_QUADRATURE@',tQuad)
                    cmd = cmd.replace('@TARGET_DEGREE@',str(tDeg))
                    cmd = cmd.replace('@NELEMS@',str(nel))
                    cmd = cmd.replace('@NVAR@',str(nvar))
                    cmd = cmd.replace('@FUNCTION_OPTS@',funcOpts['value'])
                    cmd = cmd.replace('@ADDITIONAL_OPTS@',addlOpts['value'])
                    cmd = cmd.replace('@OUTPUT_FILE@',workdir+'self.h5')
                    cmd = cmd.replace('@COMMAND@',test['cli_command'])

                    if GPU_ACCEL == 'true':
                        cmd = cmd.replace('@GPU_OPT@','--nv')
                    else:
                        cmd = cmd.replace('@GPU_OPT@','')

                    try:
                        os.makedirs(workdir)
                    except:
                        print(workdir + ' already exists. Overwriting tests')

                    with open(workdir+'test.sh','w') as f:
                      f.write(cmd)
                    os.chmod(workdir+'test.sh',0o755)

                    for partition in PARTITIONS:
                        tests["tests"].append({'command_group':test['cli_command'],
                                               'execution_command': workdir+'test.sh',
                                               'output_directory': workdir,
                                               'partition': partition})
                        ntests+=1

  with open('fluid-run.json','w')as f:          
    f.write(json.dumps(tests,indent=2))

  print('Generated {} tests for SELF'.format(ntests))

if __name__=='__main__':
    main()
