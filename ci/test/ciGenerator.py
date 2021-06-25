#!/usr/bin/python3

import json
import os

INSTALL_ROOT=os.getenv('INSTALL_ROOT')
WORKSPACE=os.getenv('WORKSPACE')
GPU_ACCEL=os.getenv('GPU_ACCEL')


def main():

  with open(INSTALL_ROOT+'/test/ci.json','r') as f:
    ci_conf = json.load(f)
        
  with open(INSTALL_ROOT+'/test/cmd.tmpl') as f:
    cmd_tmpl = f.read()


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

                    outdir = WORKSPACE+'/test/'
                    outdir += test['cli_command']+'/'
                    outdir += 'nel_{}'.format(nel)+'/'
                    outdir += 'nvar_{}'.format(nvar)+'/'
                    outdir += 'cQuad_{}'.format(cQuad)+'/'
                    outdir += 'cDeg_{}'.format(cDeg)+'/'
                    outdir += 'tQuad_{}'.format(tQuad)+'/'
                    outdir += 'tDeg_{}'.format(tDeg)+'/'
                    outdir += '{}'.format(addlOpts['name'])+'/'
                    outdir += '{}'.format(funcOpts['name'])+'/'

                    cmd = cmd_tmpl

                    cmd = cmd.replace('@PROFILER@',workdir)
                    cmd = cmd.replace('@WORKDIR@',workdir)
                    cmd = cmd.replace('@GPU_ACCEL@',GPU_ACCEL)
                    cmd = cmd.replace('@CONTROL_QUADRATURE@',cQuad)
                    cmd = cmd.replace('@CONTROL_DEGREE@',str(cDeg))
                    cmd = cmd.replace('@TARGET_QUADRATURE@',tQuad)
                    cmd = cmd.replace('@TARGET_DEGREE@',str(tDeg))
                    cmd = cmd.replace('@NELEMS@',str(nel))
                    cmd = cmd.replace('@NVAR@',str(nvar))
                    cmd = cmd.replace('@FUNCTION_OPTS@',funcOpts['value'])
                    cmd = cmd.replace('@ADDITIONAL_OPTS@',addlOpts['value'])
                    cmd = cmd.replace('@OUTPUT_FILE@',outdir+'self.h5')
                    cmd = cmd.replace('@COMMAND@',test['cli_command'])

                    os.makedirs(workdir)
                    os.makedirs(outdir)

                    with open(workdir+'test.sh','w') as f:
                      f.write(cmd)
                    os.chmod(workdir+'test.sh',0o755)

                    tests.append({'test_script':workdir+'test.sh',
                                  'status':'ready',
                                  'exit_code':999,
                                  'stdout': '',
                                  'stderr': '',
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

    with open(INSTALL_ROOT+'/tests.json','w')as f:          
      f.write(json.dumps(tests,indent=2))

if __name__=='__main__':
    main()
