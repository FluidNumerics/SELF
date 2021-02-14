#!/usr/bin/python3

import json


def load_config():

    with open('tests.json','r') as f:
        config = json.load(f)
    return config

# load_config

def expand_scalar(function, derivative="", gradient=[]):

    cmd_flags = '--function "{}" \\\n'.format(function)

    if derivative:
      cmd_flags += '--derivative "{}" \\\n'.format(derivative)

    if gradient:
      if len(gradient) == 2:
        cmd_flags += '--vector-x "{}" \\\n--vector-y "{}" \\\n'.format(gradient[0],gradient[1])

      elif len(gradient) == 3:
        cmd_flags += '--vector-x "{}" \\\n--vector-y "{}" \\\n--vector-z "{}" \\\n'.format(gradient[0],gradient[1],gradient[2])

    return cmd_flags

#END expand_scalar

def expand_vector(function, divergence="", gradient=[]):

    if len(function) == 2:
      cmd_flags = '--vector-x "{}" \\\n--vector-y "{}" \\\n'.format(function[0],function[1])

    elif len(function) == 3:
      cmd_flags = '--vector-x "{}" \\\n--vector-y "{}" \\\n--vector-z "{}" \\\n'.format(function[0],function[1],function[2])

    if divergence:
      cmd_flags += '--function "{}" \\\n'.format(divergence)

    if gradient:
      if len(gradient) == 2:
        cmd_flags += '--tensor-11 "{}" \\\n--tensor-12 "{}" \\\n--tensor-21 "{}" \\\n --tensor-22 "{}" \\\n'.format(gradient[0][0],gradient[0][1],
                                                                 gradient[1][0],gradient[1][1])

      elif len(gradient) == 3:
        cmd_flags += '--tensor-11 "{}" \\\n--tensor-12 "{}" \\\n--tensor-13 "{}" \\\n--tensor-21 "{}" \\\n--tensor-22 "{}" \\\n--tensor-23 "{}" \\\n--tensor-31 "{}" \\\n--tensor-32 "{}" \\\n--tensor-33 "{}" \\\n '.format(gradient[0][0],gradient[0][1],gradient[0][2],
                                          gradient[1][0],gradient[1][1],gradient[1][2],
                                          gradient[2][0],gradient[2][1],gradient[2][2],)

    return cmd_flags

#END expand_vector

def expand_tensor(function, divergence=[]):


    if len(function) == 2:
      cmd_flags = '--tensor-11 "{}" \\\n--tensor-12 "{}" \\\n--tensor-21 "{}" \\\n--tensor-22 "{}" \\\n '.format(function[0][0],function[0][1],
                                                               function[1][0],function[1][1])

    elif len(function) == 3:
      cmd_flags = '--tensor-11 "{}" \\\n--tensor-12 "{}" \\\n--tensor-13 "{}" \\\n--tensor-21 "{}" \\\n--tensor-22 "{}" \\\n--tensor-23 "{}" \\\n--tensor-31 "{}" \\\n--tensor-32 "{}" \\\n--tensor-33 "{}" \\\n'.format(function[0][0],function[0][1],function[0][2],
                                        function[1][0],function[1][1],function[1][2],
                                        function[2][0],function[2][1],function[2][2],)

    if len(divergence) == 2:
      cmd_flags += '--vector-x "{}" \\\n--vector-y "{}" \\\n'.format(divergence[0],divergence[1])

    elif len(divergence) == 3:
      cmd_flags += '--vector-x "{}" \\\n--vector-y "{}" \\\n--vector-z "{}" \\\n'.format(divergence[0],divergence[1],divergence[2])

    return cmd_flags

#END expand_tensor

def expand_tolerance(tol0,deg0,degree,etype):

    if etype == 'exact':
      tol = tol0
    elif etype == 'exponential':
      tol = float(tol0)*10.0**(float(deg0)-float(degree))
    return tol

#END expand_tolerance

def expand_function(fset):

    if 'scalar_function' in fset.keys():

      if 'gradient' in fset.keys():
        fcmd = expand_scalar(function = fset['scalar_function'],
                             gradient = fset['gradient'])

      elif 'derivative' in fset.keys():
        fcmd = expand_scalar(function = fset['scalar_function'],
                             derivative = fset['derivative'])

      else:
        fcmd = expand_scalar(function = fset['scalar_function'])

    elif 'vector2d_function' in fset.keys():

      if 'gradient' in fset.keys():
        fcmd = expand_vector(function = fset['vector2d_function'],
                             gradient = fset['gradient'])

      elif 'divergence' in fset.keys():
        fcmd = expand_vector(function = fset['vector2d_function'],
                             divergence = fset['divergence'])

      else:
        fcmd = expand_vector(function = fset['vector2d_function'])

    elif 'vector3d_function' in fset.keys():

      if 'gradient' in fset.keys():
        fcmd = expand_vector(function = fset['vector3d_function'],
                             gradient = fset['gradient'])

      elif 'divergence' in fset.keys():
        fcmd = expand_vector(function = fset['vector3d_function'],
                             divergence = fset['divergence'])

      else:
        fcmd = expand_vector(function = fset['vector3d_function'])

    elif 'tensor2d_function' in fset.keys():

      if 'divergence' in fset.keys():
        fcmd = expand_tensor(function = fset['tensor2d_function'],
                             divergence = fset['divergence'])

      else:
        fcmd = expand_tensor(function = fset['tensor2d_function'])

    elif 'tensor3d_function' in fset.keys():

      if 'divergence' in fset.keys():
        fcmd = expand_tensor(function = fset['tensor3d_function'],
                             divergence = fset['divergence'])

      else:
        fcmd = expand_tensor(function = fset['tensor3d_function'])

    return fcmd

#END expand_function

def main():

    config = load_config()

    ci_sh = """#!/bin/bash
COUNTER=0
trap '(( $? && ++errcount))' DEBUG

{CMDS}

echo "============================="
echo "        Test Summary         "
echo "============================="
echo " "
echo " Total Tests : $COUNTER "
echo " Fail : $errcount "
test $errcount = 0
"""


    crange = config['control_range']
    tdeg = config['target_range'][0]
    nel = config['elements_range'][0]
    nvar = config['var_range'][0]

    cmds = ""

    for gpu in config['gpu_accel']:
      for cq in config['control_quadratures']:
        for tq in config['target_quadratures']:
          for test in config['tests']:
            for fset in test['fsets']:

              if 'derivative_types' in fset.keys():
                for dt in fset['derivative_types']:
                  for cdeg in range(crange[0],crange[1]+1):
                    tol = expand_tolerance(fset['tolerance'],crange[0],cdeg,fset['error_type'])
                    cmd = "let COUNTER++ \n"
                    cmd += "${INSTALL_ROOT}/bin/self "
                    cmd +='--tolerance "{}" \\\n'.format(str(tol))
                    cmd +='--control-quadrature "{}" \\\n'.format(cq)
                    cmd +='--control-degree {} \\\n'.format(cdeg)
                    cmd +='--target-quadrature "{}" \\\n'.format(tq)
                    cmd +='--target-degree {} \\\n'.format(tdeg)
                    cmd +='--derivative-type "{}" \\\n'.format(dt)
                    cmd +='--nelements {} \\\n'.format(nel)
                    cmd +='--nvar {} \\\n'.format(nvar)
                    cmd +='--gpu-accel "{}" \\\n'.format(gpu)
                    fcmd = expand_function(fset)
                    cmd += fcmd
                    cmd += test['command']
                    cmds +=cmd + '\n\n'

              else:
                for cdeg in range(crange[0],crange[1]+1):
                  tol = expand_tolerance(fset['tolerance'],crange[0],cdeg,fset['error_type'])
                  cmd = "let COUNTER++ \n"
                  cmd += "${INSTALL_ROOT}/bin/self "
                  cmd +='--tolerance "{}" \\\n'.format(str(tol))
                  cmd +='--control-quadrature "{}" \\\n'.format(cq)
                  cmd +='--control-degree {} \\\n'.format(cdeg)
                  cmd +='--target-quadrature "{}" \\\n'.format(tq)
                  cmd +='--target-degree {} \\\n'.format(tdeg)
                  cmd +='--nelements {} \\\n'.format(nel)
                  cmd +='--nvar {} \\\n'.format(nvar)
                  cmd +='--gpu-accel "{}" \\\n'.format(gpu)
                  fcmd = expand_function(fset)
                  cmd += fcmd
                  cmd += test['command']
                  cmds +=cmd + '\n\n'


    # Write to file
    ci_sh = ci_sh.replace('{CMDS}',cmds)
    with open('ci.sh','w') as f:
        f.write(ci_sh)


if __name__ == '__main__':
    main()
