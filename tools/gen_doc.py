#!/usr/bin/python3

import json


SRC_ROOT = '../'
SRC_FILES = ['src/self/Lagrange_Class.F03']

# docs:
# - file: string
#   dependencies: list(string)
#   datatypes: 
#     - name: string
#       attributes:
#         - name: string
#           type: string
#           shape: string
#           description: list(string)
#
#       procedures: 
#         - name: string
#           visibility: string, enum(private | public)
#           type: string, enum(generic | procedure)
#           description: list(string)
#           line_bounds: list(number)
#           input :
#             - name: string
#               type: string
#               shape: string
#   

def declaration(line):

    if 'REAL' in line or 'INTEGER' in line or 'PROCEDURE' in line or 'GENERIC' in line :
        return True
    else:
        return False

#END declaration

typeOpen = False
routineOpen = False
descriptorOpen = False

for srcFile in SRC_FILES:

    doc = {"file":srcFile,
           "dependencies":[],
           "datatypes": [] }

    with open(SRC_ROOT+srcFile, "r") as f:
        line = f.readline()
        cnt = 1
        while line:
            line = f.readline()
            cnt += 1

            if typeOpen:
                if '!>' in line:
                    descriptor = line.rstrip().split(' : ')
                    if 'attribute' in descriptor[0]:
                        nameShape = descriptor[1].split('(')
                        attribute = {'name':nameShape[0],
                                     'description':[descriptor[2]],
                                     'shape':'('+nameShape[1],
                                     'type':''}
                        dataType['attributes'].append(attribute)

                elif 'END' in line:
                    doc['datatypes'].append(dataType)
                    typeOpen = False

           
                elif declaration(line.rstrip()):
                    # In this case, a declaration has been made for an attribute or procedure

                    if 'PROCEDURE' in line.rstrip() or 'GENERIC' in line.rstrip():

                        declared = line.rstrip().split(' :: ')
                        ptypeVis = declared[0].lstrip().rstrip()
                        ptype = ptypeVis.split(', ')[0]
                        pvis = ptypeVis.split(', ')[1]
                        pname = declared[1].split(' => ')[0]
                        if len(declared[1].split(' => ')) > 1:
                            pmaps = declared[1].split(' => ')[1].split(', ')
                        else:
                            pmaps = []

                        procedure = {'name': pname,
                                     'visibility': pvis, 
                                     'type': ptype,
                                     'maps_to': pmaps,
                                     'line_bounds': [],
                                     'description': [],
                                     'input':[],
                                     'output':[]}
                        dataType['procedures'].append(procedure)

                    else:
                        # This is a data attribute declaration
                        # We now want to assign the datatype
                        declared = line.rstrip().split(' :: ')
                        dtype = declared[0].lstrip().rstrip()
                        dname = declared[1].split('(')[0]
                        k = 0
                        for attr in dataType['attributes']:
                           if dname == attr['name']:
                               dataType['attributes'][k]['type'] = dtype
                               break
                           k+=1

            elif routineOpen:
                print('Routine Open Stub')

            else:
                if '!>' in line:
                    print('Descriptor found outside of routine or type')
                elif 'USE' in line:
                    dependency = line.split(' ')[1].rstrip()
                    if dependency :
                        doc['dependencies'].append(dependency)

                elif 'TYPE, PUBLIC' in line:
                    typeOpen = True
                    thisType = line.split(' :: ')[1].rstrip()
                    print('Open type '+thisType)
                    dataType = {'name':thisType,
                                'attributes':[],
                                'procedures':[]} 
                        
                   

print(json.dumps(doc,indent=2,sort_keys=True))
                    
