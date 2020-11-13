#!/usr/bin/python3

import git
import subprocess
import json


repo = git.Repo(search_parent_directories=True)
sha = repo.head.object.hexsha
print(sha)

proc = subprocess.Popen(['./test'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = proc.communicate()

if proc.returncode == 0 :
    tests = stdout.decode('utf-8').split('}')
    data = json.loads( tests[0]+'}')
    data['commit_sha'] = sha
    print( data )
else :
    print( stderr )


