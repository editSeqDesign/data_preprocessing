#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-02-16 06:17:32
LastEditors: wangruoyu
LastEditTime: 2023-02-16 06:17:51
Description: file content
FilePath: /chopchop_crispr_cdk/lambda/data_preprocessing/data_preprocessing/io_sys.py
'''
import subprocess
import sys
from os.path import exists, dirname

dirs2ps={'pyp':str(subprocess.check_output('which python'.split(' '))).replace("b'",'').replace("\\n'",''),
#'binDir':dirname(str(subprocess.check_output('which bwa'.split(' '))).replace("b'",'').replace("\\n'",'')),
'binDir':dirname(__file__)+'/bin', 
'scriptDir':dirname(__file__)+'/bin',
}   

def runbashcmd(cmd,test=False,logf=None):
    # from super_beditor.lib.global_vars import dirs2ps 

    print('before:',cmd)
    # cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    # cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    # cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
    print('after:',cmd)
    if test:
        print(cmd)
    err=subprocess.call(cmd,shell=True,stdout=logf,stderr=subprocess.STDOUT)
    print(err)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)

def is_interactive():
    # thanks to https://stackoverflow.com/a/22424821/3521099
    import __main__ as main
    return not hasattr(main, '__file__')
