import os
import stat
import subprocess

def write_script_to_file(filename, multiline_string):
    lines = multiline_string.splitlines()        
    stripped_lines = [line.strip() for line in lines]        
    stripped_text = '\n'.join(stripped_lines)        
    with open(filename, 'w') as file:
        file.write(stripped_text)
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

def bash_run(cmd):
    subprocess.run("LD_PRELOAD=/usr/lib/nsight-systems/host-linux-x64/Mesa/libGL.so:$LD_PRELOAD " + cmd,
                   shell=True,
                   cwd=os.getcwd(),
                   executable='/bin/bash')
