import os
import neoqm

ROOTDIR = os.path.dirname(os.path.dirname(neoqm.__file__))
DATADIR = os.path.join(ROOTDIR, "data")

import subprocess



def run_command(command, stdout_file=None, stderr_file=None):
    if stdout_file is None:
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    else:
        assert isinstance(stdout_file, str)
        outputf = open(stdout_file, "w")

    while True:
        line = process.stdout.readline().rstrip().decode("utf-8")
        if line == '' and process.poll() is not None:
            break

        if stdout_file:
            outputf.write(line)
        else:
            print(line)
    return process

class RuntimeTempDir():

    def __init__(self, directory):
        if directory.endswith("/"):
            directory = directory[:-1]
        self.tmpdir = tempfile.TemporaryDirectory(
            dir=os.path.dirname(directory))
        self.directory = directory

    def __enter__(self):
        return self.tmpdir

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        for fname in os.listdir(self.tmpdir.name):
            os.rename(os.path.join(self.tmpdir.name, fname),
                      os.path.join(self.directory, fname))
        del self.tmpdir

