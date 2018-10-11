import os
import gzip


def remove_safely(fn):
    if os.path.isfile(fn):
        os.remove(fn)


def check_call(command, quiet=False):
    if not quiet:
        print(repr(command))
    exitcode = os.system(command)
    assert exitcode == 0, f"Command failed: {command}"


def smart_open(filename, mode):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)
