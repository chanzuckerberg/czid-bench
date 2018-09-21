import os
import gzip
import traceback


def remove_safely(fn):
    try:
        if os.path.isfile(fn):
            os.remove(fn)
    except:
        print(f"Ignoring exception {traceback.format_exc()}")


def check_call(command, quiet=False):
    if not quiet:
        print(command)
    exitcode = os.system(command)
    assert exitcode == 0, f"Command failed: {command}"


def smart_open(filename, mode):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
