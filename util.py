import os
import gzip
import time
import subprocess


def remove_safely(fn):
    if os.path.isfile(fn):
        os.remove(fn)


def smart_open(filename, mode):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def smarter_open(filename, mode, quiet=False):
    assert mode == "r"
    if filename.startswith("s3://"):
        s3cat = f"aws s3 cp {filename} -"
        if filename.endswith(".gz"):
            s3cat += " | gzip -dc"
        if not quiet:
            print(repr(s3cat))
        # TODO: check=True?
        return subprocess.Popen(s3cat, shell=True, stdout=subprocess.PIPE)
    return smart_open(filename, mode)


def smarter_readline(f):
    if hasattr(f, "stdout"):
        return f.stdout.readline()
    return f.readline()


def chop(txt, suffix):
    assert txt.endswith(suffix)
    return txt[:-len(suffix)]


def check_call(command, capture_stdout=False, quiet=False):
    # Assuming python >= 3.5
    shell = isinstance(command, str)
    if not quiet:
        command_str = command if shell else " ".join(command)
        print(repr(command_str))
    # In Python 3.7 the subprocess.run function accepts capture_output param,
    # so setting stdout like so is only done for backward compatibility with
    # python versions >= 3.5.  That's the minimum we support.
    stdout = subprocess.PIPE if capture_stdout else None
    p = subprocess.run(command, shell=shell, check=True, stdout=stdout)
    return p.stdout.decode('utf-8') if capture_stdout else None


def check_output(command, quiet=False):
    return check_call(command, capture_stdout=True, quiet=quiet)


class ProgressTracker:

    def __init__(self, target):
        self.target = target
        self.current = 0
        self.t_start = time.time()

    def advance(self, amount):
        PESSIMISM = 2.0
        self.current += amount
        t_elapsed = time.time() - self.t_start
        t_remaining = (t_elapsed / self.current) * self.target - t_elapsed
        t_remaining *= PESSIMISM
        t_eta = self.t_start + t_elapsed + t_remaining
        t_eta_str = time.strftime("%H:%M:%S", time.localtime(t_eta))
        print(f"*** {self.current/self.target*100:3.1f} percent done, {t_elapsed/60:3.1f} minutes elapsed, {t_remaining/60:3.1f} minutes remaining, ETA {t_eta_str} ***\n")
