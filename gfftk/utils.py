import errno
import os


def readBlocks(source, pattern):
    buffer = []
    for line in source:
        try:
            line = line.decode("utf-8")
        except AttributeError:
            line = line
        if line.startswith(pattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def readBlocks2(source, startpattern, endpattern):
    buffer = []
    for line in source:
        try:
            line = line.decode("utf-8")
        except AttributeError:
            line = line
        if line.startswith(startpattern) or line.endswith(endpattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def check_inputs(inputs):
    for filename in inputs:
        if not is_file(filename):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)


def is_file(f):
    if os.path.isfile(f):
        return True
    else:
        return False


def is_gzipped(filepath):
    """Checks if a file is gzipped."""
    with open(filepath, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def is_text_file(filepath):
    """Checks if a file is likely a text file."""
    try:
        with open(filepath, "r") as f:
            f.read(4096)
        return True
    except UnicodeDecodeError:
        return False


def check_file_type(filepath):
    """Checks if a file is text, gzipped binary, or other binary."""
    if is_gzipped(filepath):
        return "gzipped binary"
    elif is_text_file(filepath):
        return "text"
    else:
        return "binary"


def which2(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def open_pipe(command, mode="r", buff=1024 * 1024):
    import signal
    import subprocess

    if "r" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL),
        ).stdout
    elif "w" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            universal_newlines=True,
            stdin=subprocess.PIPE,
        ).stdin
    return None


NORMAL = 0
PROCESS = 1
PARALLEL = 2

WHICH_BZIP2 = which2("bzip2")
WHICH_PBZIP2 = which2("pbzip2")


def open_bz2(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import bz2

        return bz2.BZ2File(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_BZIP2:
            return open_bz2(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("bzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("bzip2 >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PBZIP2:
            return open_bz2(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pbzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pbzip2 >" + filename, mode, buff)
    return None


WHICH_GZIP = which2("gzip")
WHICH_PIGZ = which2("pigz")


def open_gz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import gzip

        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("gzip -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("gzip >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pigz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pigz >" + filename, mode, buff)
    return None


WHICH_XZ = which2("xz")


def open_xz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if WHICH_XZ:
        if "r" in mode:
            return open_pipe("xz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("xz >" + filename, mode, buff)
    return None


def zopen(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    """
    Open pipe, zipped, or unzipped file automagically

    # external == 0: normal zip libraries
    # external == 1: (zcat, gzip) or (bzcat, bzip2)
    # external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if "r" in mode and "w" in mode:
        return None
    if filename.startswith("!"):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith(".bz2"):
        return open_bz2(filename, mode, buff, external)
    elif filename.endswith(".gz"):
        return open_gz(filename, mode, buff, external)
    elif filename.endswith(".xz"):
        return open_xz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
    return None
