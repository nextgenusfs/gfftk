import errno
import os
import re


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


def filter_annotations(annotations, grep=None, grepv=None):
    """Filter annotations based on grep and grepv patterns.

    Parameters
    ----------
    annotations : dict
        GFFtk annotation dictionary keyed by gene ID
    grep : list, optional
        List of patterns to keep (format: "key:pattern" or "key:pattern:flags")
    grepv : list, optional
        List of patterns to remove (format: "key:pattern" or "key:pattern:flags")

    Returns
    -------
    dict
        Filtered annotation dictionary

    Examples
    --------
    # Keep only genes with "kinase" in product
    filtered = filter_annotations(genes, grep=["product:kinase"])

    # Remove genes from specific source
    filtered = filter_annotations(genes, grepv=["source:augustus"])

    # Case-insensitive matching
    filtered = filter_annotations(genes, grep=["product:kinase:i"])
    """
    if not grep and not grepv:
        return annotations

    if grep is None:
        grep = []
    if grepv is None:
        grepv = []

    # Parse grep patterns
    grep_patterns = []
    for pattern in grep:
        grep_patterns.append(_parse_filter_pattern(pattern))

    # Parse grepv patterns
    grepv_patterns = []
    for pattern in grepv:
        grepv_patterns.append(_parse_filter_pattern(pattern))

    filtered_annotations = {}

    for gene_id, gene_data in annotations.items():
        # Check if gene should be kept based on grep patterns
        keep_gene = True

        # If grep patterns exist, gene must match at least one
        if grep_patterns:
            keep_gene = False
            for key, pattern, flags in grep_patterns:
                if _match_gene_pattern(gene_data, key, pattern, flags):
                    keep_gene = True
                    break

        # If gene passes grep filter, check grepv patterns
        if keep_gene and grepv_patterns:
            for key, pattern, flags in grepv_patterns:
                if _match_gene_pattern(gene_data, key, pattern, flags):
                    keep_gene = False
                    break

        if keep_gene:
            filtered_annotations[gene_id] = gene_data

    return filtered_annotations


def _parse_filter_pattern(pattern_str):
    """Parse a filter pattern string into components.

    Parameters
    ----------
    pattern_str : str
        Pattern in format "key:pattern" or "key:pattern:flags"

    Returns
    -------
    tuple
        (key, pattern, flags) where flags is a string of regex flags
    """
    parts = pattern_str.split(":", 2)
    if len(parts) < 2:
        raise ValueError(
            f"Invalid filter pattern: {pattern_str}. Expected format: 'key:pattern' or 'key:pattern:flags'"
        )

    key = parts[0]
    pattern = parts[1]
    flags = parts[2] if len(parts) > 2 else ""

    return key, pattern, flags


def _match_gene_pattern(gene_data, key, pattern, flags):
    """Check if a gene matches a specific pattern.

    Parameters
    ----------
    gene_data : dict
        Gene annotation data
    key : str
        Key to search in (e.g., 'product', 'source', 'note')
    pattern : str
        Regular expression pattern to match
    flags : str
        Regex flags (i=ignore case, m=multiline, etc.)

    Returns
    -------
    bool
        True if pattern matches, False otherwise
    """
    # Convert flags string to re flags
    re_flags = 0
    if "i" in flags.lower():
        re_flags |= re.IGNORECASE
    if "m" in flags.lower():
        re_flags |= re.MULTILINE
    if "s" in flags.lower():
        re_flags |= re.DOTALL

    # Get the value(s) to search
    if key not in gene_data:
        return False

    value = gene_data[key]

    # Handle different value types
    search_strings = []
    if isinstance(value, str):
        search_strings = [value]
    elif isinstance(value, list):
        # Flatten nested lists and convert to strings
        for item in value:
            if isinstance(item, list):
                search_strings.extend([str(x) for x in item if x])
            else:
                if item:  # Skip empty/None values
                    search_strings.append(str(item))
    else:
        search_strings = [str(value)]

    # Check if pattern matches any of the search strings
    try:
        compiled_pattern = re.compile(pattern, re_flags)
        for search_str in search_strings:
            if compiled_pattern.search(search_str):
                return True
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{pattern}': {e}")

    return False
