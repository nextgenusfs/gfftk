import os

from .utils import readBlocks, zopen


def go_term_dict():
    """Parses the Gene Ontology file and returns GO lookup dictionary, version, and date.

    The go.obo.gz file is included in the GFFtk package, it can be updated in site-packages/gfftk/data/go.obo.gz

    Returns
    -------
    go : dict of dict
        dictionary of go_term: { "name": str, "namespace": str}
    format_version : str
        version of the go.obo file
    release_date : str
        release data of the go.obo file

    """
    goref = os.path.join(os.path.join(os.path.dirname(__file__)), "data", "go.obo.gz")
    go = {}
    format_version = ""
    release_date = ""
    with zopen(goref) as infile:
        for i, x in enumerate(readBlocks(infile, "[Term]")):
            if i == 0:
                for item in x:
                    if item.startswith("format-version: "):
                        format_version = item.split("format-version: ", 1)[-1].rstrip()
                    elif item.startswith("data-version: releases/"):
                        release_date = item.split("data-version: releases/", 1)[-1].rstrip()
            go_term, go_name, go_space = (None,) * 3
            for item in x:
                if item.startswith("id: GO"):
                    go_term = item.split("id: ", 1)[-1].rstrip()
                elif item.startswith("name: "):
                    go_name = item.split("name: ", 1)[-1].rstrip()
                elif item.startswith("namespace: "):
                    go_space = item.split("namespace: ", 1)[-1].rstrip()
                elif item.startswith("\n"):
                    break
            if go_term and go_name and go_space:
                go[go_term] = {"name": go_name, "namespace": go_space}
    return go, format_version, release_date
