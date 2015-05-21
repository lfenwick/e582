#!/usr/bin/env python
"""strip outputs from an IPython Notebook

Opens a notebook, strips its output, and writes the outputless version to the original file.

Useful mainly as a git pre-commit hook for users who don't want to track output in VCS.

This does mostly the same thing as the `Clear All Output` command in the notebook UI.

Adapted from rom https://gist.github.com/minrk/6176788 to work with
git filter driver

https://github.com/petered/plato/blob/fb2f4e252f50c79768920d0e47b870a8d799e92b/notebooks/config/strip_notebook_output
"""
import sys

#You may need to do this for your script to work with GitX or Tower:
#sys.path.append("/Users/chris/anaconda/envs/conda/lib/python2.7/site-packages")

try:
    from IPython.nbformat import v4
except ImportError:
    raise Exception("Failed to import the latest IPython while trying to strip output "
        "from your notebooks.  Either run venv/bin/activate to enter your virtual env, or update "
        "the IPython version on your machine (sudo pip install -U ipython)")
from IPython.nbconvert.preprocessors.clearoutput import ClearOutputPreprocessor

def strip_output(nb):
    """strip the outputs from a notebook object"""
    stripout=ClearOutputPreprocessor()
    nb,res=stripout.preprocess(nb,{})
    return nb

if __name__ == '__main__':
    nb = v4.reads(sys.stdin.read())
    nb = strip_output(nb)
    sys.stdout.write(v4.writes(nb))
