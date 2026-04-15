# Build script for the Cob compositional-loop Cython extension.
# Usage:   python3 setup_cob_ext.py build_ext --inplace
# Run from the KhT/ directory.

from setuptools import setup
from Cython.Build import cythonize

setup(
    name="cob_ext",
    ext_modules=cythonize(
        "cob_ext.pyx",
        compiler_directives={"language_level": "3"},
    ),
)
