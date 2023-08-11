import setuptools
from os import path
import re

PKG_NAME = "squigualiser"
MOD_NAME = "src"

# add readme to long description as that's what pypi sees
with open("README.md", "r") as f:
    long_description = f.read()

# get version from file rather than here so change isn't in this file
__version__ = ""
exec(open("{}/_version.py".format(MOD_NAME)).read())

# User can set version of ont-pyguppy-client-lib to match guppy version
with open(path.join(path.abspath(path.dirname(__file__)), "requirements.txt")) as f:
    install_requires = [p.strip() for p in f]


setuptools.setup(
    name=PKG_NAME,
    version=__version__,
    url="https://github.com/hiruna72/squigualiser",
    author="Hiruna Samarakoon",
    author_email="h.samarakoon@garvan.org.au",
    maintainer='Hiruna Samarakoon',
    maintainer_email='h.samarakoon@garvan.org.au',
    description="Visualise ONT raw signals ",
    license='MIT',
    keywords=['nanopore', 'slow5', 'signal'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    python_requires=">=3.8",
    install_requires=install_requires,
    setup_requires=["numpy"],
    entry_points={
        "console_scripts": [
            "squigualiser={}.__init__:main".format(MOD_NAME),
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
