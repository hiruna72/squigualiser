import setuptools
from os import path
import re

PKG_NAME = "ideal-goggles"
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
    url="https://github.com/hiruna72/ideal-goggles",
    author="Hiruna Samarakoon",
    author_email="h.samarakoon@garvan.org.au",
    description="Visualise ONT raw signals ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    python_requires=">=3.7",
    install_requires=install_requires,
    setup_requires=["numpy"],
    entry_points={
        "console_scripts": [
            "ideal_goggles={}.__init__:main".format(MOD_NAME),
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)