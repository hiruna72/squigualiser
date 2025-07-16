#!/bin/bash

die (){
    echo >&2 "$@"
    exit 1
}

PYTHON_VERSION="python3.9"
PY_VENV=""
ARCH=$(uname -m)
OS=$(uname -s)

echo "O/S:${OS} architecture:${ARCH} python:${PYTHON_VERSION}"

if [ "${OS}" == "Linux"  ] && [ "${ARCH}" == "x86_64" ];
then
    apt-get update || die "apt-get update failed"
    apt install wget gcc make zlib1g-dev git -y || die "system tools install failed"
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.9.23+20250712-x86_64-unknown-linux-gnu-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.9.23+20250712-x86_64-unknown-linux-gnu-install_only.tar.gz || die "untar python failed"
elif [[ "${OS}" == "Darwin" && ( "${ARCH}" == "arm64" || "${ARCH}" == "aarch64" ) ]];
then
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.9.23+20250712-aarch64-apple-darwin-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.9.23+20250712-aarch64-apple-darwin-install_only.tar.gz || die "untar python failed"
elif [ "${OS}" == "Darwin"  ] && [ "${ARCH}" == "x86_64" ];
then
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.9.23+20250712-x86_64-apple-darwin-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.9.23+20250712-x86_64-apple-darwin-install_only.tar.gz || die "untar python failed"
else
    die "Unsupported O/S ${OS} or architecture ${ARCH} for packaging."
fi

python/bin/${PYTHON_VERSION} -m venv ${PY_VENV} || die "create venv failed"
source ${PY_VENV}/bin/activate || die "sourcing venv failed"
pip install --upgrade pip || die "upgrade pip failed"
export CC=gcc
export HTSLIB_CONFIGURE_OPTIONS="--enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no"

if [[ "$1" == "test_pypi" ]]; then
    echo "Installing from Test PyPI"
    pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple squigualiser --pre || die "test.pip install squigualiser failed"
else
    echo "Installing from PyPI"
    pip install squigualiser --no-cache || die "pip install squigualiser failed"
fi

find ./ -name __pycache__ -type d | xargs rm -r || die "removing pycache failed"
mv ${PY_VENV}/bin/squigualiser python/bin/ || die "moving squigualiser to bin failed"
cp -r ${PY_VENV}/lib/${PYTHON_VERSION}/site-packages/* python/lib/${PYTHON_VERSION}/site-packages/ || die "copying site-packages failed"

if [ "${OS}" == "Linux"  ]; then
    sed -i "s/${PY_VENV}\/bin\/${PYTHON_VERSION}/\/usr\/bin\/env ${PYTHON_VERSION}/g" python/bin/squigualiser  || die "changing headerline failed"
elif [ "${OS}" == "Darwin"  ]; then
    sed -i '' "1s/.*/#\!\/usr\/bin\/env ${PYTHON_VERSION}/" python/bin/squigualiser || die "changing headerline failed"
fi

git clone --depth 1 --branch main https://github.com/hiruna72/squigualiser.git  || die "Failed to clone squigualiser"
cp -r squigualiser/docs python || die "docs copy failed"
cp squigualiser/test/package/squigualiser python || die "script copy failed" 
cp squigualiser/LICENSE python || die "license copy failed"
cp squigualiser/README.md python || die "readme copy failed"

rm -rf squigualiser || die "remove cloned dir failed"

mv python/ squigualiser || die "renaming python to squigualiser failed"

tar zcvf squigualiser.tar.gz squigualiser/ || die "tar balling squigualiser failed"

# if user arg "docker" is provided, copy tarball to host directory
if [[ "$2" == "docker" ]]; then
    echo "copying tar file to host directory"
    cp squigualiser.tar.gz /host/ || die "copying tar file to host"
fi
