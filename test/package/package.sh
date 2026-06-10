#!/bin/bash

die (){
    echo >&2 "$@"
    exit 1
}

# python3.9 is EOL; pandas 2.3+ dropped py3.9 support and has no pre-built wheels for it
PYTHON_VERSION="python3.10"
PY_VENV="squig-venv"
ARCH=$(uname -m)
OS=$(uname -s)
REPO_LINK="https://github.com/hiruna72/squigualiser.git"
BRANCH="dev"
TOOL="squigualiser"

echo "O/S:${OS} architecture:${ARCH} python:${PYTHON_VERSION}"

if [ "${OS}" == "Linux"  ] && [ "${ARCH}" == "x86_64" ];
then
    apt-get update || die "apt-get update failed"
    apt install wget gcc make zlib1g-dev git -y || die "system tools install failed"
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.10.18+20250712-x86_64-unknown-linux-gnu-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.10.18+20250712-x86_64-unknown-linux-gnu-install_only.tar.gz || die "untar python failed"
elif [[ "${OS}" == "Darwin" && ( "${ARCH}" == "arm64" || "${ARCH}" == "aarch64" ) ]];
then
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.10.18+20250712-aarch64-apple-darwin-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.10.18+20250712-aarch64-apple-darwin-install_only.tar.gz || die "untar python failed"
elif [ "${OS}" == "Darwin"  ] && [ "${ARCH}" == "x86_64" ];
then
    wget https://github.com/indygreg/python-build-standalone/releases/download/20250712/cpython-3.10.18+20250712-x86_64-apple-darwin-install_only.tar.gz || die "python wget failed"
    tar xf cpython-3.10.18+20250712-x86_64-apple-darwin-install_only.tar.gz || die "untar python failed"
else
    die "Unsupported O/S ${OS} or architecture ${ARCH} for packaging."
fi

python/bin/${PYTHON_VERSION} -m venv ${PY_VENV} || die "create venv failed"
source ${PY_VENV}/bin/activate || die "sourcing venv failed"
pip install --upgrade pip || die "upgrade pip failed"
export CC=gcc
export HTSLIB_CONFIGURE_OPTIONS="--enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no"

git clone --depth 1 --branch ${BRANCH} ${REPO_LINK} || die "Failed to clone ${TOOL}"

if [[ "$2" == "test_pypi" ]]; then
    echo "Installing from source (${BRANCH} branch)"
    pip install ./${TOOL}/ --no-cache || die "pip install ${TOOL} from source failed"
else
    echo "Installing from PyPI"
    pip install ${TOOL} --no-cache || die "pip install ${TOOL} failed"
fi

find ./ -name __pycache__ -type d | xargs rm -r || die "removing pycache failed"
mv ${PY_VENV}/bin/${TOOL} python/bin/ || die "moving ${TOOL} to bin failed"
cp -r ${PY_VENV}/lib/${PYTHON_VERSION}/site-packages/* python/lib/${PYTHON_VERSION}/site-packages/ || die "copying site-packages failed"

if [ "${OS}" == "Linux"  ]; then
    sed -i "s/${PY_VENV}\/bin\/${PYTHON_VERSION}/\/usr\/bin\/env ${PYTHON_VERSION}/g" python/bin/${TOOL}  || die "changing headerline failed"
elif [ "${OS}" == "Darwin"  ]; then
    sed -i '' "1s/.*/#\!\/usr\/bin\/env ${PYTHON_VERSION}/" python/bin/${TOOL} || die "changing headerline failed"
fi

cp -r ${TOOL}/docs python || die "docs copy failed"
cp ${TOOL}/test/package/${TOOL} python || die "script copy failed"
cp ${TOOL}/LICENSE python || die "license copy failed"
cp ${TOOL}/README.md python || die "readme copy failed"

rm -rf ${TOOL} || die "remove cloned dir failed"

LATEST_TAG=$(git ls-remote --tags ${REPO_LINK} | cut -d/ -f3 | grep -v '\^{}' | sort -V | tail -n1)
OS_NAME="linux"
if [ "${OS}" == "Darwin"  ]; then
    OS_NAME="macos"
fi
ARCH_NAME="x86_64"
if [ "${ARCH}" == "arm64"  ] || [ "${ARCH}" == "aarch64" ]; then
    ARCH_NAME="arm64"
fi
TOOL_NAME=${TOOL}-${LATEST_TAG}
TAR_NAME=${TOOL}-${LATEST_TAG}-${ARCH_NAME}-${OS_NAME}-binaries.tar.gz
echo "TOOL_NAME: ${TOOL_NAME}"
echo "TAR_NAME: ${TAR_NAME}"

mv python/ ${TOOL_NAME} || die "renaming python to ${TOOL_NAME} failed"

tar zcvf ${TAR_NAME} ${TOOL_NAME}/ || die "tar balling ${TOOL_NAME} failed"

# if user arg "docker" is provided, copy tarball to host directory
if [[ "$1" == "docker" ]]; then
    echo "copying tar file to host directory"
    cp ${TAR_NAME} /host/ || die "copying tar file to host"
fi