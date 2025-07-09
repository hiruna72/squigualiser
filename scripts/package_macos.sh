#!/bin/bash

die (){
    echo >&2 "$@"
    exit 1
}

# wget https://github.com/indygreg/python-build-standalone/releases/download/20230726/cpython-3.8.16+20230726-aarch64-apple-darwin-install_only.tar.gz || die "python wget failed"
# tar xf cpython-3.8.16+20230726-aarch64-apple-darwin-install_only.tar.gz || die "untar python failed"

# Download Python standalone build for macOS x86_64
wget https://github.com/indygreg/python-build-standalone/releases/download/20230726/cpython-3.8.16+20230726-x86_64-apple-darwin-install_only.tar.gz || die "python wget failed"
tar xf cpython-3.8.16+20230726-x86_64-apple-darwin-install_only.tar.gz || die "untar python failed"

python/bin/python3.8 -m venv squig-venv || die "create venv failed"
source squig-venv/bin/activate || die "sourcing venc failed"
pip3 install --upgrade pip || die "upgrade pip failed"
export HTSLIB_CONFIGURE_OPTIONS="--enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no"


if [[ "$1" == "test_pypi" ]]; then
    echo "Installing from Test PyPI"
    pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple squigualiser --pre || die "test.pip install squigualiser failed"
else
    echo "Installing from PyPI"
    pip install squigualiser --no-cache || die "pip install squigualiser failed"
fi

find ./ -name __pycache__ -type d | xargs rm -r || die "removing pycache failed"
mv squig-venv/bin/squigualiser python/bin/ || die "moving squigualiser to bin failed"
cp -r squig-venv/lib/python3.8/site-packages/* python/lib/python3.8/site-packages/ || die "copying site-packages failed"
sed -i '' "1s/.*/#\!\/usr\/bin\/env python3.8/" python/bin/squigualiser
git clone --depth 1 --branch main https://github.com/hiruna72/squigualiser.git  || die "Failed to clone squigualiser"

cp -r squigualiser/docs python || die "docs copy failed"
cp squigualiser/scripts/squigualiser python || die "script copy failed" 
cp squigualiser/LICENSE python || die "license copy failed"
cp squigualiser/README.md python || die "readme copy failed"

rm -rf squigualiser || die "remove cloned dir failed"

mv python/ squigualiser || die "renaming python to squigualiser failed"
tar zcvf squigualiser.tar.gz squigualiser/ || die "tar balling squigualiser failed"
