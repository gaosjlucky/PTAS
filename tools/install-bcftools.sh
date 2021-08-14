#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

pushd "${SCRIPT_DIR}" > /dev/null

set -e
set -o pipefail

echo "Installing bcftools..."

INSTALL_DIR="${SCRIPT_DIR}/bcftools-1.9"
MAKE_DIR="${SCRIPT_DIR}/.."

rm -rf "${INSTALL_DIR}"
mkdir -p "${INSTALL_DIR}"

tar -xvf bcftools-1.9.tar.bz2 -C "${SCRIPT_DIR}"

pushd "${INSTALL_DIR}" > /dev/null

./configure --prefix=${MAKE_DIR}
make install

popd > /dev/null

echo "Finished installing bcftools."

popd > /dev/null
