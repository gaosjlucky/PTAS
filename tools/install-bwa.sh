#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

pushd "${SCRIPT_DIR}" > /dev/null

set -e
set -o pipefail

echo "Installing bwa..."

INSTALL_DIR="${SCRIPT_DIR}/bwa-0.7.17"
MAKE_DIR="${SCRIPT_DIR}/.."

rm -rf "${INSTALL_DIR}"
mkdir -p "${INSTALL_DIR}"

tar -xvf bwa-0.7.17.tar.bz2 -C "${SCRIPT_DIR}"

pushd "${INSTALL_DIR}" > /dev/null

make
ln -s ${INSTALL_DIR}/bwa ${MAKE_DIR}/bin/

popd > /dev/null

echo "Finished installing bwa."

popd > /dev/null
