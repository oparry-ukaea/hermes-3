#!/bin/env bash
REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )

spack env activate . -p

# build_dir="$REPO_ROOT/builds/release"
# spack build-env hermes-3%gcc cmake "$REPO_ROOT" -B "$build_dir" -DBOUT_DOWNLOAD_SUNDIALS=ON
# spack build-env hermes-3%gcc cmake --build "$build_dir" -j 8 

build_dir="$REPO_ROOT/builds/debug"
spack build-env hermes-3%gcc cmake "$REPO_ROOT" -DCMAKE_BUILD_TYPE=Debug -B "$build_dir" -DBOUT_DOWNLOAD_SUNDIALS=ON
spack build-env hermes-3%gcc cmake --build "$build_dir" -j 8 

# cd "$build_dir"
# ctest
# cd -