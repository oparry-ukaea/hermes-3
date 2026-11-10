#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <EXECUTABLE_ABS_PATH>" >&2
    exit 1
fi

executable_path=$1

petsc_path=$(ldd "$executable_path" |
    awk '$1 ~ /^libpetsc[.]so/ && $2 == "=>" && $3 ~ /^\// { print $3; exit }')

if [[ -z "$petsc_path" ]]; then
    echo "PETSc library path not found" >&2
    exit 1
fi

strumpack_path=$(ldd "$petsc_path" |
    awk '$1 ~ /^libstrumpack[.]so/ && $2 == "=>" && $3 ~ /^\// { print $3; exit }')

if [[ -z "$strumpack_path" ]]; then
    echo "STRUMPACK library path not found for PETSc: $petsc_path" >&2
    exit 1
fi

if [[ "$strumpack_path" =~ /strumpack-[^/]+-([a-z0-9]{32})/ ]]; then
    strumpack_hash=${BASH_REMATCH[1]}
else
    echo "Spack hash not found in STRUMPACK library path: $strumpack_path" >&2
    exit 1
fi

printf 'Linked strumpack version:\n'
spack find -v "strumpack/$strumpack_hash"
