#!/bin/bash


root_dir="$(dirname "$(realpath "$0")")"
cd "$root_dir"

log_file="build.log"
err_file="build.err"
filtered_err_file="build_filtered.err"

rm -f "$log_file" "$err_file" "$filtered_err_file"
make clean > /dev/null && make > "$log_file" 2>"$err_file"
grep -e "Duplicate C++ declaration" -e "^Declaration is " -v "$err_file" >| "$filtered_err_file"

if [ -s "$filtered_err_file" ]; then
    cat "$filtered_err_file"
    exit 1
else
    echo "No warnings or errors (filtering out breathe duplicates)."
    exit 0
fi