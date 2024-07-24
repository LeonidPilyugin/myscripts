#!/bin/bash

if [[ $1 == "-h" ]]; then
    echo "Deletes all subdirs named 'trajectory' in '$1'"
    exit 0
fi

if [[ ! -d $1 ]]; then
    echo Directory "$1" does not exist
    exit 1
fi

for directory in $(find "$1" -type d); do
	if [[ $(basename "$directory") == "trajectory" ]]; then
		rm -rf "$directory"
	fi
done
