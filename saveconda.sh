#!/bin/bash

DIRECTORY=~/conda-dumps

mkdir -p $DIRECTORY

for env in $(conda env list | tail -n +3 | cut -d' ' -f2- | xargs); do
    conda env export -n "$env" > "$DIRECTORY/$(basename ${env}).yaml"
done
