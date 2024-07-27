#!/bin/bash

if [[ $1 == "-h" ]]; then
    echo "Make all scripts executable"
    exit 0
fi

find "$(dirname "$0")" -type f \( -iname \*.sh -o -iname \*.py \) -exec chmod +x {} \;
