#!/bin/bash

if [[ $1 == "-h" ]]; then
    echo "Make all scripts executable"
    exit 0
fi

find "$(dirname "$0")" -iname "*.sh" -exec chmod +x {} \;
