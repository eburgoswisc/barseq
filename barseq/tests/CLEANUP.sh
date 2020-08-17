#!/usr/bin/env bash

DUMP="dump/"

# look for empty dir
if [ "$(ls -A $DUMP)" ]; then
     rm -r "./dump"
     mkdir dump
else
    exit 0
fi