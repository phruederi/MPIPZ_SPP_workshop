#!/bin/bash

cat $1 | awk '/>/ {print $0, "barcodelabel="$1} !/>/ {print $0}' | \
    sed 's/=>/=/g;s/_[0-9]*$/;/g;s/ /;/g' \
    >> $2

