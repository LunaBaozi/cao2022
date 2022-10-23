#!/bin/bash

score_file=$(basename $1)_combined.sc

find $1 -name '*.sc' -exec cat {} \; 2>/dev/null 2>/dev/null | grep description | head -n1 > $score_file
find $1 -name '*.sc' -exec cat {} + | grep -v description | grep -v SEQUENCE >> $score_file

