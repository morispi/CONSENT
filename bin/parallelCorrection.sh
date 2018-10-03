#!/bin/bash

parallel -j "$8" ./bin/correctOneRead.sh "$2" "$3" "$4" "$5" "$6" "$7" 1 "$1" :::: "$1"/ListAlignments >> "$9"
