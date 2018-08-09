#!/bin/bash

#while read line; do
	#./bin/GoRExtractor -a "$2/""$line" -d TEST12/RawLongReads/ -s 3 -l 100 -k 7 -c 3 -f 3 -m 50 -j 1 >> gr.fa
	#./bin/GoRExtractor -a "$2/""$line" -d TEST12/RawLongReads/ -s 3 -l 300 -k 9 -c 7 -f 3 -m 33 -j 1 >> gr.fa
	#./bin/GoRExtractor -a "$2/""$line" -d TEST12/RawLongReads/ -s 3 -l 500 -k 8 -c 7 -f 3 -m 17 -j 1 >> gr.fa
#	./bin/GoRExtractor -a "$2/""$line" -d TEST12/RawLongReads/ -s 3 -l 50 -k 7 -c 3 -f 3 -m 10 -j 1 >> ALLREADS.fasta
#	echo "$1"
	#parallel -j "$3" ./correctOneRead.sh 5 50 6 3 5 10 1 "$2" :::: "$1" >> ALLREADS.fasta #92% id
	#parallel -j "$3" ./correctOneRead.sh 5 50 7 2 5 10 1 "$2" :::: "$1" >> ALLREADS.fasta #93% id
	#echo "$1 $2 $3 $4 $5 $6 $7 $8 $9"
	parallel -j "$8" ./bin/correctOneRead.sh "$2" "$3" "$4" "$5" "$6" "$7" 1 "$1" :::: "$1"/ListAlignments >> "$9"
	#parallel -j "$3" ./correctOneRead.sh 5 75 7 3 10 10 1 "$2" :::: "$1" >> ALLREADS.fasta
#done < "$1"
