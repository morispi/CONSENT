#!/bin/bash

minSupport="$1"
windowSize="$2"
merSize="$3"
commonKMers="$4"
solid="$5"
windowOverlap="$6"
tmpdir="$7"
LRSCf="$8"
cluster="$9"

#Retrieve the cluster's sequences
$LRSCf/bin/retrieveClusterSequences.py $tmpdir/Clusters/$cluster $tmpdir/RawLongReads/ $tmpdir/ClustersSequences/$cluster"_ref" $tmpdir/ClustersSequences/$cluster"_reads"

#Align the cluster's sequences
$LRSCf/minimap2/minimap2 -w5 --dual=yes -P --no-long-join -f10000000 $tmpdir/ClustersSequences/$cluster"_reads" $tmpdir/ClustersSequences/$cluster"_ref" > $tmpdir/ClustersAlignments/"dual_$cluster"

#Sort the cluster's alignments
sort -T $tmpdir/ -k1,1 -k3,3n -k4,4n -k8,8n -k9,9n $tmpdir/ClustersAlignments/"dual_$cluster" > $tmpdir/ClustersAlignments/$cluster

#Run the correction process with the aligned cluster
$LRSCf/bin/LRSelfCorrection -a $tmpdir/ClustersAlignments/$cluster -d $tmpdir/RawLongReads/ -s $minSupport -l $windowSize -k $merSize -c $commonKMers -f $solid -m $windowOverlap -j 1

#Remove temporary files
rm $tmpdir/ClustersSequences/$cluster"_ref" $tmpdir/ClustersSequences/$cluster"_reads" $tmpdir/ClustersAlignments/"dual_$cluster" $tmpdir/ClustersAlignments/$cluster
