for f in ../PgSA/toConsensus/*; do ./poa -read_fasta $f -pir $f.consensus -preserve_seqorder -do_progressive -hb -hbmin 0.001 blosum80.mat; done
