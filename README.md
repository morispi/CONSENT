# CONSENT

CONSENT (sCalable self-cOrrectioN of long reads with multiple SEquence alignmeNT) is a self-correction method for long reads.
It works by, first, computing overlaps between the long reads, in order to define an alignment pile (i.e. a set of overlapping reads used for
correction) for each read. Each read's alignment pile is then further divided into smaller windows, that are corrected idependently.
First, a multiple alignment strategy is used in order to compute consensus. Then, this consensus is further polished with a local de Bruijn
graph, in order to get rid of the remaining errors.
Additionally to error correction, CONSENT can also perform contigs polishing.

Requirments
--------------

  - A Linux based operating system.
  - Python3.
  - A g++ version supporting C++11
  - fpa (at least v0.4) accessible through your PATH environment variable (https://github.com/natir/fpa)
  
Installation
--------------

Clone the CONSENT repository, along with its submodules with:

  ```bash
  git clone --recursive https://github.com/morispi/CONSENT
  ```

Then run the install.sh script:

  ```bash
  ./install.sh
  ```
  
Running CONSENT
--------------

### Self-correction

To run CONSENT for long reads self-correction, run the following command:

`./CONSENT-correct --in longReads.fasta --out result.fasta --type readsTechnology`

  - longReads.fasta:	fasta file of long reads to correct, with one sequence per line.
  - result.fasta:		fasta file where to output the corrected long reads.
  - readsTechnology:	Indicate whether the long reads are from PacBio (--type PB) or Oxford Nanopore (--type ONT)


### Polishing

To run CONSENT for contigs polishing, run the followning command:

`./CONSENT-polish --contigs contigs.fasta --reads longReads.fasta --out result.fasta`

  - contigs.fasta:		fasta file of contigs to polish, with one sequence per line.
  - longReads.fasta:	fasta file of long reads to use for polishing, with one sequence per line.
  - result.fasta:		fasta file where to output the polished contigs.

### Options

      --windowSize INT, -l INT:      Size of the windows to process. (default: 500)
      --minSupport INT, -s INT:      Minimum support to consider a window for correction. (default: 4)
      --maxSupport INT, -S INT:      Maximum number of sequences to include into a window. (default: 1,000)
      --maxMSA INT, -M:              Maximum number of sequences to include into the MSA. (default: 150)
      --merSize INT, -k INT:         k-mer size for chaining and polishing. (default: 9)
      --solid INT, -f INT:           Minimum number of occurrences to consider a k-mer as solid during polishing. (default: 4)
      --anchorSupport INT, -c INT:   Minimum number of sequences supporting (Ai) - (Ai+1) to keep the two anchors in the chaining. (default: 8)
      --minAnchors INT, -a INT:      Minimum number of anchors in a window to allow consensus computation. (default: 2)
      --windowOverlap INT, -o INT:   Overlap size between consecutive windows. (default: 50)
      --nproc INT, -j INT:           Number of processes to run in parallel (default: number of cores).
      --minimapIndex INT, -m INT:    Split minimap2 index every INT input bases (default: 500M).
      --tmpdir STRING, -t STRING:    Path where to store the temporary overlaps file (default: working directory, as Alignments_dateTimeStamp.paf).
      --help, -h:                    Print this help message.

Notes
--------------

CONSENT has been developed and tested on x86-64 GNU/Linux.          
Support for any other platform has not been tested.

Authors
--------------

Pierre Morisse, Camille Marchet, Antoine Limasset, Arnaud Lefebvre and Thierry Lecroq.

Reference
--------------

Pierre Morisse, Camille Marchet, Antoine Limasset, Thierry Lecroq, Arnaud Lefebvre; CONSENT: Scalable self-correction of long reads with multiple sequence alignment, bioRxiv 546630; doi: https://doi.org/10.1101/546630 

Contact
--------------

You can report problems and bugs to pierre[dot]morisse2[at]univ-rouen[dot]fr
