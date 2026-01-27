## This folder contains code for generating Fig. 4.

If a figure requires both Python and R code, please run the Python scripts first to generate intermediate files, and then run the R scripts to generate the final figures.

Some of the required raw data files are provided in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data

Additional required files are provided within this folder.

### Extra notes

#### Fig4a:
After running the Fig4a.py, a FASTA file containing the selected somatic and progeny readouts will be generated.

The next step is to perform multiple sequence alignment and cell lineage tree construction using the following commands:

mafft --thread 50 ConsensusSeq_Plant1.fasta > ConsensusSeq_Multiple_Alignment_Plant2.fasta &

FastTree -nt ConsensusSeq_Multiple_Alignment_Plant2.fasta > ConsensusSeq_Tree_Plant2.nwk &

The resulting Newick tree file (.nwk) can be visualized using iTOL (Interactive Tree Of Life): https://itol.embl.de/

#### Fig4d:
After running the Fig4d.py, a FASTA file containing the selected somatic and progeny readouts will be generated.

The next step is to perform multiple sequence alignment and cell lineage tree construction using the following commands:

mafft --thread 50 ConsensusSeq_Plant3.fasta > ConsensusSeq_Multiple_Alignment_Plant3.fasta &

FastTree -nt ConsensusSeq_Multiple_Alignment_Plant3.fasta > ConsensusSeq_Tree_Plant3.nwk &

The resulting Newick tree file (.nwk) can be visualized using iTOL (Interactive Tree Of Life): https://itol.embl.de/
