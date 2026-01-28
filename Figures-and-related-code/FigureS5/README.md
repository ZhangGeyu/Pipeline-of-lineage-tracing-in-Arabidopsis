## This folder contains code for generating Fig. S5.

If a figure requires both Python and R code, please run the Python scripts first to generate intermediate files, and then run the R scripts to generate the final figures.

Some of the required raw data files are provided in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data

Additional required files are provided within this folder.

### Extra notes

#### FigS5a, FigS5c, and FigS5e:
After running FigS5a.py, FigS5c.py, and FigS5e.py, a file named "character_matrix.csv" will be generated, containing the mutation information of both somatic and progeny readouts. 

In addition, "mutation_prior.csv" will be generated, which shows the rate of each mutation type at different positions.

The next step is to perform cell lineage tree construction using the following commands:

python ../startle/scripts/nj.py character_matrix.csv --output tree.newick &

startle large character_matrix.csv mutation_prior.csv tree.newick --output startle &

The resulting Newick tree file (.nwk) can be visualized using iTOL (Interactive Tree Of Life): https://itol.embl.de/
