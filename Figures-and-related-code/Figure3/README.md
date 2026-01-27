## This folder contains code for generating Fig. 3.

If a figure requires both Python and R code, please run the Python scripts first to generate intermediate files, and then run the R scripts to generate the final figures.

Some of the required raw data files are provided in: Pipeline-of-lineage-tracing-in-Arabidopsis/Data

Additional required files are provided within this folder.

### Extra notes

#### Fig3a:
After running the Fig3a.py, a TXT file containing the progeny mutations from different branches will be generated.

The Venn diagram can be visualized using the online tool at: https://www.bic.ac.cn/test/venn/#/
#### Fig3b:
Fig3b consists of two parts. 

(1) generation of pseudo-germline readouts and identification of mutations shared between pairs of branches (Fig3b-1.py and Fig3b-1.R); 

(2) generation of pseudo-germline readouts and identification of mutations shared by all three branches(Fig3b-2.py and Fig3b-2.R).
