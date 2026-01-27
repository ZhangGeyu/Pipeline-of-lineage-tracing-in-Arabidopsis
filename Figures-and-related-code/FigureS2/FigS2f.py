from Bio import Phylo    # Import Phylo module for phylogenetic tree analysis

# Define function: filter phylogenetic tree nodes by branch length
def filter_by_branch_length(nwk_file, threshold=0):
    tree = Phylo.read(nwk_file, "newick")
    filtered = []
    for clade in tree.get_nonterminals():
        # skip the nodes without bootstrap value
        if clade.confidence is None:
            continue
        # extract branch length
        branch_length = clade.branch_length if clade.branch_length else 0
        # remain nodes which branch length > 0
        if branch_length > threshold:
            filtered.append(clade.confidence)
    return filtered

# Raw file: Tree constructed in Fig2a
bootstrap_list = filter_by_branch_length("ConsensusSeq_Tree_Plant1.nwk",threshold=0)
bootstrap_df = pd.DataFrame({'bootstrap':bootstrap_list})

#Clean file
bootstrap_df.to_csv('Tree_bootstrap.txt', sep = '\t', index = False)
