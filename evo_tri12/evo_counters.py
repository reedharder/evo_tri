#evo_counters.py: these are functions called by hit_plotter2D() and hit_plotter3d() function

#they are stripped down evolutionary triangulation and gene finding functions,
#repurposed specifically for counting number of SNPs or genes found for a given threshold configuration


import numpy as np

#parameters are equivalent to evo_triangulator, except genenum: a boolean, if True (if genes will also be found), also return a list of SNPs found for gene finding
def stripped_down_evo_triangulator(fst_dict, pops, operators, thresholds, genenum):

    #key population combos to corresponding operator/threshold indices
    pair_key={0: pops[0] + '_' + pops[1], 1: pops[0] + '_' + pops[2], 2: pops[1] + '_' + pops[2]}

    #get arrays of SNPS for each population pair, reduced by Fst thresholds and selected operator
    reduced_snp_arrays_list=[]
    for i in xrange(0,3):
        current_func=operators[i]
        #below, snp/fst array for relevant population pair has rows filtered by Fst threshold and specified operator (current_func()),
        #and column of remaining snps is returned
        with np.errstate(invalid='ignore'): #ignore warnings when comparing divide-by-zero Fsts, these will not return SNP hits
            reduced_snp_arrays_list.append(fst_dict[pair_key[i]][current_func(fst_dict[pair_key[i]][:,1],thresholds[i]),0])

    #get set of overlapping snps between reduced snp sets
    reduced_overlap = reduce(lambda x,y: np.intersect1d(x,y,True), reduced_snp_arrays_list)

    #return list of snps if number of genes in vicinity will be calculated
    if genenum == True:
        snplist=reduced_overlap
    else:
        snplist==[]
    #return number of overlapping SNPs found, and list of SNPs if requested
    return [len(reduced_overlap), snplist]



#parameters are equivalent to genefinder(). Returns number of genes found for a given SNP list
def stripped_down_genefinder(current_overlap_snps, loc_dict, bp_range, genenames, genenums):

    bp_range=abs(bp_range)

    #select arbitrary pop to get rs numbers and corresponding positions
    locs=loc_dict[loc_dict.keys()[0]][:,[0,1,2]]

    #filter location data such that only locations of overlapping snps remain
    locs=locs[np.searchsorted(locs[:,0],current_overlap_snps)]

    #array of location of lower and upper bound to count as a "hit"
    gene_bounds=np.column_stack((genenums[:,0],genenums[:,1]-bp_range,genenums[:,2]+bp_range))

    #find genes...
    gene_boolean=np.squeeze(np.zeros((len(genenums[:,0]),1),dtype=bool)) #initialize boolean array for genes found
    n = len(locs[:,0]) #get number of snps to find genes in vicinity of
    for i in xrange(0,n):
        #get boolean array of genes for which share a chromosome with current snp, and current snp falls within
        new_genes=np.squeeze(reduce(lambda x,y: np.logical_and(x,y),[(gene_bounds[:,0]==locs[i,1]),(gene_bounds[:,1]<=locs[i,2]),(gene_bounds[:,2]>=locs[i,2])]))
        #if new genes are found, update gene boolean
        if new_genes.any()==True:
            gene_boolean=np.logical_or(new_genes,gene_boolean)

    #force gene_boolean to be 1d
    gene_boolean=np.squeeze(gene_boolean)

    #get array of genes found
    genes_found=genenames[gene_boolean]

    #return number of genes found
    return len(genes_found)