
#Main function file: contains primary functions for algorithm data flow, from downloading or loading files, to parsing files and calculating Fsts,
#to running evolutionary triangulation algorithm, to finding genes in vicinty of SNPs found by evolutionary triangulation algorithm,
#to graphing number of SNPs and genes found over ranges of Fst cutoffs.

#This package requires numpy (and matplotlib for hit_plotter2D() and hitplotter3D()

from os import listdir, getcwd, chdir, path
from glob import glob
from string import split, count
import numpy as np

'''
evo_hapmapDLer: function for externally downloading hapmap allele_freq files for an arbitrary number of populations

    parameters:
        phase: String that specifies hapmap phase, e.g. '2009-01_phaseIII'
        pops: list of populations as 3 letter hapmap abbreviations, e.g. 'CEU', or 'YRI'
        chrs: Enter string of desired chromosomes (1-22,X,Y,M) to be downloaded for all populations seperated by commas, no spaces.
            Use "autosome" to select all autosomal chromosomes (e.g. "autosome,X,Y").
            Use colons to indicate ranges (e.g. "1,4:10,15,Y")
        file_folder: output folder pathname in which to dump downloaded files.

'''

def evo_hapmapDLer(phase='2009-01_phaseIII', pops=['CEU','YRI','GIH'], chrs='autosome,X',
    file_folder="current"):

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    #get parsing and downloading functions
    from chromosome_parser import chromeparse
    from hapmap_tools import get_hapmap

    #get list of chromosomes from user input
    chr_list=chromeparse(chrs)[1]

    #download hapmap files for specifed populations and chromosomes
    get_hapmap(chr_list, phase, pops, file_folder)

    return None


'''
evo_text2numpy: converts custom data in text format into a .npz (numpy archive) file, for efficient processing

parameters:
    pop1,pop2,pop3: hapmap population 3 letter abbrev, or custom population abbreviation.
    customFst: True if custom FST data will be provided
    no_popdata: if True, no population data file will be created
    pop[1-3]filenames: text file name strings for each population to be combined in numpy array archive. Rows  of text file should be each SNP, columns
        should be rs number, chromosome number, snp position, allele frequency, and sample size, in that order, seperated by white space.
        If certain data not required for calculation (for example, if custom FST values will be provided), fill in place holder column with -1's.
        SNPs given do not need to be completely overlapping.
    pop[12-13-23]FSTfilenames: text file name strings for each population to be combined in numpy array archive. Rows should be each SNP,
        columns should be rs number and associated Fst, seperated by white space. SNPs given do not need to be completely overlapping.
    file_folder: input and output folder pathname.
    out_compressed: output file will be compressed
    pop_out: file name (without .npz extension) of population raw data
    FST_out: file name (without .npz extension) of population FST data
'''

def evo_text2numpy(pop1="CEU", pop2="YRI", pop3="GIH", customFst=False, no_popdata=False,
    pop1filename='pop1textdata.txt', pop2filename='pop2textdata.txt', pop3filename='pop3textdata.txt',
    pop12filename='pop12textdata.txt',pop13filename='pop13textdata.txt', pop23filename='pop23textdata.txt',
    file_folder="current", out_compressed=False, pop_out="pop_archive1", FST_out="FST_archive1"):


    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)

    import numpy as np

    if no_popdata==False:
        #save population raw data
        x=np.loadtxt(pop1filename)
        y=np.loadtxt(pop2filename)
        z=np.loadtxt(pop3filename)
        dict_args={pop1:x,pop2:y,pop3:z}

        if out_compressed == False:
            np.savez(pop_out,**dict_args)
        else:
            np.savez_compressed(pop_out,**dict_args)


    #if specified, also save custom FST data in the same manner
    if customFst==True:
        x=np.loadtxt(pop12filename)
        y=np.loadtxt(pop13filename)
        z=np.loadtxt(pop23filename)
        dict_args={pop1+'_'+pop2:x,pop1+'_'+pop3:y,pop2+'_'+pop3:z}

        if out_compressed == False:
            np.savez(FST_out,**dict_args)
        else:
            np.savez_compressed(FST_out,**dict_args)




'''evo_tri_data: accesses and sets up data, calculates pairwise Fst values for overlaping SNPs. Saves a dictionary of 2d numpy arrays
   for each population combo with columns as rs numbers and corresponding Fst for each SNP common to all three populations. Also optionally saves
   a dictionary of 2d numpy arrays for each population with columns as rs number, chromosome number (X->23,Y->24,M->25), snp position, allele frequency,
   and sample size, in that order, for each SNP common to all three populations. See parameters for details.

parameters:
    datasource: 1 is hapmap online (will prompt for local DB location and hapmap phase), 2 is local hapmap txt allele_freq text files,
        3 is local custom data (see below for required format), with or without Fst and gene location data.

    data3customFst: Only relevant if data source 3 is selected. If True, will use custom calculated Fst values, see below.
                    else, will calculate Fst from given data. Only relevant if data source 3 is selected.
    no_popdata: Only relevant if data source 3 is selected and custom Fst will be used. If true, only SNP/precalculated Fst is taken, allowing for faster
                calculation. Genefinding will not be possible if this option is used.
    custom_popdata, custom_FSTdata: Only relevant if data source 3 is selected. File names for preprocessed custom data in .npz files.
                    To convert text data files to .npz, see evo_text2numpy()
                    custom_popdata file should be .npz dicionary of 3 numpy matrices of floats, one for each population, keyed with '[population abbreviation]'.
                        rows should be SNPs, columns should be rs number, chromosome number, snp position, allele frequency, and sample size, in that order.
                        Chromosomes X,Y and mitochondrial (M) should be coded as 23, 24 and 25 respectively.
                        If certain data  are not required for calculation (for example, if custom FST values will be provided), fill in place holder column with -1's.
                        SNPs given do not need to be completely overlapping.
                    If it is being used, custom_FSTdata file should be be .npz dicionary of 3 numpy matrices of floats, one for each population combo (1-2, 1-3, 2-3),
                    keyed with '[population 1 abbreviation] + '[population 2 abbreviation]'. Rows should be each SNP, columns should be rs number and
                    associated Fst, SNPs given do not need to be completely overlapping.

                    Data required:
                    If custom Fst values are provided, only Fsts and matching rs-number file are required
                    If unbiased Fst will be calculated, at least rs-numbers, allele frequencies, and sample sizes for each population are required.
                    If uncorrected Fst will be calculated, only rs numbers and allele frequencies are required.
                    If genefinding will be used, chromosome and snp location data for each population is required/
                    Each .txt file should be a column of numerical values corresponding to the order of rs numbers for that population provided.
    pop1,pop2,pop3: hapmap population 3 letter abbrev, or custom population abbreviation.
    chrs: Only relevant with datasource 1 and datasource 2. Enter string of desired chromosomes (1-22,X,Y,M) seperated by commas, no spaces. Use "autosome" to select all autosomal chromosomes (e.g. "autosome,X,Y"). Use colons to indicate ranges (e.g. "1,4:10,15,Y")
    phase: only relevant if hapmap files will be downloaded (i.e. datasource 1). String that specifies hapmap phase, e.g. '2009-01_phaseIII'
    pop_out, FST_out: string with output file name (.npz extension not necessary). For raw population data and calculated Fst data respectively.
    out_compressed: if true, output files are compressed
    file_folder: directory with files, for datasource 2 and 3; directory for storing downloaded file database for datasource 1
'''


def evo_tri_data(datasource=2,
    data3customFst=False, no_popdata=False, custom_popdata='pop_archive1.npz', custom_FSTdata='FST_archive1.npz',
    pop1="CEU", pop2="YRI", pop3="GIH", chrs="autosome,X", unbiasedFst=True,
    phase='2009-01_phaseIII', pop_out="pop_data1", FST_out="fst_data1", out_compressed=False, file_folder="current"):

    import numpy as np

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    # list of populations
    pops =[pop1,pop2,pop3]

    #parse chromosome input (see chromosome_parser module for details)
    from chromosome_parser import chromeparse
    chr_regexes = chromeparse(chrs)[0]

    if data3customFst==True: #i.e., if custom Fst values be used...
    #return arrays that include only SNPs shared between the three population pairs
        chdir(file_folder)
        print "Formating custom Fst data..."
        #load  custom data file
        Fst_data=np.load(custom_FSTdata)
        for i in Fst_data:
            #sort array by first column (i.e. sort by rs number)
            Fst_data[i]=Fst_data[i][Fst_data[i][:,0].argsort()]

        #filter arrays such that only SNPs common to all three populations remain
            #get array of common SNPs
        common_snps = reduce(lambda x,y: np.intersect1d(x,y,True),(i[:,0] for i in (Fst_data[pop1+'_'+pop2],Fst_data[pop1+'_'+pop3],Fst_data[pop2+'_'+pop3])))
            #use common SNP array to get index of these SNPs in each Fst data array and appropriately filter it
        for i in Fst_data:
            Fst_data[i] = Fst_data[i][np.searchsorted(Fst_data[i][:,0],common_snps)]

        #save population data and fst data as .npz file, compressed or uncompressed as specified
        if out_compressed == False:
            np.savez(FST_out,**Fst_data)
        else:
            np.savez_compressed(FST_out,**Fst_data)





    #exit if custom Fst data will suffice and no pop data will be used
    if (no_popdata == True) and (data3customFst == True):

        return None

    else: # i.e. if population data will be used...

        #if data needs to be downloaded from hapmap,
        if datasource==1:
            chr_groups = chromeparse(chrs)[1] #get list of chromosomes to be retrieved
            from hapmap_tools import get_hapmap
            get_hapmap(chr_groups,phase,pops,file_folder)



        # set up data from local hapmap allele .txt or .txt.gz files
        if datasource==1 or datasource==2:
            import re #regular expressions tools
            import numpy as np
            import gzip #for .gz files
            from hapmap_tools import allelefreq_readline, chr_as_float #function to parse hapmap allele frequency files, and function to code chromosomes numerically

            #go to specified directory
            chdir(file_folder)
            #initialize lists for files for each population
                #for each population, lists are: rs ids, chromosomes, positions, allele frequencies, and sample sizes, respectively
            pop_data={pop1:[], pop2:[], pop3:[]}

            print "Parsing HapMap files..."
            #grab relevant file names for each population, extract relevant data from each file
            for i in pop_data:
                #get list of relevant file names for given population
                glob_string='allele_freqs_chr*_%s_*.txt*' % i
                allele_freq_files = glob(glob_string)
                #select from list those files for appropriate chromosomes
                allele_freq_files_chrs = [x for x in allele_freq_files if any( re.match(pattern, x) for pattern in chr_regexes)] #TEST, THEN REMOVE _chr VARIABLE
                #get data line by line for each of these file

                for filename in allele_freq_files_chrs:
                    #get chromosome from file name, to avoid reading and converting it from each line; see chr_as_float() defined above
                    print "parsing %s" % filename
                    chr_float=chr_as_float(filename)
                    #open file and read...
                    try:
                        with gzip.open(filename) as freq_file: #assume hapmap file is .gz
                            next(freq_file) #skip header
                            for line in freq_file:
                                #parse line with allelefreq_readline()
                                #append returned 1d numpy array to list of 1d arrays for current population
                                pop_data[i].append(allelefreq_readline(line,chr_float))
                    except IOError: #if IOError, take as uncompressed text file
                        with open(filename) as freq_file:
                            next(freq_file) #skip header
                            for line in freq_file:
                                #parse line with allelefreq_readline()
                                #append returned 1d numpy array to list of 1d arrays for current population
                                pop_data[i].append(allelefreq_readline(line,chr_float))

                #convert list of 1d arrays to 2d numpy array, for each population
                pop_data[i]=np.array(pop_data[i])
                #sort array by first column (i.e. sort by rs number)
                pop_data[i]=pop_data[i][pop_data[i][:,0].argsort()]


        #set up data from local preprocessed custom data
        if datasource == 3:
            chdir(file_folder)
            #load  custom data file
            print "Loading custom population data..."
            np.load(custom_popdata)
            for i in pop_data:
                #sort array by first column (i.e. sort by rs number)
                pop_data[i]=pop_data[i][pop_data[i][:,0].argsort()]



        #filter arrays such that only SNPs common to all three populations remain
        #get array of common SNPs
        print "Filtering SNPs..."
        common_snps = reduce(lambda x,y: np.intersect1d(x,y,True),(i[:,0] for i in (pop_data[pop1],pop_data[pop2],pop_data[pop3])))
        #use common SNP array to index of these SNPs in each population data array and appropriately filter it
        for i in pop_data:
            pop_data[i] = pop_data[i][np.searchsorted(pop_data[i][:,0],common_snps)]

        #save population data and fst data as .npz file, compressed or uncompressed as specified
        if out_compressed == False:
            np.savez(pop_out,**pop_data)
        else:
            np.savez_compressed(pop_out,**pop_data)


        # if Fsts need to be calculated...
        if data3customFst == False:
            #Calculate pairwise Fst between each population pair, then save
            import fst_calc #import fst functions

            print "Calculating Fsts..."

            if unbiasedFst == True: #based on user parameters, select FST formula to use
                FST = fst_calc.CaiWeirFst2pop #use unbiased FST
                Fstpop1pop2 = FST(pop_data[pop1][:,4],pop_data[pop2][:,4],pop_data[pop1][:,3],pop_data[pop2][:,3]) #input sample size and frequency vectors for pops 1 and 2, calculate pairwise FST
                Fstpop1pop3 = FST(pop_data[pop1][:,4],pop_data[pop3][:,4],pop_data[pop1][:,3],pop_data[pop3][:,3]) # repeat for other population pairs
                Fstpop2pop3 = FST(pop_data[pop2][:,4],pop_data[pop3][:,4],pop_data[pop2][:,3],pop_data[pop3][:,3])
            else:
                FST = fst_calc.WrightFst #else, use traditional formulation
                Fstpop1pop2 = FST(pop_data[pop1][:,3],pop_data[pop2][:,3]) #input sample size and frequency vectors for pops 1 and 2, calculate pairwise FST
                Fstpop1pop3 = FST(pop_data[pop1][:,3],pop_data[pop3][:,3]) # repeat for other population pairs
                Fstpop2pop3 = FST(pop_data[pop2][:,3],pop_data[pop3][:,3])
            #build dictionary with numpy arrays of snps and Fsts, for each population combo
            snplist=pop_data[pop1][:,0] #get ordered list of rs numbers
            Fst_data = {pop1+'_'+pop2:np.column_stack((snplist,Fstpop1pop2)), pop1+'_'+pop3:np.column_stack((snplist,Fstpop1pop3)), pop2+'_'+pop3:np.column_stack((snplist,Fstpop2pop3))} #Fst dictionary

            #save population data and fst data as .npz file, compressed or uncompressed as specified
            if out_compressed == False:
                np.savez(FST_out,**Fst_data)
            else:
                np.savez_compressed(FST_out,**Fst_data)

    print "Done"

    return None


'''evo_triangulator(): implements evolutionary triangulation algorithm, saving SNPs found to .npy and (optionally) text files.

parameters:
    fst_file: name of dictionary containing Fst data for pairings of three populations (such as Fst_out from evo_tri_data() )
    pops: list of population abbreviations for population 1, population 2, and population 3, in that order
    ptile_cutoff: if True, will consider cuttoffs entered below as percentiles (in decimal form, e.g. 95% = .95) rather than as absolute cuttoffs. Default to False.
    pop12lim,pop13lim,pop23lim: Fst threshold (string with operator [>,>,>=,<=] followed by Fst between 0 and 1) for pop1-pop2 Fsts, pop1-pop3 Fsts, and pop2-pop3 Fsts respectively
    snps2screen: true prints snps found to screen
    snps2txt: true prints snps to txt file
    snps_filename: name for snps txt file and .npy file
    file_folder: directory from which to load and save files

    '''

def evo_triangulator(fst_file='fst_data1.npz', pops=['CEU','YRI','GIH'], ptile_cutoff=False, pop12lim=">=.45", pop13lim=">=.45", pop23lim="<=.05", snps2screen=True, snps2txt=False, snps_filename="evotri_snps",
    file_folder="current"):

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)
    import numpy as np

    #load fst data
    fst_dict=np.load(fst_file)

    lims = [pop12lim,pop13lim,pop23lim]
    op_inputs = [] #list of operators
    fst_bounds = [] #and list of corresponding Fst values
    #parse bounds inputs
    for lim in lims:
        #make sure operator is valid
        if lim[0] != ">" and lim[0] != "<":
            print "invalid operator in given bounds"
            import sys
            sys.exit()
        #seperate operator from Fst (Fst should start at second or third char of input string
        if lim[1] == '=':
            op=lim[:2]
            fst_bound=float(lim[2:])
        else:
            op=lim[:1]
            fst_bound=float(lim[1:])
        #make list of operators and thresholds
        op_inputs.append(op)
        fst_bounds.append(fst_bound)

    #key population combos to corresponding operator/threshold indices
    pair_key={0: pops[0] + '_' + pops[1], 1: pops[0] + '_' + pops[2], 2: pops[1] + '_' + pops[2]}

    #if cutoffs are to be percentiles, get corresponding absolute cutoffs for each
    if ptile_cutoff==True:
        for i in xrange(0,3):
            fst_bounds[i]=np.percentile(fst_dict[pair_key[0]][:,1],100*fst_bounds[i])

    #dictionary of possible operators
    import operator #
    ops = {">": operator.gt,
         "<": operator.lt,
         ">=": operator.ge,
         "<=": operator.le}

    #get three operators as list of functions
    op_func=[]
    for i in xrange(0,3):
        op_func.append(ops[op_inputs[i]])


    #get arrays of SNPS for each population pair, reduced by Fst thresholds and selected operator
    reduced_snp_arrays_list=[]
    for i in xrange(0,3):
        current_func=op_func[i]
        #below, snp/fst array for relevant population pair has rows filtered by Fst threshold and specified operator (current_func()),
        #and column of remaining snps is returned
        with np.errstate(invalid='ignore'): #ignore warnings regarding comparison to NaN values, these will return false and thus not return SNPs
            reduced_snp_arrays_list.append(fst_dict[pair_key[i]][current_func(fst_dict[pair_key[i]][:,1],fst_bounds[i]),0])

    #get set of overlapping snps between reduced snp sets
    reduced_overlap = reduce(lambda x,y: np.intersect1d(x,y,True), reduced_snp_arrays_list)
    #get overlap size
    overlap_shape=len(reduced_overlap)

    #print number of overlapping SNPs found to screen
    print '%s overlapping SNPs found.' % (str(overlap_shape))

    #save overlapping SNPs as .npy file
    np.save(snps_filename, reduced_overlap)

    #save overlapping SNPs as text file if requested
    if snps2txt == True:
        fn=snps_filename + '.txt'
        np.savetxt(fn,reduced_overlap,fmt='%.d')

    #print overlapping SNPs to screen if requested
    if snps2screen == True:
        print "overlapping SNPs:"
        for i in reduced_overlap:
            print int(i)

    return None

''' gene_finder: finds genes in regions of overlapping SNPs found (within specified range)

parameters:
    snplist: .npy file storing array of overlapping SNPs
    custom_loc: if True, uses .npy file storing 3 column 2d numpy array, with rs numbers, corresponding chromosomes, and corresponding locations, in that order, seperated
    by white space. Allows use of custom location file, not generated by evo_tri_data.
    loc_data: string with .npy or .npz file name (without extension), with file containing SNP location data. If custom_loc==True, provide 3 column file (see above).
        Otherwise, use .npz file in format of pop_out from evo_tri_data().
    pop: only relevant if custom_loc is false. Specify popuation in loc_data from which to take SNP locations. If population is not given or not found,
        popuation data will be taken from first population in loc_data.
    bp_range: integer specifing how many bases from a SNP a gene must be to be considered a hit.
    custom_gene: if True, use custom gene text file. File should contain 4 columns seperated by whitespace: chromosome number, gene start, gene end, gene name.
        Chromosomes X,Y and mitochondrial (M) should be coded as 23, 24 and 25 respectively.
        If False, will use default genes on file.
    genes2text: if True, will save text file of genes found
    genes2screen: if True, will print genes found to screen
    display_chr: if True, will display and save chromosome number with genes
    genes_filename: name for .npy and .txt (if requested) file, where gene data is saved.
    file_folder: folder from which to save and load files. If using default gene list, genefinder_default_genelocs.npy and genefinder_default_genenames.npy must be in this file

'''

def genefind_local(snplist='evotri_snps.npy', custom_loc=False, loc_data='pop_data1', pop='first', bp_range=100000, custom_gene=False, custom_genefile='custom_genes.txt',
    genes2txt=False, genes2screen=True, genes_filename='genes1', display_chr=False, file_folder="current"):

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)
    import numpy as np

    bp_range=abs(bp_range)

    #get positions of possible SNPs in a certain population
    if custom_loc == True: # simply load if 3 column custom positions array is provided
        fname = loc_data + '.npy'
        locs = np.load(fname)
    else: #otherwise, take positions from a pop_data-stayle array such as that from evo_tri_data()
        fname = loc_data + '.npz'
        loc_dict = np.load(fname)
        try:
            #try getting rs numbers and corresponding positions from requested pop...
            locs=loc_dict[pop][:,[0,1,2]]
        except KeyError:
            #if pop not given or found, select arbitrary pop to get rs numbers and corresponding positions
            locs=loc_dict[loc_dict.keys()[0]][:,[0,1,2]]

    if custom_gene==True: #load custom gene file if requested
        genenums=np.loadtxt(custom_genefile, dtype=int, usecols=(0,1,2)) #get genes' chromosome, start and end as integers
        genenames=np.loadtxt(custom_genefile, dtype=str, usecols=(3)) #get genes' names as strings
    else: #otherwise use default data
        import evo_tri
        chdir(path.dirname(evo_tri.__file__))
        genenums=np.load('genefinder_default_genelocs.npy')
        genenames=np.load('genefinder_default_genenames.npy')
        chdir(file_folder)

    #load list of overlapping SNPs
    overlap_snps=np.load(snplist)

    #filter location data such that only locations of overlapping snps remain
    locs=locs[np.searchsorted(locs[:,0],overlap_snps)]

    #array of chromosome location of lower and upper bound to count as a "hit" for each gene
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

    #if user does not want to save chromosome number, print and save information without chromosome number
    if display_chr == False:
        #get array of genes found
        genes_found=genenames[gene_boolean]
        num_genes=len(genes_found)

        #print number of genes found to screen
        print '%s genes found.' % (num_genes)

        #if genes have been found, save and display according to user specifications
        if (num_genes > 0):
            #save gene names as .npy file
            np.save(genes_filename, genes_found)

            #save gene names as text file if requested
            if genes2txt == True:
                fn=genes_filename + '.txt'
                np.savetxt(fn,genes_found, fmt='%.s')

            #print overlapping SNPs to screen if requested
            if genes2screen == True:
                print "genes found:"
                for i in genes_found:
                    print i

    #otherwise, print and save genes with chromosome number
    else:
        #get array of genes found
        genes_found=genenames[gene_boolean]

        #get array of chromosomes corresponding to above as strings
        chrs_of_genes=genenums[:,0][gene_boolean].astype(int).astype(str)

        #get number of genes found
        num_genes=len(genes_found)

        #print number of genes found to screen
        print '%s genes found.' % (num_genes)

        #if genes have been found, save and display according to user specifications
        if (num_genes > 0):
            #save gene names as .npy file
            save_stack=np.column_stack((chrs_of_genes,genes_found))
            np.save(genes_filename, genes_found)

            #for text file and printing to screen, replace chromosome codes 23, 24, and 25 with X, Y, and M respectively, then join chromosomes with genes
            if (genes2txt ==  True) or (genes2screen ==  True) :
                chrs_of_genes[chrs_of_genes=='23']='X'
                chrs_of_genes[chrs_of_genes=='24']='Y'
                chrs_of_genes[chrs_of_genes=='25']='M'

                save_stack=np.column_stack((chrs_of_genes,genes_found))

            #save gene names as text file if requested
            if genes2txt == True:
                fn=genes_filename + '.txt'
                np.savetxt(fn,save_stack, fmt='%.s')

            #print overlapping SNPs to screen if requested
            if (genes2screen == True):
                print "genes found:"
                print "chr  gene"
                for i in xrange(0,num_genes):
                    print chrs_of_genes[i] + '  ' + genes_found[i]

    return None



'''
ncbi2py: function for getting genes in the vicinity of a list of snps, such as those generated by evo_triangulator
    parameters:
        snplist: numpy array or list of snp numbers to find nearby genes for
        bp_range: range of base pairs on either side of snp in which to search for genes
        data_verbose: if True,  list of genes for each snp will be a list of dictionaries of various additional gene data:
            gene name, chromosome, description, aliases, gene start position, gene stop position,summary
            keyed as, respectively:
                'Name','Chr','Description','Alias','Start','Stop','Summary'
        gene2screen: if True, will print summary of genes found to screen
        file_folder: directory in which to save gene dataframes
        genefile: file name to save dataframe
    returns:
        with data_verbose == True: a list of two dataFrames, first one with simple lists of genes as final entry for each snp, second with lists of dictionaries with gene information for each snp (as described above)
        with data_verbose == False: just returns simple data frame


'''



def genefind_ncbi(snplist=[123434,12343557,2342342],bp_range=100000, data_verbose=True, gene2screen=True, file_folder="current", genefile="geneDF.p"):

    import re
    import urllib2
    import xml.etree.cElementTree as ET
    import time
    import numpy as np
    import pandas as pd
    import cPickle as pickle

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)

    #convert list/np.array to list of strings
    snplist=[str(int(i)) for i in snplist]
    bp_range=abs(bp_range)
    #Get SNP location according to NCBI
    #make snp list comma seperated string
    snp_call=",".join(snplist)
    #create url for NCBI snp id data
    snp_handle='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id={}&rettype=docset&tool=evotri&email=reedharder@dartmouth.edu'.format(snp_call)

    #download file
    snp_file=urllib2.urlopen(snp_handle)

    #read and close snp file
    snp_data=snp_file.read()
    snp_file.close()
    #pause as to not irritate NCBI
    time.sleep(.334)

    # parse snp file for location information for each snp, put into snp_matrix
    snp_table=snp_data.split('\n\n')
    # initialize snp data list
    snp_rows=[]
    for snp in snp_table:
        id_match=re.search("SNP_ID=.*?$",snp,flags=re.MULTILINE)
        if id_match != None:
            snpid= int(id_match.group(0)[7:])
            loc_match=re.search("POSITION=.*?$",snp,flags=re.MULTILINE).group(0)[9:].split(":") #find chromosome and position in string
            #create dictionary to use as DataFrame row: with snp id, chr, position, lower bound of gene search, upper bound of gene search, number of gene hits and a gene list as entries,
            snp_data_dict={'snp':snpid,'chr':loc_match[0],'loc':int(loc_match[1]),'lbound':max(int(loc_match[1])-bp_range,0),'ubound':int(loc_match[1])+bp_range, 'geneHits':0, 'genes':[]}
            #append dictionary to list of snp rows
            snp_rows.append(snp_data_dict)

    #create snp data frome
    snpdf=pd.DataFrame(snp_rows)

    # create copy to use as verbose DF
    if data_verbose == True:
        snpdfVerbose=pd.DataFrame.copy(snpdf)

    #initialize list to keep track of genes found and counter of number of genes found
    genecount=0
    genesfound=[]

    #fill dictionary
    for i in range(0,len(snp_rows)):
        #generate esearch region call using chr and range for each snp
        term="{}[chr]+AND+{}:{}[chrpos]+AND+human[orgn]+AND+alive[prop]".format(snpdf['chr'][i],snpdf['lbound'][i],snpdf['ubound'][i])
        html="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}&retmax=30000&tool=evotri&email=reedharder@dartmouth.edu".format(term)
        #get data from esearch call
        gene_file=urllib2.urlopen(html)
        #parse xml, extract gene id's
        root=ET.parse(gene_file).getroot()
        gene_file.close()
        time.sleep(.334)
        gene_idlist=[] #initialize list of genes
        for child in root.find('IdList'):
            gene_idlist.append(child.text) #append each gene Id in IdList to gene_list
        #convert to comma seperated string
        gene_idlist=",".join(gene_idlist)

        #look up gene info from gene id's
        gene_datafile=urllib2.urlopen("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}&retmax=30000&tool=evotri&email=reedharder@dartmouth.edu".format(gene_idlist))
        root = ET.parse(gene_datafile).getroot()
        gene_datafile.close()
        time.sleep(.334)

        current_snp_genes=[] #list to store genes for this snp
        if data_verbose== True:
            gene_data_list=[] #initialize a list of dictionaries of gene data as well
        #loop through gene id's
        for docsum in root:
            #append gene names to list
            try:
                current_snp_genes.append(docsum.find("./Item/[@Name='Name']").text)
                if data_verbose ==True:
                    #get various other data on gene
                    Name=docsum.find("./Item/[@Name='Name']").text
                    Chrom=docsum.find("./Item/[@Name='Chromosome']").text
                    Description=docsum.find("./Item/[@Name='Description']").text
                    OtherAliases=docsum.find("./Item/[@Name='OtherAliases']").text
                    Start=docsum.find(".//Item/[@Name='ChrStart']").text
                    Stop=docsum.find(".//Item/[@Name='ChrStop']").text
                    Summary=docsum.find("./Item/[@Name='Summary']").text
                    #create dictionary of relevant gene data append to list for this snp
                    genedict={'Name':Name,'Chr':Chrom,'Description':Description,'Alias':OtherAliases,'Start':Start,'Stop':Stop,'Summary':Summary}
                    gene_data_list.append(genedict)
            except AttributeError:
                pass

        snpdf['geneHits'][i]=len(current_snp_genes) #add number of gene hits for this snp to data frame
        genecount=genecount + len(current_snp_genes) # add number of genes found for this snp to grand total
        genesfound = genesfound + current_snp_genes # add list of genes found to grand list
        snpdf['genes'][i]=current_snp_genes #add list of genes found for this snp to dataFrame
        #add number of hits and gene data to Verbose dataFrame if specified
        if data_verbose == True:
            snpdfVerbose['geneHits'][i]=len(current_snp_genes)
            snpdfVerbose['genes'][i]=gene_data_list

    #print overall data to screen if specified

    if gene2screen == True:
        print("{} hits, {} unique genes found".format(genecount, len(set(genesfound)))) # number of total and unique hits
        print("\n")
        print(" unique genes:")
        print("\n")
        for gene in set(genesfound): #print each unique gene
            print(gene)

    if genecount > 0:
        #return simple DataFrame or both simple an verbose, depending on specification
        cols=['snp','chr','loc','lbound','ubound','geneHits','genes'] #order of columns
        if data_verbose == False:
            DF = snpdf[cols]
        if data_verbose == True:
            DF = [snpdf[cols],snpdfVerbose[cols]]
        pickle.dump(DF, open(genefile, "w"))
    else:
        DF= None

    return DF #return DataFrame(s)

'''
hit_plotter2D: function to plot number of hits (genes and/or snps) that the evolutionary triangulator finds against a changing Fst on one
of the population-pair axes

parameters:
    backend: enter string to select matplotlib backend for generating figure
    granularity: how many points to plot between an Fst of 0 and 1
    fst_file: Fst data to perform evolutionary tringulation data on. Should be a name of dictionary containing Fst data for pairings of three populations (such as Fst_out from evo_tri_data() )
    pops: list of population abbreviations for population 1, population 2, and population 3, in that order
    pop12lim,pop13lim,pop23lim: Fst threshold (string with operator [>,>,>=,<=] followed by Fst between 0 and 1) for pop1-pop2 Fsts, pop1-pop3 Fsts, and pop2-pop3 Fsts respectively
        One and only one of these thresholds should be an operator followed by character X. This will be the axis along which the value changes.
    loc_data: if plotting gene hits, provide a string with .npy or .npz file name (without extension), with file containing SNP location data, in format of pop_out from evo_tri_data().
    bp_range: if plotting gene hits, provide an integer specifing how many bases from a SNP a gene must be to be considered a hit.
    file_folder: folder from which to save and load files. Default to current directory (input: "current")

'''

def hit_plotter2D(backend="TkAgg", granularity=10, fst_file='fst_data1.npz', pops=['CEU','YRI','GIH'], plot_snps=True, plot_genes=True,
    pop12lim=">=X", pop13lim=">=.45", pop23lim="<=.05",
    loc_data='pop_data1', bp_range=100000, file_folder="current"):

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)

    #import numpy and matplotlib with chosen backend
    import numpy as np
    import matplotlib
    matplotlib.interactive(True)
    matplotlib.use(backend)
    import matplotlib.pyplot as plt

    #load fst data
    fst_dict=np.load(fst_file)

    #vector of FSTs to get number of SNPs and/or genes found at
    xvect = np.linspace(0,1,granularity)

    #find which population pair will be the variable axis. Will take first threshold 'X', default will be pop1-pop2
    lims=[pop12lim,pop13lim,pop23lim]
    axis=0
    for i in xrange(0,3):
        #check if
        if lims[i][-1] == 'X':
            axis=i
            break

    #take operator for variable threshold as the charachters given before 'X'
    variable_operator=lims[axis][:-1]

#generate list of operators and thresholds for each population pair
    op_inputs = [] #list of operators
    fst_bounds = [] #and list of corresponding Fst values
    #parse bounds inputs
    for i in xrange(0,3):
        lim = lims[i]
        # if given population pair will be static, parse input normally
        if i!= axis:
            #make sure operator is valid
            if lim[0] != ">" and lim[0] != "<":
                print "invalid operator in given bounds"
                import sys
                sys.exit()
            #seperate operator from Fst (Fst should start at second or third char of input string
            if lim[1] == '=':
                op=lim[:2]
                fst_bound=float(lim[2:])
            else:
                op=lim[:1]
                fst_bound=float(lim[1:])
            #make list of operators and thresholds
            op_inputs.append(op)
            fst_bounds.append(fst_bound)
        #if the axis is variable, append placeholder value for threshold
        else:
            op_inputs.append(variable_operator)
            fst_bounds.append(-1)

    #dictionary of possible operators
    import operator #
    ops = {">": operator.gt,
         "<": operator.lt,
         ">=": operator.ge,
         "<=": operator.le}

    #get three operators as list of functions
    op_func=[]
    for i in xrange(0,3):
        op_func.append(ops[op_inputs[i]])

    # import stripped down evolutionary triangulation and gene finder functions
    import evo_counters

# generate y1 vector of number of SNP hits for x vector of Fsts. Generate y2 vector of gene hits if requested

    #initialize y vector(s)
    xlen= granularity
    y1=np.zeros(xlen)
    if plot_genes:
        y2=np.zeros(xlen)
        #load population data
        fname = loc_data + '.npz'
        loc_dict = np.load(fname)
        #preload gene data
        import evo_tri
        chdir(path.dirname(evo_tri.__file__))
        genenums=np.load('genefinder_default_genelocs.npy')
        genenames=np.load('genefinder_default_genenames.npy')
        chdir(file_folder)

    #fill y vector(s)
    for i in xrange(0,xlen):
        #insert variable Fst value into variable into variable axis
        fst_bounds[axis]=xvect[i]
        # run evolutiory triangulation with specified fst data, list of pops, list of operators, list of Fst thresholds,
        #and whether or not to return a list of SNPs found for subsequent gene finding
        snps_out=evo_counters.stripped_down_evo_triangulator(fst_dict,pops,op_func,fst_bounds,plot_genes)
        #insert number of snps found into y1 vector
        y1[i]=snps_out[0]

        if plot_genes:
            #if also plotting genes, run gene finder with output snp list, specified location data, bp range, and gene location and name data
                #insert output (number of genes found) into y2 vector
            y2[i]=evo_counters.stripped_down_genefinder(snps_out[1],loc_dict,bp_range,genenames, genenums)

#plot number of SNP and/or gene hits

    #get strings of pop combos
    pop_combos=[pops[0]+'-'+pops[1],pops[0]+'-'+pops[2],pops[1]+'-'+pops[2]]
    #get indices of pop combos that are static
    static_index=[x for x in range(0,3) if x!=axis]
    #generate appropriate plot title and xlabel
    title=pop_combos[static_index[0]]+' '+'Fst'+op_inputs[static_index[0]]+str(fst_bounds[static_index[0]])+' & '+pop_combos[static_index[1]]+' '+'Fst'+op_inputs[static_index[1]]+str(fst_bounds[static_index[1]])
    xlab=pop_combos[axis]+' '+'Fst'+op_inputs[axis]

    #plot appropriate combination of data
    if plot_snps and plot_genes:

        plt.figure(1)
        plt.subplot(211)
        plt.plot(xvect, y1, 'bo')
        plt.xlabel(xlab)
        plt.ylabel('snp hits')
        plt.title(title)

        plt.subplot(212)
        plt.plot(xvect, y2, 'ro')
        plt.xlabel(xlab)
        plt.ylabel('gene hits')
        plt.show()

    elif plot_snps and not plot_genes:

        plt.figure(1)
        plt.plot(xvect, y1, 'bo')
        plt.xlabel(xlab)
        plt.ylabel('snp hits')
        plt.title(title)
        plt.show()

    elif plot_genes and not plot_snps:

        plt.figure(1)
        plt.plot(xvect, y2, 'ro')
        plt.xlabel(xlab)
        plt.ylabel('gene hits')
        plt.title(title)
        plt.show()


    return None







'''
hit_plotter3D: function to plot number of hits (genes and/or snps) that the evolutionary triangulator finds against a changing Fst on two
of the population-pair axes. This function can take a long time to run, depending on granularity.

parameters:
    backend: enter string to select matplotlib backend for generating figure
    granularity: how many points to plot between an Fst of 0 and 1
    fst_file: Fst data to perform evolutionary tringulation data on. Should be a name of dictionary containing Fst data for pairings of three populations (such as Fst_out from evo_tri_data() )
    pops: list of population abbreviations for population 1, population 2, and population 3, in that order
    pop12lim,pop13lim,pop23lim: Fst threshold (string with operator [>,>,>=,<=] followed by Fst between 0 and 1) for pop1-pop2 Fsts, pop1-pop3 Fsts, and pop2-pop3 Fsts respectively
        One of these thresholds should be an operator followed by character X. This will be the first axis along which the value changes.
        Another of these thresholds shoul be an operator followed by character Y. This will be the first axis along which the value changes.
    loc_data: if plotting gene hits, provide a string with .npy or .npz file name (without extension), with file containing SNP location data, in format of pop_out from evo_tri_data().
    bp_range: if plotting gene hits, provide an integer specifing how many bases from a SNP a gene must be to be considered a hit.
    gene_maxZ: max of z-axis range for number of genes
    snp_maxZ: max of z-axis range for number of genes
    file_folder: folder from which to save and load files.

'''




def hit_plotter3D(backend="TkAgg", granularity=10, fst_file='fst_data1.npz', pops=['CEU','YRI','GIH'], plot_snps=True, plot_genes=True,
    pop12lim=">=X", pop13lim=">=Y", pop23lim="<=.05", gene_maxZ=100, snp_maxZ=300,
    loc_data='pop_data1', bp_range=100000, file_folder="current"):

    #use current directory if another is not specified
    if file_folder == "current":
        file_folder=getcwd()

    chdir(file_folder)

    #import numpy and matplotlib with chosen backend
    import numpy as np
    import matplotlib
    matplotlib.interactive(True)
    matplotlib.use(backend)
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    #load fst data
    fst_dict=np.load(fst_file)

    #matrices of FSTs to get number of SNPs and/or genes found at
    x = np.linspace(0,1,granularity)
    y = np.linspace(0,1,granularity)
    xx, yy = np.meshgrid(x,y)

    #find which population pairs will be the variable axes.
    lims=[pop12lim,pop13lim,pop23lim]
    axis1=0
    for i in xrange(0,3):
        #check if
        if lims[i][-1] == 'X':
            axis1=i
            break

    axis2=1
    for i in xrange(0,3):
        #check if
        if lims[i][-1] == 'Y':
            axis2=i
            break

    #take operator for variable threshold as the charachters given before 'X'
    variable_operator1=lims[axis1][:-1]
    variable_operator2=lims[axis2][:-1]

#generate list of operators and thresholds for each population pair
    op_inputs = [] #list of operators
    fst_bounds = [] #and list of corresponding Fst values
    #parse bounds inputs
    for i in xrange(0,3):
        lim = lims[i]
        # if given population pair Fst will be variable, append placeholder value for threshold
        if i==axis1:
            op_inputs.append(variable_operator1)
            fst_bounds.append(-1)
        elif i==axis2:
            op_inputs.append(variable_operator2)
            fst_bounds.append(-1)
        #if population pair Fst will be static, parse normally
        else:
            #record static index
            static_index=i
            #make sure operator is valid
            if lim[0] != ">" and lim[0] != "<":
                print "invalid operator in given bounds"
                import sys
                sys.exit()
            #seperate operator from Fst (Fst should start at second or third char of input string)
            if lim[1] == '=':
                op=lim[:2]
                fst_bound=float(lim[2:])
            else:
                op=lim[:1]
                fst_bound=float(lim[1:])
            #make list of operators and thresholds
            op_inputs.append(op)
            fst_bounds.append(fst_bound)



    #dictionary of possible operators
    import operator #
    ops = {">": operator.gt,
         "<": operator.lt,
         ">=": operator.ge,
         "<=": operator.le}

    #get three operators as list of functions
    op_func=[]
    for i in xrange(0,3):
        op_func.append(ops[op_inputs[i]])

    # import stripped down evolutionary triangulation and gene finder functions
    import evo_counters

# generate Z1 matrix  of number of SNP hits for x vector of Fsts. Generate z2 vector of gene hits if requested

    #initialize y vector(s)
    Z1=np.zeros([granularity,granularity])
    if plot_genes:
        Z2=np.zeros([granularity,granularity])
        #load population data
        fname = loc_data + '.npz'
        loc_dict = np.load(fname)
        #load gene data
        import evo_tri
        chdir(path.dirname(evo_tri.__file__))
        genenums=np.load('genefinder_default_genelocs.npy')
        genenames=np.load('genefinder_default_genenames.npy')
        chdir(file_folder)

    #fill y vector(s)
    for i in xrange(0,granularity):
        for j in xrange(0,granularity):
            #insert variable Fst value into variable into variable axis
            fst_bounds[axis1]=xx[i,j]
            fst_bounds[axis2]=yy[i,j]
            # run evolutiory triangulation with specified fst data, list of pops, list of operators, list of Fst thresholds,
            #and whether or not to return a list of SNPs found for subsequent gene finding
            snps_out=evo_counters.stripped_down_evo_triangulator(fst_dict,pops,op_func,fst_bounds,plot_genes)
            #insert number of snps found into y1 vector
            Z1[i,j]=snps_out[0]

            if plot_genes:
                #if also plotting genes, run gene finder with output snp list, specified location data, and bp range
                    #insert output (number of genes found) into y2 vector
                Z2[i,j]=evo_counters.stripped_down_genefinder(snps_out[1],loc_dict,bp_range,genenames,genenums)

#plot number of SNP and/or gene hits

    #get strings of pop combos
    pop_combos=[pops[0]+'-'+pops[1],pops[0]+'-'+pops[2],pops[1]+'-'+pops[2]]
    #generate appropriate plot title and xlabel
    title=pop_combos[static_index]+' '+'Fst'+op_inputs[static_index]+str(fst_bounds[static_index])
    xlab=pop_combos[axis1]+' '+'Fst'+op_inputs[axis1]
    ylab=pop_combos[axis2]+' '+'Fst'+op_inputs[axis2]


    #plot appropriate combination of data
    if plot_snps and plot_genes:

        fig=plt.figure()
        ax = fig.add_subplot(211, projection='3d')
        ax.plot_wireframe(xx, yy, Z1)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_zlabel('snp hits')
        ax.set_title(title)
        ax.auto_scale_xyz([0, 1], [0,1], [0, snp_maxZ])


        ax = fig.add_subplot(212, projection='3d')
        ax.plot_wireframe(xx, yy, Z2)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_zlabel('gene hits')
        ax.auto_scale_xyz([0, 1], [0,1], [0, gene_maxZ])

        plt.show()

    elif plot_snps and not plot_genes:

        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(xx, yy, Z1)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_zlabel('snp hits')
        ax.set_title(title)
        ax.auto_scale_xyz([0, 1], [0,1], [0, snp_maxZ])

        plt.show()

    elif plot_genes and not plot_snps:

        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(xx, yy, Z2)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_zlabel('gene hits')
        ax.set_title(title)
        ax.auto_scale_xyz([0, 1], [0,1], [0, gene_maxZ])

        plt.show()


    return None

