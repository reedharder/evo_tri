#chromosome_parser.py: functions for parsing user input chromosome strings,
# including colons to indicate ranges (e.g. 1:10 for chrs 1-10), 'autosome' to
#indicate all autosomal chromosomes.
# Example input string:
#    '1,2,10:14,X' would return list ['1','2','10','11','12','13','14','X'] and relevant regular expressions for file searching


from string import split, count

# function for recursively flattening lists, configured to work in python 3 as well
def flatten(mylist):
    for x in mylist:
        if hasattr(x, '__iter__') and not isinstance(x, str):
            for y in flatten(x):
                yield y
        else:
            yield x

#main parsing function
def chromeparse(chrs):

     #select chrs to include in analysis
    chrs = chrs.translate(None, ' ') #remove any spaces (BEWARE, won't work in Python3) CHANGE PERHAPS
    chr_groups=split(chrs,',') #split input into entries by comma
    #parse input
    for i in range(0,len(chr_groups)):
        S=chr_groups[i]
        Slist=S.split(':') #get start and end of range from form "start:end", if relevant to entry
        #if the start and end ranges are reasonable chromosome numbers, make list entry a list of integers in that range
        if len(Slist)==2 and (int(Slist[0]) >= 1 and int(Slist[0]) <= 22) and (int(Slist[1]) >= 1 and int(Slist[1]) <= 22) and (Slist[1] >= Slist[0]):
            chr_groups[i]=range(int(Slist[0]),int(Slist[1])+1)
        #if user specifies full autosome in entry, make list entry a list of all integers that are valid chromosome labels
        if chr_groups[i] =='autosome':
            chr_groups[i]=range(1,23)
    #destructure list and eliminate duplicates
    chr_groups = list(set(flatten(chr_groups)))
    #convert elements of list to string
    chr_groups = map(str,chr_groups)
    #eliminate invalid chrs from list by making list of valid entries and finding intersection of this and the entered list
    pos_chrs=map(str,[item for sublist in [range(1,23),'X','Y','M'] for item in sublist])
    chr_groups = list(set(chr_groups) & set(pos_chrs))
    #from chr_groups, i.e chrs specified by user, build list of reg. expression patterns to match desired files
    chr_regexes = [('^allele_freqs_chr' + x) for x in chr_groups ]

    #return regular expressions for relevant chromosomes and list of chromosomes
    return [chr_regexes,chr_groups]

