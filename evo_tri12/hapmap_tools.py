# hammap_tools.py: tools for downloading and parsing hapmap files

import numpy as np

'''
allelefreq_readline(): function to parse line of allele_freq hapmap file, takes line file, outputs 1d numpy array or relevant data
'''

def allelefreq_readline(line,chr_float):
    elements=line.rstrip().split(' ') #remove end of line character and split into list by spaces
    try:
        rs_id  = float(elements[0][2:]) #try removing the first to chars of the returned rs string (expected to remove "rs")
    except ValueError: #for abnormal SNP names, remove all non-digit characters to leave rs number only
        import re
        rs_id=float(re.sub("[^0-9]", "", elements[0]))
    #gather data as 1d numpy array with cols as rs number, chromosome, position, allele frequency, and sample size, in that order
    return np.array([rs_id, chr_float, float(elements[2]), float(elements[11]), float(elements[16])])

'''
get_hapmap(): function to download hapmap files for selected populations and chromosomes. Takes list of chromosomes as output by chromeparse(),
hapmap phase
'''

def get_hapmap(chr_groups,phase,pops,savedir):
#connect to Hapmap FTP
    from ftplib import FTP
    import os
    ftp=FTP(host='ftp.hapmap.org', timeout=30)
    ftp.login()
    #try getting list of files from selected phase
    try: #probable name of directory for given phase
        hapmap_dir='hapmap/frequencies/%s/fwd_strand/non-redundant' % phase
        filenames = ftp.nlst(hapmap_dir)
    except: #if fails, try other possible directory configuration
        hapmap_dir = 'hapmap/frequencies/%s' % phase
        filenames = ftp.nlst(hapmap_dir)

    #create list of regular expressions for getting all filenames relevant to user's chromosome and population selection
    import re
    pop_chr_regexes = [('allele_freqs_chr' + x + '_' + y) for x in chr_groups for y in pops]
    allele_freq_files=[]
    #get reduced list of matching filenames
    for pattern in pop_chr_regexes:
        match=False
        #check if pop/chr combo is in list of filenames
        for fn in filenames:
            if re.search(pattern, fn):
                allele_freq_files.append(fn) # if so, append to reduced list
                match=True
        #if no files are found for pop/chr combo, alert user
        if match==False:
            file_meta=pattern.split('_')
            print "File not found for %s %s" % (file_meta[3],file_meta[2])
    #download files
    for filename in allele_freq_files:
        #alert user to download
        currentfile=filename.split('/')
        print "Downloading %s" % currentfile[-1]
        try: #attempt to save file in specified location, with same filename
            saveas=savedir + os.sep + currentfile[-1]
            ftp.retrbinary('RETR ' + filename, open(saveas, 'wb').write)
            print "Success"
        except: #alert user if download fails
            import sys
            print "Failure, error:", sys.exc_info()[0]
    ftp.close()

    return None

#chr_as_float: function for converting from hapmap chromosome coding to numerical coding, getting chromosome number as float from hapmap file name

def chr_as_float(filestring):
    #file name will be of format 'allele_freqs_chr[chromosome]...'
    chr_string=filestring[16]
    try:
        #for numerical chromosomes, simply convert
        chr_num=float(chr_string)
    except ValueError:
        #for X,Y, and mitochondrial chromosome, code as 23, 24 and 25 respectively:
        if chr_string=='X':
            chr_num=23.0
        elif chr_string=='Y':
            chr_num=24.0
        elif chr_string=='M':
            chr_num=25.0
    return chr_num
