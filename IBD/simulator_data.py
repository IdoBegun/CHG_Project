import os, string, sys
from global_params import numChrom, personListFilename, DEBUG
from common import compute_allele_correlation, compute_allele_frequencies
from math import log

################################################################################
#                             simplify_input_data                              #
################################################################################

def simplify_input_data(filename, snpCount, outDir):
    '''
    Input:
    filename - Name of the file containing haplotypes input data
    snpCount - A list of number of SNPs per chromosome
    outDir - Name of the directory which will contain all output files
    
    Output:
    None - fills outDir with files in the format of:
        Filename - chrom_<chrom number>
        Content -  Each line contains genotype data for the relevant person
        for the relevant chromosome
    Also creates a file named "personList" which contains a list of all
    person IDs
    '''
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    inputFileHandle = open(filename, 'r')
    x = inputFileHandle.read(1)
    st = ""
    chrom = 0
    counter = None
    name = None
    personCounter = 0
    outputFileHandles = []
    for chrom in range(numChrom):
        outputFileHandles.append(open(outDir + "/chrom_" + str(chrom + 1),'w'))
        
    personListFileHandle = open(outDir + "/" + personListFilename, 'w')
    
    # Read each char separately, because reading the entire line crashes the program
    while x:
        if x in string.whitespace:
            if len(st) == 2:
                if chrom > numChrom: # Sanity
                    print "Error: simplify_input_data: chromosome number mismatch"
                    sys.exit();
                outputFileHandles[chrom - 1].write(st)
                counter += 1
                if counter == snpCount[chrom - 1]:
                    outputFileHandles[chrom - 1].write('\n')
                    chrom += 1
                    counter = 0
                else:
                    outputFileHandles[chrom - 1].write(' ')
            else:
                #if personCounter != 0:
                #    print "Done."
                name = st
                counter = 0
                chrom = 1
                if personCounter != 0:
                    personListFileHandle.write(' ')
                    if DEBUG:
                        print "Done."
                personListFileHandle.write(name)
                personCounter += 1
                if DEBUG:
                    print "--> --> simplify_input_data: Started reading person %s (%s)..." % (name, personCounter),
            st = ""
        else:
            st += x
        x = inputFileHandle.read(1)
    
    personListFileHandle.write('\n')
    print "Done."
    personListFileHandle.close()
    for fileHandle in outputFileHandles:
        fileHandle.close()
    inputFileHandle.close()
    
################################################################################
#                             simplify_snp_data                                #
################################################################################

def simplify_snp_data(filename, outDir):
    '''
    Input:
    filename - Name of the file containing haplotypes input data
    outDir - Name of the directory which will contain all output files
    
    Output:
    None - fills outDir with files in the format of:
        Filename - chrom_<chrom number>
        Content -  Each line contains SNP data for the relevant SNP:
        name and offset
    '''
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    inputFileHandle = open(filename, 'r')
    inputFileHandle.readline()
    
    outputFileHandles = []
    # TODO: make sure we compute data for all chromosomes
    for chrom in range(numChrom):
        outputFileHandles.append(open(outDir + "/chrom_" + str(chrom + 1),'w'))
    
    
    for line in inputFileHandle:
        splitLine = line.split()
        chrom = int(splitLine[1][3:])
        
        snpName = splitLine[0]
        snpOffset = splitLine[2]
        
        outputFileHandles[chrom - 1].write(snpName + ' ' + snpOffset + '\n')
    
    for fileHandle in outputFileHandles:
        fileHandle.close()
    
    inputFileHandle.close()

################################################################################
#                            read_person_list_file                             #
################################################################################

def read_person_list_file(inputDir):
    '''
    Input:
    inputDir - Name of a directory containing the personList file

    Output:
    A list containing all person IDs in the personList file
    '''
    
    fileName = inputDir + '/' + personListFilename
    fileHandle = open(fileName, 'r')
    line = fileHandle.readline()
    res = line.split() 
    
    fileHandle.close()
    return res


################################################################################
#                                load_person_list                              #
################################################################################

def load_person_list(directory, filename):
    '''
    Input:
    directory - Directory which includes the personList file
    
    Output:
    List of person IDs
    '''
    
    filename = directory + '/' + filename
    fileHandle = open(filename, 'r')
    
    res = []
    for line in fileHandle:
        res.append(line.strip())
    
    return res

