import os, string, sys
import global_params
from common import *
from math import log

################################################################################
#                             simplify_input_data                              #
################################################################################

def simplify_input_data():
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
    filename = global_params.inputDataFile
    outDir = global_params.inputDataDirectory
    
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
    for chrom in range(global_params.numChrom):
        outputFileHandles.append(open(outDir + '/' + \
                                      global_params.inputDataPrefix \
                                      + str(chrom + 1),'w'))
        
    personListFileHandle = open(outDir + "/" + \
                                global_params.personListFilename, 'w')
    
    # Read each char separately, because reading the entire line
    # crashes the program
    while x:
        if x in string.whitespace:
            if len(st) == 2:
                if chrom > global_params.numChrom: # Sanity
                    print "Error: simplify_input_data: chromosome number " \
                            "mismatch"
                    sys.exit();
                outputFileHandles[chrom - 1].write(st)
                counter += 1
                if counter == global_params.snpCount[chrom - 1]:
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
                    if global_params.DEBUG:
                        print "Done"
                personListFileHandle.write(name)
                personCounter += 1
                if global_params.DEBUG:
                    print "    --> simplify_input_data: Started reading " \
                            "person %s (%s)..." % (name, personCounter),
            st = ""
        else:
            st += x
        x = inputFileHandle.read(1)
    
    personListFileHandle.write('\n')
    if global_params.DEBUG:
        print "Done"
    personListFileHandle.close()
    for fileHandle in outputFileHandles:
        fileHandle.close()
    inputFileHandle.close()
    
################################################################################
#                             simplify_snp_data                                #
################################################################################

def simplify_snp_data():
    '''
    Input:
    None
    
    Output:
    Returns a list of number of SNPs in each chromosome
    
    Note:
    Fills snpDataDirectory with files in the format of:
        Filename - chrom_<chrom number>
        Content -  Each line contains SNP data for the relevant SNP:
        name and offset
    '''
    
    outDir = global_params.snpDataDirectory
     
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    inputFilename = global_params.snpInfoFile
    inputFileHandle = open(inputFilename, 'r')
    inputFileHandle.readline()
    
    outputFileHandles = []
    for chrom in range(global_params.numChrom):
        outputFileHandles.append(open(outDir + '/' + \
                                      global_params.snpDataPrefix + \
                                      str(chrom + 1),
                                      'w'))
    
    res = [0 for x in range(global_params.numChrom)]
    
    for line in inputFileHandle:
        splitLine = line.split()
        chrom = int(splitLine[1][3:])
        
        snpName = splitLine[0]
        snpOffset = splitLine[2]
        res[chrom - 1] += 1
        outputFileHandles[chrom - 1].write(snpName + ' ' + snpOffset + '\n')
    
    for fileHandle in outputFileHandles:
        fileHandle.close()
    
    inputFileHandle.close()
    
    return res


################################################################################
#                            read_person_list_file                             #
################################################################################

def read_person_list_file():
    '''
    Input:
    None

    Output:
    A list containing all person IDs in the personList file
    '''
    
    fileName = global_params.inputDataDirectory + '/' + \
                global_params.personListFilename
    fileHandle = open(fileName, 'r')
    line = fileHandle.readline()
    res = line.strip().split() 
    
    fileHandle.close()
    return res


