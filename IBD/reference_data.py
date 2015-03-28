# This file contains functions relevant only to the reference data

import os, sys
import global_params
from common import compute_allele_correlation, compute_allele_frequencies

DEBUG = global_params.DEBUG

################################################################################
#                              simplify_ref_data                               #
################################################################################

def simplify_ref_data(popName):
    '''
    Input:
    popName - Name of the population
    
    Output:
    None - fills outDir with files in the format of:
        Filename - <population>_<chromosome number>
        Content -  Haplotype data for the relevant chromosome, each line
                    represents a different person
    
    '''
    outDir = global_params.refDataDirectory
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    filename = popName + global_params.populationFilenameSuffix
    fileHandle = open(filename, 'r')
    firstLine = fileHandle.readline().split()
    numHaps = len(firstLine) - 3
    chrom = 1
    count = 0
    data = [[] for y in range(numHaps)]
    
    #TODO: debug?
    #print "--> simplify_ref_data: Started reading data for chromosome %s..." % (chrom),
    
    for line in fileHandle:
        splitLine = line.split()
        if chrom != int(splitLine[1][3:]): # Sanity
            print "Error: simplify_ref_data: chromosome number mismatch."
            print "chrom = %d, int(splitLine[2][3:]) = %d" % (chrom, int(splitLine[1][3:]))
            sys.exit();
        
        count += 1
        for i in range(numHaps):
            data[i].append(splitLine[i + 3])
        
        if count == global_params.snpCount[chrom - 1]:
            outFileName = popName + "_" + str(chrom)
            outputFileHandle = open(outDir + "/" + outFileName, 'w');
            for hap in data:
                outputFileHandle.write(' '.join(hap) + '\n')
            outputFileHandle.close()
            if DEBUG:
                if chrom != 1:
                    print "Done."
                print "--> --> simplify_ref_data: Started reading data for " \
                        "chromosome %s/%s..."  \
                        "" % (chrom, global_params.numChrom),
            chrom += 1
            count = 0
            del data
            data = [[] for y in range(numHaps)]
    
    print "Done."
    fileHandle.close()

################################################################################
#                             compute_LD_ind_windows                           #
################################################################################
def compute_LD_ind_windows(hapData):
    '''
    Input:
    hapData - A list of haplotypa data (strings)
    
    Output:
    A list containing windows in which both start and end indices have LD
    value of at most epsilon
    '''
    
    numPop = len(hapData)
    hapLen = len(hapData[0][0])
    
    alleleFreqs = compute_allele_frequencies(hapData)
    res = []
    
    i = 0
    j = 1
    while (j < hapLen):
        viable = True
        if (((j - 1) - i + 1) < global_params.minInd):
            j += 1
            continue

        for k in range(numPop):
            alleleCorr = compute_allele_correlation(hapData[k], i, j)
        
            D = abs(alleleCorr - (alleleFreqs[k][i] * alleleFreqs[k][j]))
            
            if (D > global_params.indEpsilon):
                viable = False
                break
            
        if viable:
            if global_params.DEBUG:
                print "        --> Adding window %s-%s" % (str(i), str(j - 1))
            res.append([i,j - 1])
            i = j
            j = i + 1
        else:
            j += 1
    
    if ([i, j - 1] != res[-1]):
        res.append([i, j - 1])
    
    # In case the last window is smaller than minInd, merge last two windows
    if ((res[-1][1] - res[-1][0] + 1) < global_params.minInd):
        res[-2][1] = res[-1][1]
        res = res[:-1]
    
    return res
    
################################################################################
#                               load_LD_ind_windows                            #
################################################################################
def load_LD_ind_windows(directory, hapData):
    '''
    Input:
    directory - Directory to save the data to or load from
    hapData - A list of haplotypa data (strings)
    epsilon - Maximum threshold for LD score
    minInd - Minimum number of SNPs in each window
    
    Output:
    List of LD windows. The list is loaded from a file if it exists, otherwise
    the file is created
    '''
    filename = directory + '/ld_windows'
    if not os.path.exists(filename):
        print "    --> Computing LD windows..."
        res = compute_LD_ind_windows(hapData)
        
        fileHandle = open(filename, 'w')
        for window in res:
            fileHandle.write(str(window[0]) + ' ' + str(window[1]) + '\n')
            
        fileHandle.close()
    
    else:
        print "    --> LD independent windows already computed, loading..."
        res = []
        fileHandle = open(filename, 'r')
        
        for line in fileHandle:
            splitLine = line.split()
            start = int(splitLine[0])
            end = int(splitLine[1])
            res.append((start, end))
        
    print "    --> Done"
    return res
    
################################################################################
#                               get_snp_offsets                                #
################################################################################

def get_snp_offsets(chrom):
    '''
    Input:
    inputDir - Name of the directory containing SNP data files
    chrom - Chromosome number

    Output:
    A list containing the offset of each SNPs in the given chromosome
    '''
    
    fileHandle = open(global_params.snpDataDirectory + '/' + \
                      global_params.snpDataPrefix + str(chrom),'r')
    res = [0 for x in range(global_params.snpCount[chrom - 1])]
    count = 0
    
    for line in fileHandle:
        splitLine = line.split()
        offset = int(splitLine[1])
        res[count] = offset
        count += 1
    
    fileHandle.close()
    return res