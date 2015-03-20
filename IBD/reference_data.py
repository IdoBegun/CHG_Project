# This file contains functions relevant only to the reference data

import os, sys
from global_params import numChrom


################################################################################
#                              simplify_ref_data                               #
################################################################################

def simplify_ref_data(filename, snpCount, outDir):
    '''
    Input:
    filename - Name of the file containing haplotypes reference data
    snpCount - A list of number of SNPs per chromosome
    outDir - Name of the directory which will contain all output files
    
    Output:
    None - fills outDir with files in the format of:
        Filename - <population>_<chromosome number>
        Content -  Haplotype data for the relevant chromosome, each line
                    represents a different person
    
    '''
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    fileHandle = open(filename, 'r')
    populationName = filename.split('.')[0]
    firstLine = fileHandle.readline().split()
    numPersons = len(firstLine) - 3
    chrom = 1
    count = 0
    data = [[] for y in range(numPersons)]
    
    #print "--> simplify_ref_data: Started reading data for chromosome %s..." % (chrom),
    
    for line in fileHandle:
        splitLine = line.split()
        if chrom != int(splitLine[1][3:]): # Sanity
            print "Error: simplify_ref_data: chromosome number mismatch."
            print "chrom = %d, int(splitLine[2][3:]) = %d" % (chrom, int(splitLine[1][3:]))
            sys.exit();
        
        count += 1
        for i in range(numPersons):
            data[i].append(splitLine[i + 3])
        
        if count == snpCount[chrom - 1]:
            outFileName = populationName + "_" + str(chrom)
            outputFileHandle = open(outDir + "/" + outFileName, 'w');
            for person in data:
                outputFileHandle.write(' '.join(person) + '\n')
            outputFileHandle.close()
            if chrom != 1:
                print "Done."
            print "--> --> simplify_ref_data: Started reading data for "\
                    "chromosome %s/%s..." % (chrom, numChrom),
            chrom += 1
            count = 0
            del data
            data = [[] for y in range(numPersons)]
    
    print "Done."
    fileHandle.close()
    #print "Done."

################################################################################
#                          compute_allele_frequencies                          #
################################################################################

def compute_allele_frequencies(hapData):
    '''
    Input:
    hapData - A list of haplotypa data (strings)
    
    Output:
    A list containing allele frequencies. Each entry represents the frequency
    of the allele represented by '1'
    '''
    
    numPop = len(hapData)
    s = len(hapData[0][0])
    
    res = [[] for x in range(numPop)]
    for i in range(s):
        for j in range(numPop):
            #for hap in hapData:
            #    print i, hap[i]
            ones = [int(hap[i]) for hap in hapData[j]]
        
            res[j].append(sum(ones))
    
    for i in range(numPop):
        res[i] = [float(x) / len(hapData[i]) for x in res[i]]
    
    return res

################################################################################
#                          compute_allele_correlation                          #
################################################################################

def compute_allele_correlation(hapData, ind1, ind2):
    '''
    Input:
    hapData - A list of haplotypa data (strings)
    
    Output:
    Probability that both indices have the allele represented by '1'
    '''
    count = 0
    for hap in hapData:
        if (hap[ind1] == '1' and hap[ind2] == '1'):
            count += 1
    
    res = float(count) / len(hapData)
    
    return res

################################################################################
#                               compute_LD_windows                             #
################################################################################

def compute_LD_windows(hapData, epsilon):
    '''
    Input:
    hapData - A list of haplotypa data (strings)
    epsilon - Maximum threshold for LD score
    
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
        for k in range(numPop):
            alleleCorr = compute_allele_correlation(hapData[k], i, j)
        
            D = abs(alleleCorr - (alleleFreqs[k][i] * alleleFreqs[k][j]))
            
            if (D > epsilon):
                viable = False
                break
            
        if viable:
            res.append((i,j))
            i = j + 1
            j = i + 1
        else:
            j = j + 1
    
    if (((i, j - 1)) != res[-1]):
        res.append((i, j - 1))
    
    return res
    
################################################################################
#                                 load_LD_windows                               #
################################################################################

def load_LD_windows(directory, hapData, epsilon):
    '''
    Input:
    directory - Directory to save the data to or load from
    hapData - A list of haplotypa data (strings)
    epsilon - Maximum threshold for LD score
    
    Output:
    List of LD windows. The list is loaded from a file if it exists, otherwise
    the file is created
    '''
    print "--> --> --> Computing LD windows...",
    filename = directory + '/ld_windows'
    if not os.path.exists(filename):
        res = compute_LD_windows(hapData, epsilon)
        
        fileHandle = open(filename, 'w')
        for window in res:
            fileHandle.write(str(window[0]) + ' ' + str(window[1]) + '\n')
            
        fileHandle.close()
        print "Done"
    
    else:
        res = []
        fileHandle = open(filename, 'r')
        
        for line in fileHandle:
            splitLine = line.split()
            start = int(splitLine[0])
            end = int(splitLine[1])
            res.append((start, end))
        
        print "Skipping"
    
    return res
    
################################################################################
#                               get_snp_offsets                                #
################################################################################

def get_snp_offsets(inputDir, chrom, snpCount):
    '''
    Input:
    inputDir - Name of the directory containing SNP data files
    chrom - Chromosome number
    snpCount - A list containing how many SNPs are in each chromosome

    Output:
    A list containing the offset of each SNPs in the given chromosome
    '''
    
    fileHandle = open(inputDir + "/chrom_" + str(chrom),'r')
    res = [0 for x in range(snpCount[chrom - 1])]
    count = 0
    
    for line in fileHandle:
        splitLine = line.split()
        offset = int(splitLine[1])
        res[count] = offset
        count += 1
    
    fileHandle.close()
    return res