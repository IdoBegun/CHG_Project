import os
import global_params
from simulator_data import *

################################################################################
#                             count_snps_in_chrom                              #
################################################################################

def count_snps_in_chrom():
    '''
    Input:
    None

    Output:
    A list containing the number of SNPs in each chromosome
    '''
    
    dirName = global_params.snpDataDirectory
    if not os.path.exists(dirName):
        os.makedirs(dirName)
        print "    --> Processing SNP data..."
        res = simplify_snp_data()
        print "    --> Done"
    else:
        print "    --> SNP data already processed, loading..."
        res = [0 for x in range(global_params.numChrom)]
        for chrom in range(global_params.numChrom):
            filename = dirName + '/' + global_params.snpDataPrefix + \
                         str(chrom + 1)
            fileHandle = open(filename,'r')
            while fileHandle.readline():
                res[chrom] += 1
        
            fileHandle.close()
        print "    --> Done"
        
    return res

################################################################################
#                          read_translated_chrom_data                          #
################################################################################

def read_translated_chrom_data(filename):
    '''
    Input:
    filename - Name of the file containing the data
    
    Output:
    A list containing data from the file (haplotype or genotype data). Each
    entry represents a different person or haplotype
    '''
    
    fileHandle = open(filename, 'r')
    res = []
    for line in fileHandle:
        res.append(line.strip()) 
    
    fileHandle.close()
    return res

################################################################################
#                          compute_allele_frequencies                          #
################################################################################

def compute_allele_frequencies(hapData):
    '''
    Input:
    hapData - A 2d list of haplotype data (strings)
                1st index - population number
                2nd index - haplotype number
    
    Output:
    A list containing allele frequencies. Each entry represents the frequency
    of the allele represented by '1'
    '''
    
    numPop = len(hapData)
    hapLength = len(hapData[0][0])
    
    res = [[] for x in range(numPop)]
    for i in range(hapLength):
        for j in range(numPop):
            snpVals = [int(hap[i]) for hap in hapData[j]]
        
            res[j].append(sum(snpVals))
    
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