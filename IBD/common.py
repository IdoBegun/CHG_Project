from global_params import numChrom

################################################################################
#                             count_snps_in_chrom                              #
################################################################################

def count_snps_in_chrom(inputDir):
    '''
    Input:
    inputDir - Name of the directory containing SNP data files for
                each chromosome

    Output:
    A list containing the number of SNPs in each chromosome
    '''

    res = [0 for x in range(numChrom)]
    for chrom in range(numChrom):
        filename = inputDir + "/chrom_" + str(chrom + 1)
        fileHandle = open(filename,'r')
        while fileHandle.readline():
            res[chrom] += 1
        
        fileHandle.close()
        
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
    entry represents a different persson or haplotype
    '''
    
    fileHandle = open(filename, 'r')
    res = []
    for line in fileHandle:
        res.append(line.strip()) 
    
    fileHandle.close()
    return res
