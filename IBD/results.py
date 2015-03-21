import sys, subprocess

################################################################################
#                              read_ibd_results                                #
################################################################################

def read_ibd_results(inputDir, chrom, numBlocks):
    '''
    Input:
    inputDir - Name of a directory containing the personList file
    chrom - Chromosome number
    numBlocks - Number of blocks for the relevant chromosome

    Output:
    A list containing all person IDs in the personList file
    '''
    
    fileName = inputDir + "/results_" + str(chrom)
    fileHandle = open(fileName, 'r')
    res = [[] for x in range(numBlocks)]
    count = 0
    
    for line in fileHandle:
        splitLine = line.split()
        for pair in splitLine:
            splitPair = pair.split('-')
            first = int(splitPair[0])
            second = int(splitPair[1])
            res[count].append((first, second))
        
        count += 1
    
    # Sanity
    if count != numBlocks:
        print "Error: read_ibd_results: Inconsistent number of blocks in "\
                "IBD results for file chromosome %s" % chrom
        sys.exit();

################################################################################
#                                 compute_ibd                                  #
################################################################################

def compute_ibd(executableName, chrom, numWindows, numHaps, epsilon, \
                blockSize, threshold, maxDiff):
    '''
    Input:
    executableName - Name of the executable to run
    chrom - Chromosome number
    numWindows - Number of windows in the relevant chromosome
    numHaps - Number of haplotypes
    epsilon - Error estimation
    blockSize - Number of blocks in a window
    threshold - Threshold for the block score
    maxDiff - Maximum difference between IBD in two populations

    Output:
    None.
    '''

    command = executableName + ' ' + str(chrom) + ' ' + str(numWindows) + \
                str(numHaps) + ' ' + str(epsilon) + ' ' + str(blockSize) + \
                str(threshold) + ' ' + str(maxDiff)
    
    output = subprocess.check_output(command, shell=True)
    
    if output:
        print output
