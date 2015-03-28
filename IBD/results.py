import sys, subprocess
import global_params
 
from reference_data import get_snp_offsets

################################################################################
#                              read_ibd_results                                #
################################################################################

def read_ibd_results(inputDir, chrom):
    '''
    Input:
    inputDir - Name of a directory containing the personList file
    chrom - Chromosome number
    numBlocks - Number of blocks for the relevant chromosome

    Output:
    A list containing all person IDs in the personList file
    '''
    
    fileName = inputDir + "/chrom" + str(chrom)
    fileHandle = open(fileName, 'r')
    res = []
    
    for line in fileHandle:
        res.append([])
        splitLine = line.split()
        for pair in splitLine:
            splitPair = pair.split('-')
            first = int(splitPair[0]) // 2
            second = int(splitPair[1]) // 2
            res[-1].append((first, second))
    
    return res

################################################################################
#                                 compute_ibd                                  #
################################################################################

def compute_ibd(workingDir, chrom, numWindows, numHaps):
    '''
    Input:
    workingDir - Path of the working directory
    chrom - Chromosome number
    numWindows - Number of windows in the relevant chromosome
    numHaps - Number of haplotypes

    Output:
    None.
    '''

    command = global_params.ibdExe + ' ' + \
                str(global_params.DEBUG) + ' ' + \
                workingDir + ' ' + \
                str(chrom) + ' ' + \
                str(numWindows) + ' ' + \
                str(numHaps) + ' ' + \
                str(global_params.ibdEpsilon) + ' ' + \
                str(global_params.blockSize) + ' ' + \
                str(global_params.ibdThreshold) + ' ' + \
                str(global_params.maxDiff)
    
    if global_params.DEBUG:
        print "    --> Executing:" + command
    
    output = subprocess.check_output(command, shell=True)
    
    print "    --> Done."
    
    if output:
        print output

################################################################################
#                              process_results                                 #
################################################################################

def process_results(resList, chrom, personList, windowList):
    '''
    Input:
    resList - 2d list. 1st index is the block number, 2nd index is
        the number of pairs who have IBD in said block
    chrom - Chromosome number
    personList - Person list

    Output:
    A 3d list:
        1st index - first person index
        2nd index - second person index
        3rd index - number of SNPs ranges in which both persons have IBD
    
    The value of each entry is a pair of indices - indicating IBD between these
    indices
    '''
    numPersons = len(personList)
    offsets = get_snp_offsets(chrom)
    
    ibdBlocks = [[[] for x in range(numPersons)] for y in range(numPersons)]
    
    numBlocks = len(resList)
    
    for blockInd in range(numBlocks):
        for pair in resList[blockInd]:
            first = pair[0]
            second = pair[2]
            ibdBlocks[first][second].append(blockInd)
    
    res = [[[] for x in range(numPersons)] for y in range(numPersons)]
    for person1 in range(numPersons):
        for person2 in range(numPersons):
            if ibdBlocks[person1][person2]:
                start = -2
                end = start
                for block in ibdBlocks[person1][person2]:
                    if (block == end + 1):
                        end = block
                    else:
                        res[person1][person2].append((start,end))
                        start = block
                        end = start
                
                startOffset = offsets[windowList[start][0]]
                endOffset = offsets[windowList[end + global_params.blockSize][1]]
                res[person1][person2].append((startOffset, endOffset))
    
    return res

################################################################################
#                              export_results                                  #
################################################################################

def export_results(fileHandle, chrom, personList, ibdList):
    '''
    Input:
    fileHandle - File handle to the file to write the results in
    chrom - Chromosome number
    personList - Person list
    ibdList - List of all IBD windows

    Output:
    None - write the results to the file
    '''
    
    numPersons = len(personList)
    
    for i in range(numPersons):
        for j in range(numPersons):
            if ibdList[i][j]:
                personName1 = personList[i]
                personName2 = personList[j]
                for pair in ibdList[i][j]:
                    resString = personName1 + '\t' + personName2 + '\tchr' + \
                                str(chrom) + '\t' + str(pair[0]) + '\tchr' + \
                                str(chrom) + '\t' + str(pair[1]) + '\n'
                    
                    fileHandle.write(resString)
    