import os, sys

################################################################################
#                              create_translator                               #
################################################################################

def create_translator(inputDir, outputDir, snpCount, populationNames, chrom):
    '''
    Input:
    inputDir - Name of directory containing haplotypes reference data in files
                named <population_name>_<chromosome number>
    outputDir - Name of directory to store the translator data
    snpCount - A list of number of SNPs per chromosome
    populationNames - List containing names of populations
    chrom - Chromosome number
    
    Output:
    None- creates a file named chrom_<chrom> which contains all possible values
            for a SNP in the given chromosome. The first value will be
            translated as 0, and the second as 1
    '''
    
    hap = [[] for y in range(snpCount[chrom - 1])]
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    outputFilename = outputDir + '/chrom_' + str(chrom)
    
    for name in populationNames:
        inputFilename = inputDir + '/' + name + '_' + str(chrom)
        fileHandle = open(inputFilename,'r')
        
        for line in fileHandle:
            splitLine = line.split()
            for i in range(len(splitLine)):
                if splitLine[i] not in hap[i]:
                    hap[i].append(splitLine[i])
        
        fileHandle.close()
    
    outputFileHandle = open(outputFilename, 'w');
    
    for i in range(len(hap)):
        if len(hap[i]) > 2: # Sanity
            print "Error: create_translator: SNP %s of chromosome %s has " \
                    "more than 2 possible values: %s" % (i + 1, chrom, hap[i])
            sys.exit();
        
        outputFileHandle.write(' '.join(hap[i]) + '\n')
    
    outputFileHandle.close()

################################################################################
#                         load_trans_dictionary_hap                            #
################################################################################

def load_trans_dictionary_hap(transDir, snpCount, chrom):
    '''
    Input:
    transDir - Name of the directory containing all translation dictionaries
    snpCount - A list of number of SNPs per chromosome
    chrom - Chromosome number
    
    Output:
    A list of dictionaries. Each entry represents a single SNP in the relevant
    chromosome. The dictionary contains data for each haplotype value.
    '''
    
    inputFilename = transDir + '/chrom_' + str(chrom)
    res = [dict() for y in range(snpCount[chrom - 1])]
    count = 0
    
    fileHandle = open(inputFilename,'r')
    for line in fileHandle:
        hapVals = line.split()
        zippedVals = zip(hapVals, range(len(hapVals)))
        
        for key,val in zippedVals:
            res[count][key] = val
        
        count += 1
    
    return res

################################################################################
#                         load_trans_dictionary_gen                            #
################################################################################

def load_trans_dictionary_gen(hapDict):
    '''
    Input:
    hapDict - A list of dictionaries. Each entry represents a single SNP in the
    relevant chromosome. The dictionary contains data for each haplotype value.
    
    Output:
    A list of dictionaries. Each entry represents a single SNP in the relevant
    chromosome. The dictionary contains data for each genotype value.
    
    Note:
    Some dictionaries may contain a single value.
    '''
    
    length = len(hapDict)
    res = [dict() for y in range(length)]
    
    for i in range(length):
        for x in hapDict[i]:
            for y in hapDict[i]:
                res[i][x + y] = hapDict[i][x] + hapDict[i][y]
    
    return res    

################################################################################
#                             translate_ref_data                               #
################################################################################

def translate_ref_data(inputDir, outputDir, hapDict, populationName, chrom):
    '''
    Input:
    inputDir - Name of directory containing haplotypes reference data in files
                named <population_name>_<chromosome number>
    outputDir - Name of directory to store the translated data
    hapDict - List of dictionaries, containing translations to each SNP
    populationName - Population name
    chrom - Chromosome number
    
    Output:
    None- creates a file identical to the <populationName>_<chrom> file,
            in which the data is translated according to hapDict
    '''
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    inputFilename = inputDir + '/' + populationName + '_' + str(chrom)
    outputFilename = outputDir + '/' + populationName + '_' + str(chrom)
    inputFileHandle = open(inputFilename, 'r')
    outputFileHandle = open(outputFilename, 'w')
    
    length = len(hapDict)
    
    for line in inputFileHandle:
        splitLine = line.split()
        
        for i in range(length):
            translated = str(hapDict[i][splitLine[i]])
            outputFileHandle.write(translated)
        
        outputFileHandle.write('\n')    
            
    inputFileHandle.close()
    outputFileHandle.close()  

################################################################################
#                            translate_input_data                              #
################################################################################

def translate_input_data(inputDir, outputDir, chrom, genDict):
    '''
    Input:
    inputDir - Name of directory containing genotype input data
    outputDir - Name of directory to store the translated data
    chrom - Chromosome number
    genDict - Genotype translation dictionary
    
    Output:
    None- creates a file identical to the original file, where all the genotype
            data is translated according to genDict
    '''
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    filename = "chrom_" + str(chrom)
    inputFilename = inputDir + '/' + filename
    outputFilename = outputDir + '/' + filename
    inputFileHandle = open(inputFilename, 'r')
    outputFileHandle = open(outputFilename, 'w')
    
    length = len(genDict)
    
    for line in inputFileHandle:
        splitLine = line.split()
        for i in range(length):
            translated = str(genDict[i][splitLine[i]])
            outputFileHandle.write(translated)
        
        outputFileHandle.write('\n')
            
    inputFileHandle.close()
    outputFileHandle.close()  
