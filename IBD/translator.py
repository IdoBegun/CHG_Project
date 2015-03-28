import os, sys

import global_params

################################################################################
#                              create_translator                               #
################################################################################

def create_translator(chrom):
    '''
    Input:
    chrom - Chromosome number
    
    Output:
    A list of dictionaries - each entry represents a SNP
    
    Note:
    Creates a file named chrom_<chrom> which contains all possible values
    for a SNP in the given chromosome. The first value will be
    translated as 0, and the second as 1
    '''
    
    hap = [[] for y in range(global_params.snpCount[chrom - 1])]
    res = [dict() for y in range(global_params.snpCount[chrom - 1])]
    
    dirName = global_params.translationDirectory
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    
    outputFilename = dirName + '/' + global_params.translationaPrefix + \
                        str(chrom)
    
    for name in global_params.populationNames:
        inputFilename = global_params.inputDataDirectory + '/' + name + '_' + \
                        str(chrom)
        fileHandle = open(inputFilename,'r')
        
        for line in fileHandle:
            splitLine = line.split()
            for i in range(len(splitLine)): #len(splitLine) = numOfSNPs
                if splitLine[i] not in hap[i]:
                    hap[i].append(splitLine[i])
        
        fileHandle.close()
    
    outputFileHandle = open(outputFilename, 'w');
    
    for i in range(len(hap)):
        if len(hap[i]) > 2: # Sanity
            print "Error: create_translator: SNP %s of chromosome %s has " \
                    "more than 2 possible values: %s" % (i + 1, chrom, hap[i])
            sys.exit();
        
        zippedVals = zip(hap[i], range(len(hap[i])))
        for key,val in zippedVals:
            res[i][key] = val
        
        outputFileHandle.write(' '.join(hap[i]) + '\n')
        
    outputFileHandle.close()
    return res

################################################################################
#                         load_trans_dictionary_hap                            #
################################################################################

def load_trans_dictionary_hap(chrom):
    '''
    Input:
    chrom - Chromosome number
    
    Output:
    A list of dictionaries. Each entry represents a single SNP in the relevant
    chromosome. The dictionary contains data for each SNP value.
    
    Note:
    Output should be in the form of:
    [{'A':0, 'C':1}, {'G':0, 'T':1}]
    '''
    filename = global_params.translationDirectory + '/' + \
                global_params.translationaPrefix + str(chrom)
    
    if not os.path.exists(filename):
        print "    --> Creating translator for chromosome %s..." % str(chrom)
        res = create_translator(chrom)
        print "    --> Done"
    else:
        print "    --> Translator for chromosome %s already exists, " \
                "loading..." % str(chrom)
        res = [dict() for y in range(global_params.snpCount[chrom - 1])]
        count = 0
    
        fileHandle = open(filename,'r')
        for line in fileHandle:
            hapVals = line.split()
            zippedVals = zip(hapVals, range(len(hapVals)))
        
            for key,val in zippedVals:
                res[count][key] = val
        
            count += 1
        print "    --> Done"
    
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
    
    Note:
    Output should be in the form of:
    [{'AA':0, 'AC':1, 'CA':1, 'CC':2}]
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

def translate_ref_data(hapDict, populationName, chrom):
    '''
    Input:
    hapDict - List of dictionaries, containing translations to each SNP
    populationName - Population name
    chrom - Chromosome number
    
    Output:
    None- creates a file identical to the <populationName>_<chrom> file,
            in which the data is translated according to hapDict
    '''
    
    print "    --> Translation reference data for chromosome %s for " \
            "population %s..." % (str(chrom), populationName)
    inputDir = global_params.refDataDirectory
    outputDir = global_params.translatedRefDataDirecotry
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    inputFilename = inputDir + '/' + populationName + '_' + str(chrom)
    outputFilename = outputDir + '/' + populationName + '_' + str(chrom)
    inputFileHandle = open(inputFilename, 'r')
    outputFileHandle = open(outputFilename, 'w')
    
    length = len(hapDict) # Number of SNPs in the chromosome
    
    for line in inputFileHandle:
        splitLine = line.split()
        
        for i in range(length):
            translated = str(hapDict[i][splitLine[i]])
            outputFileHandle.write(translated)
        
        outputFileHandle.write('\n')    
            
    inputFileHandle.close()
    outputFileHandle.close()
    
    print "    --> Done"

################################################################################
#                            translate_input_data                              #
################################################################################

def translate_input_data(chrom, genDict):
    '''
    Input:
    chrom - Chromosome number
    genDict - Genotype translation dictionary
    
    Output:
    None- creates a file identical to the original file, where all the genotype
            data is translated according to genDict
    '''
    inputDir = global_params.inputDataDirectory
    outputDir = global_params.translatedInputDataDirecotry
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    inputFilename = inputDir + '/' + global_params.inputDataPrefix + str(chrom)   
    outputFilename = outputDir + '/' + global_params.translatedDataPrefix + \
                        str(chrom)
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
