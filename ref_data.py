#!/usr/bin/python

import sys
import re
import string
import os


################################################################################
#                               Global Variable                                #
################################################################################

# Total number of chromosomes
numChrom = 22

# Filename of a person list from the input data
personListFilename = "personList"

# Reference data files
ceuFile = "ceu.haps.long"
yriFile = "yri.haps.long"
populationNames = ["ceu", "yri"]

# Input data files
snpInfoFile = "snpinfo1"
inputDataFile = "haplotypes1"

# Directories
inputDataDirectory = "input_data"
refDataDirectory = "ref_data"
snpDataDirectory = "snp_data"
translationDirectory = "translation"
translatedRefDataDirecotry = "ref_data_trans"
translatedInputDataDirecotry = "input_data_trans"
processedDataDirectory = "processed_data"

# Parameters
epsilon = 0.00001

phase = 1

# Debug
DEBUG = True
chromsToCompute = 1

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
    print "--> --> Computing LD windows...",
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
#                                    MAIN                                      #
################################################################################


print "################################################################################"
print "Phase %s: Processing reference data files" % phase
''
print "--> Loading SNP data..." ,
if os.path.exists(snpDataDirectory):
    print "Skipping"
else:
    simplify_snp_data(snpInfoFile, snpDataDirectory)
    print "Done"

snpCount = count_snps_in_chrom(snpDataDirectory)

print "--> Preprocessing reference data...",
if os.path.exists(refDataDirectory):
    print "Skipping"
else:
    print "\n    --> Processing European reference data..." ,
    simplify_ref_data(ceuFile, snpCount, refDataDirectory)
    if not DEBUG:
        print "\nDone"
    print "    --> Processing African reference data..." ,
    simplify_ref_data(yriFile, snpCount, refDataDirectory)
    if not DEBUG:
        print "\nDone"

print "--> Creating translation dictionary..." ,
if os.path.exists(translationDirectory):
    print "Skipping"
else:
    for i in range(numChrom):
        if DEBUG:
            print "--> --> Translating reference data for chromosome "\
            "%s/%s..." % (i+1, numChrom) ,
        create_translator(refDataDirectory, translationDirectory, snpCount,\
                           populationNames, i + 1)
        if DEBUG:
            print "Done."
    if not DEBUG:
        print "\nDone"

print "--> Translating reference data..."
if os.path.exists(translatedRefDataDirecotry):
    print "Skipping"
else:
    for i in range(numChrom):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpCount,\
                                            i + 1)
        if DEBUG:
            print "--> --> Translating reference data for chromosome "\
                    "%s/%s..." % (i+1, numChrom) ,
        for name in populationNames:
            translate_ref_data(refDataDirectory, translatedRefDataDirecotry,\
                               hapDict, name, i + 1)
        if DEBUG:
            print "Done."
    if not DEBUG:  
        print "\nDone"    

phase += 1

print "################################################################################"
print "Phase %s: Processing input data files" % phase

print "--> Preprocessing input data...",
if os.path.exists(inputDataDirectory):
    print "Skipping"
else:
    simplify_input_data(inputDataFile, snpCount, inputDataDirectory)
    if not DEBUG:
        print "\nDone."

#personList = read_person_list_file(inputDataDirectory)

print "--> Translating input data..." ,
if os.path.exists(translatedInputDataDirecotry):
    print "Skipping"
else:
    for i in range(numChrom):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpCount,\
                                            i + 1)
        genDict = load_trans_dictionary_gen(hapDict)
        translate_input_data(inputDataDirectory, translatedInputDataDirecotry,
                            i + 1, genDict)
        
    print "Done"


phase += 1

print "################################################################################"
print "Phase %s: Processing data for each chromosome" % phase

if chromsToCompute == 0:
    chromsToCompute = numChrom

if not os.path.exists(processedDataDirectory):
    os.makedirs(processedDataDirectory)

for chrom in range(chromsToCompute):
    chromProcessedDirectory = processedDataDirectory + '/' + str(chrom + 1)
    if not os.path.exists(chromProcessedDirectory):
        os.makedirs(chromProcessedDirectory)
    print "--> Creating windows for chromosome %s..." % (chrom + 1)
    hapData = []
    for name in populationNames:
        filename = translatedRefDataDirecotry + '/' + name + '_' + \
                    str(chrom + 1)
        refHaps = read_translated_chrom_data(filename)
        hapData.append(refHaps)
    
    ld_windows = load_LD_windows(chromProcessedDirectory, hapData, epsilon)
    
    


print "################################################################################"
