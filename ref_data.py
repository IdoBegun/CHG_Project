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
translationDirectory = "translation"
translatedRefDataDirecotry = "ref_data_trans"
translatedInputDataDirecotry = "input_data_trans"

phase = 1

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
        Filename - <person_id>_<chromosome number>
        Content -  Haplotype data for the relevant person for the relevant
        chromosome
    Also creates a file named "personList" which contains a list of all
    person IDs
    '''
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    fileHandle = open(filename, 'r')
    x = fileHandle.read(1)
    st = ""
    outputFileHandle = None
    chrom = 0
    counter = None
    name = None
    personCounter = 0
    personListFileHandle = open(outDir + "/" + personListFilename, 'w')
        
    while x:
        if x in string.whitespace:
            if len(st) == 2:
                if chrom > numChrom: # Sanity
                    print "Error: simplify_input_data: chromosome number mismatch"
                    sys.exit();
                outputFileHandle.write(st)
                counter += 1
                if counter == snpCount[chrom - 1]:
                    outputFileHandle.write('\n')
                    outputFileHandle.close()
                    chrom += 1
                    counter = 0
                    if chrom <= numChrom:
                        outFileName = name + "_" + str(chrom)
                        outputFileHandle = open(outDir + "/" + outFileName,'w')
                else:
                    outputFileHandle.write(' ')
            else:
                #if personCounter != 0:
                #    print "Done."
                name = st
                counter = 0
                chrom = 1
                if personCounter != 0:
                    personListFileHandle.write(' ')
                personListFileHandle.write(name)
                personCounter += 1
                #print "--> simplify_input_data: Started reading person %s..." % (personCounter),
                outFileName = name + "_" + str(chrom)
                outputFileHandle = open(outDir + "/" + outFileName,'w')
            st = ""
        else:
            st += x
        x = fileHandle.read(1)
    
    personListFileHandle.write('\n')
    personListFileHandle.close()
    outputFileHandle.close()
    #print "Done."
    fileHandle.close()

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
            #print "Done."
            #print "--> simplify_ref_data: Started reading data for chromosome %s..." % (chrom),
            chrom += 1
            count = 0
            del data
            data = [[] for y in range(numPersons)]
    
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
            if i != (length - 1):
                outputFileHandle.write(' ')
        
        outputFileHandle.write('\n')    
            
    inputFileHandle.close()
    outputFileHandle.close()  

################################################################################
#                            translate_input_data                              #
################################################################################

def translate_input_data(inputDir, outputDir, filename, genDict):
    '''
    Input:
    inputDir - Name of directory containing genotype input data
    outputDir - Name of directory to store the translated data
    filename - Name of the file to translate
    genDict - Genotype translation dictionary
    
    Output:
    None- creates a file identical to the original file, where all the genotype
            data is translated according to genDict
    '''
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    inputFilename = inputDir + '/' + filename
    outputFilename = outputDir + '/' + filename
    inputFileHandle = open(inputFilename, 'r')
    outputFileHandle = open(outputFilename, 'w')
    
    length = len(genDict)
    
    line = inputFileHandle.readline()
    splitLine = line.split()
        
    for i in range(length):
        translated = str(genDict[i][splitLine[i]])
        outputFileHandle.write(translated)
        if i != (length - 1):
            outputFileHandle.write(' ')
        
    outputFileHandle.write('\n')    
            
    inputFileHandle.close()
    outputFileHandle.close()  

################################################################################
#                              read_snp_data_file                              #
################################################################################

def read_snp_data_file(filename):
    '''
    Input:
    filename - Name of the file containing SNP data

    Output:
    A list containing the number of SNPs in each chromosome
    '''
    
    fileHandle = open(filename,'r')
    fileHandle.readline()
    res = [0 for x in range(numChrom)]
    
    for line in fileHandle:
        splitLine = line.split()
        currChrom = int(splitLine[1][3:])
        res[currChrom - 1] += 1 
    
    fileHandle.close()
    return res

################################################################################
#                              read_person_list_file                              #
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
#                                    MAIN                                      #
################################################################################


print "################################################################################"
print "Phase %s: Processing reference data files" % phase

print "--> Loading SNP data..." ,
snpInfo = read_snp_data_file(snpInfoFile)
print "Done"

print "--> Preprocessing reference data...",
if os.path.exists(refDataDirectory):
    print "Skipping"
else:
    print "\n    --> Processing European reference data..." ,
    simplify_ref_data(ceuFile, snpInfo, refDataDirectory)
    print "Done"
    print "    --> Processing African reference data..." ,
    simplify_ref_data(yriFile, snpInfo, refDataDirectory)
    print "Done"

print "--> Creating translation dictionary..." ,
if os.path.exists(translationDirectory):
    print "Skipping"
else:
    for i in range(numChrom):
        create_translator(refDataDirectory, translationDirectory, snpInfo,\
                           populationNames, i + 1)
    print "Done"

print "--> Translating reference data..." ,
if os.path.exists(translatedRefDataDirecotry):
    print "Skipping"
else:
    for i in range(numChrom):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpInfo,\
                                            i + 1)
        for name in populationNames:
            translate_ref_data(refDataDirectory, translatedRefDataDirecotry,\
                               hapDict, name, i + 1)
        
    print "Done"    

phase += 1

print "################################################################################"
print "Phase %s: Processing input data files" % phase

print "--> Preprocessing input data...",
if os.path.exists(inputDataDirectory):
    print "Skipping"
else:
    simplify_input_data(inputDataFile, snpInfo, inputDataDirectory)
    print "Done"

personList = read_person_list_file(inputDataDirectory)

print "--> Translating input data..." ,
if os.path.exists(translatedInputDataDirecotry):
    print "Skipping"
else:
    for i in range(numChrom):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpInfo,\
                                            i + 1)
        genDict = load_trans_dictionary_gen(hapDict)
        for person in personList:
            filename = person + '_' + str(i + 1)
            translate_input_data(inputDataDirectory,\
                                 translatedInputDataDirecotry, filename,\
                                 genDict)
        
    print "Done"


phase += 1

print "################################################################################"
