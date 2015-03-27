
import sys, re, string, os, math, subprocess
from IBD.global_params import *
from IBD.reference_data import *
from IBD.common import *

################################################################################
#                               compute_windows                                #
################################################################################

def compute_windows(inDir, hapData, epsilon, numGen, chrom, minInd, snpCount, outDir, outFile):
    '''
    Input:
    inDir - name of directory with LD windows list
    hapData - a list of haplotypes (strings) of reference data
    epsilon - maximum threshold for LD score
    numGen - number of generations of the admixture society
    chrom - number of chromosome
    minInd - Minimum number of SNPs in each window
    snpCount - List containing number of SNPs in each chromosome
    outDir - name of directory to store the new windows
    outFile - name of file to store the new windows
    
    Output:
    A list of new windows considering: LD windows, number of SNPs in window and window size.
    '''
    expectedRecomb = pow(10, 8)/numGen
    maxWindowLen = expectedRecomb/recombFact
    res = []
    winLD = load_LD_windows(inDir, hapData, epsilon, minInd)
    offsets = get_snp_offsets(snpDataDirectory, chrom, snpCount) 

    for [winLDStart, winLDEnd] in winLD:
        posStart = winLDStart
        posEnd = posStart + 1
        while (posEnd <= winLDEnd):
            while ((posEnd - posStart) < maxWindowSNPs ) and \
                  (offsets[posEnd] - offsets[posStart] <= maxWindowLen) and \
                  (posEnd <= winLDEnd):
                posEnd += 1
            res.append([posStart, posEnd - 1])
            posStart = posEnd
            posEnd = posStart + 1
        res[-1][1] = winLDEnd
        if (res[-1][1] - res[-1][0] < minInd):
            res[-2][1] = winLDEnd
            res = res[:-1]
        
    if not os.path.exists(outDir + '/' + str(chrom)):
        os.makedirs(outDir + '/' + str(chrom))
    outputPath = outDir + '/' + str(chrom) + '/' + outFile
    outputFileHandler = open(outputPath, 'w')
    for [winStart, winEnd] in res:
        outputFileHandler.write(str(winStart) + ' ' + str(winEnd) + '\n')
    outputFileHandler.close()    
    return res


################################################################################
#                               read_windows                                #
################################################################################

def read_windows(inDir, inFile):
    '''
	Input:
	inDir - name of directory with windows list
	inFile - name of file with windows list

	Output:
	winList - a list of windows with limited amount of SNPs 	and length. 
	'''
    path = inDir + '/' + inFile
    fileHandler = open(path, 'r')
    winList = []
    for line in fileHandler:
        split = line.split()
        winStart = split[0]
        winEnd = split[1]
        winList.append([winStart, winEnd])

	return winList

################################################################################
#                            load_ref_data_windows                             #
################################################################################

def load_ref_data_windows(hapData, chrom, windowList, outDir):
    '''
    Input:
    hapData - list of haplotypes (strings) of reference data
    chrom - chromosome number
    windowList - list of windows to work by
    outDir - name of directory to store the reference data

    Output:
    None - creates a file for each window with reference data where each line represents one haplotype.
    '''
    
    nPop = len(hapData)
    for pop in range(nPop):
        for iWindow in range(len(windowList)):
            if DEBUG:
                    print "--> --> load_ref_data_windows: creating file for window %s (pop %s)..." % (iWindow, pop + 1)

            outputPath = outDir + '/chrom' + str(chrom) + '/pop' + str(pop + 1) + '/win' + str(iWindow) + '.txt'
            outFileHandler = open(outputPath, 'w')
            [winStart, winEnd] = windowList[iWindow]
            nHaps = len(hapData[pop])
            for iHap in range(nHaps):
                for iSnp in range(winStart, winEnd + 1):
                    outFileHandler.write(hapData[pop][iHap][iSnp])
                outFileHandler.write('/n')
            outFileHandler.close()
################################################################################
#                            create_beagle_ref_data                            #
################################################################################

def create_beagle_ref_data(hapData, chrom, windowList, outDir, outFile):
    '''
    Input:
    hapData - list of haplotypes (strings) of reference data
    chrom - chromosome number
    windowList - list of windows to work by
    outDir - name of directory to store the reference data
    outFile - name of file to store the reference data
    
    Output:
    None - creates a file for each window with reference data.  Data is stored in VCF format for BEAGLE phasing algorithm.
    ''' 
    outputDir = outDir + '/' + str(chrom)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
    for iWindow in range(len(windowList)):
        if DEBUG:
                    print "--> --> create_beagle_ref_data: creating file for window %s..." % (iWindow)

        filename = outputDir + '/' + outFile + str(iWindow) + '.vcf'
        outputFileHandler = open(filename, 'w')
        [winStart, winEnd] = windowList[iWindow]

        
        outputFileHandler.write('##fileformat=VCFv4.1\n')
        outputFileHandler.write('##FORMAT=<ID=GT,Number=1,Type=String,' \
                                'Description="Genotype">\n')
        outputFileHandler.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')

        nHaps = len(hapData[0]) + len(hapData[1])
        for iPerson in range((nHaps / 2)):
            outputFileHandler.write('\tPERSON' + str(iPerson + 1))
        outputFileHandler.write('\n')

        hapLen = len(hapData[0][0])
        for iSnp in range(hapLen):
        # snpName = snpLocation since we're phasing small windows
            outputFileHandler.write('%s\t%s\t%s\tT\tC\t100\tPASS\t.\tGT' % \
                                    str(chrom), str(iSnp + 1), str(iSnp + 1))
            for pop in hapData:
                for hapNum in range(len(pop)):
                    if not hapNum % 2:
                        outputFileHandler.write('\t' + pop[hapNum][iSnp])
                    else:
                        outputFileHandler.write('|' + pop[hapNum][iSnp])
            
            outputFileHandler.write('\n')       
        
        
        
        outputFileHandler.close()


################################################################################
#                          create_beagle_sim_data                              #
################################################################################

def create_beagle_sim_data(genData, chrom, windowList, outDir, outFile):
    '''
    Input:
    genData: genotype data (strings) of the admixtured population
    chrom - chromosome number
    windowList - list of windows dividing the chromosome to work by
    outDir - name of directory to store the admixture genotype
    outFile - name of file to store the admixture genotype

    Output:
    None - creates a file for each window containing the admixture genotype data. Data is stored in VCF format for BEAGLE phasing algorithm.
    '''
    outputDir = outDir + '/' + str(chrom)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
    for iWindow in range(len(windowList)):
        if DEBUG:
                    print "--> --> create_beagle_sim_data: creating file for window %s..." % (iWindow)

        filename = outputDir + '/' + outFile + str(iWindow) + '.vcf'
        outputFileHandler = open(filename, 'w')
        [winStart, winEnd] = windowList[iWindow]

        
        outputFileHandler.write('##fileformat=VCFv4.1\n')
        outputFileHandler.write('##FORMAT=<ID=GT,Number=1,Type=String,\
                                Description="Genotype">\n')
        outputFileHandler.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
        for iPerson in range(len(genData)):
            outputFileHandler.write('\tPERSON' + str(iPerson + 1))
        outputFileHandler.write('\n')

        genLen = len(genData[0])
        for iSnp in range(genLen):
            outputFileHandler.write('%s\t%s\t%s\tT\tC\t100\tPASS\t.\tGT' % \
                                    str(chrom), str(iSnp + 1), str(iSnp + 1))
            for gen in genData:
                if gen[iSnp] == '0':
                    outputFileHandler.write('\t0/0')
                if gen[iSnp] == '1':
                    outputFileHandler.write('\t0/1')
                if gen[iSnp] == '2':
                    outputFileHandler.write('\t1/1')
            
            outputFileHandler.write('\n')       
        
        
        
        outputFileHandler.close()
        
        command = 'gzip ' + filename
        
        output = subprocess.check_output(command, shell=True)
        if output:
            print output
        
        
################################################################################
#                           load_beagle_phased_data                            #
################################################################################

def load_beagle_phased_data(chrom, windowList, inDir, outDir):
    '''
    Input:
    chrom - chromosome number
    windowList - list of windows dividing the chromosome to work by

    inDir - name of directory contains the admixture haplotype in VCF format
    outDir - name of directory to store the admixture haplotype

    Output:
    None - creates a file for each window containing the admixture haplotype data. Each row represents one haplotype. 
    '''
    # count windows
    numWindows = len(windowList)
    for iWindow in range(numWindows):
        if DEBUG:
            print "--> --> load_beagle_phased_data: reading file for window %s..." % (iWindow)

        inPath = inDir + 'win' + str(iWindow) + '.vcf'
        inputFileHandler = open(inPath, 'r')
        # read data from beagle file
        hapData = [] # format: row = SNP, column = haplotype
        lineCounter = 0
        for line in inputFileHandler:
            if lineCounter >= 10:
                splittedLine = line.split()
                snpData = []
                nPerson = len(splittedLine) - 9
                for iPerson in range(nPerson):
                    personData = splittedLine[9 + iPerson]
                    snpData.append(personData[0])
                    snpData.append(personData[2])
                hapData.append(snpData)
                
            lineCounter = lineCounter + 1
   
        outPath = outDir + 'win' + str(iWindow) + '.txt'
        outputFileHandler = open(outPath, 'w')
        nHaplotypes = len(hapData[0])
        nSnps = len(hapData)
        for iHap in range(nHaplotypes):
            for iSNP in range(nSnps):
                outputFileHandler.write(hapData[iSNP][iHap])
            outputFileHandler.write('\n')
        outputFileHandler.close()

################################################################################
#                           compute_generation                                 #
################################################################################

def compute_generation(chrom, populationNames, snpCount, winDir, translatedRefDataDirecotry):
    '''
    Input:
    chrom - chromosome number
    populationNames - List of population names
    snpCount - List of number of SNPs in each chromosome
    winDir - name of directory contains a list of LD windows
    translatedRefDataDirecotry - Name of the directory containing the translated reference data

    Output:
    Number of generations passed since admixture.
    '''

    filename = winDir + '/' + str(chrom) + '/ld_windows'
    LDSnps = []
    inputFileHandler = open(filename, 'r')
    for line in inputFileHandler:
        splitLine = line.split()
        LDSnps.append(splitLine[0])
    
    hapData = get_hap_data(LDSnps, chrom, populationNames, processedDataDirectory, windowListFile, phasedDirectory, phasedWindowFile, translatedRefDataDirecotry)
    offsets = get_snp_offsets(snpDataDirectory, chrom, snpCount)

    alleleFreq = compute_allele_frequencies(hapData)    
    genVec = []
    for iSnp in range(len(LDSnps)):
        for jSnp in range(iSnp + 1, len(LDSnps)):
            corr = compute_allele_correlation(hapData[2], iSnp, jSnp)
            d = offsets[LDSnps[jSnp] - 1] - offsets[LDSnps[iSnp] - 1]
            tmp = 4 * (corr - alleleFreq[2][iSnp] * alleleFreq[2][jSnp]) / ((alleleFreq[0][iSnp] - alleleFreq[1][iSnp]) * (alleleFreq[0][jSnp] - alleleFreq[1][jSnp]))
            n = math.log(tmp) / (-d)            
            genVec.append(n)
            
    return sum(x for x in genVec) / len(genVec)
            

################################################################################
#                                 get_hap_data                                         #
################################################################################

def get_hap_data(snpList, chrom, populationNames, winDir, winFile, subWinDir, subWinFile, translatedRefDataDirecotry):
    '''
    Input:
    snpList - list of relevant SNPs
    chrom - chromosome number
    populationNames - List of population names
    winDir - name of directory contains a list of windows dividing the chromosome
    winFile - name of file contains a list of windows dividing the chromosome
    subWinDir - name of directory contains admixture population haplotype data stored by windows
    subWinFile - name of file contains admixture population haplotype data stored by windows
    translatedRefDataDirecotry - Name of the directory containing the translated reference data

    Output:
    a list with haplotype data in the relevant SNPs in all populations.
    '''

    hapData = []
    # get the relevant (independent) SNPsfrom the reference data in each haplotype 
    for name in populationNames:
        filename = translatedRefDataDirecotry + '/' + name + '_' + str(chrom)
        fileHandle = open(filename, 'r')
        res = []
        for line in fileHandle:
            strippedLine = line.strip()
            hapString = ''
            for iSnp in snpList:
                hapString = hapString + strippedLine[iSnp]
            res.append(hapString)
    
        fileHandle.close()
        hapData.append(res)

    # get the relevant (independent) SNPsfrom the simulator data in each haplotype
    
    # get list of relevant sub windows (sub windows that contain SNPs from snpList)
    winPath = winDir + '/' + str(chrom) + '/' + winFile
    fileHandler = open(winPath, 'r')
    nSubWin = 0
    relSubWins = []
    for line in fileHandler:
        splitLine = line.split()
        if splitLine[0] in snpList:
            relSubWins.append(nSubWin)
        nSubWin += 1
    fileHandler.close()
    
    res = []
    for subWin in relSubWins:
        filename = subWinDir + '/' + str(chrom) + '/' + subWinFile + str(subWin)
        fileHandler = open(filename, 'r')
        if not res: # subWin is the first window
            for line in fileHandler:
                res.append(line[0])
        else:
            nLine = 0
            for line in fileHandler:
                res[nLine] = res[nLine] + line[0]
                nLine = nLine + 1
    fileHandler.close()
    hapData.append(res)
    return hapData
