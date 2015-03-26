
import sys, re, string, os, math, subprocess

################################################################################
#                               compute_windows                                #
################################################################################

def compute_windows(inDir, hapData, epsilon, numGen, chrom, outDir, outFile):
    '''
    Input:
    inDir - name of directory with LD windows list
    hapData - a list of haplotypes (strings) of reference data
    epsilon - maximum threshold for LD score
    numGen - number of generations of the admixture society
    chrom - number of chromosome
    outDir - name of directory to store the new windows
    outFile - name of file to store the new windows
    
    Output:
    A list of new windows considering: LD windows, number of SNPs in window and window size.
    '''
    expectedRecomb = pow(10, 8)/numGen
    maxWindowLen = expectedRecomb/recombFact
    res = []
    winLD = load_LD_windows(inDir, hapData, epsilon)
    snpCount = count_snps_in_chrom(snpDataDirectory) 
    if not os.path.exists(outDir + '/' + str(chrom)):
        os.makedirs(outDir + '/' + str(chrom))
            
    outputPath = outDir + '/' + str(chrom) + '/' + outFile
    outputFileHandler = open(outputPath, 'w')
    for [winLDStart, winLDEnd] in winLD:
        offsets = get_snp_offsets(snpDataDirectory, chrom, snpCount)
        posStart = winLDStart
        posEnd = posStart
        while (posEnd <= winLDEnd):
            while (posEnd - posStart < maxWindowSNPs ) & \
                  (offsets[posEnd] - offsets[posStart-1] <= maxWindowLen) & \
                  (posEnd != winLDEnd):
                posEnd = posEnd + 1

            res.append([posStart, posEnd])
        # IdoB: Indentation?
        outputFileHandler.write(str(posStart) + ' ' + str(posEnd) + '\n')
            posStart = posEnd  + 1
            posEnd = posStart
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

            outputPath = outDir + '/chrom' + str(chrom) + '/pop' + str(pop + 1) + '/win' + str(iWindow)
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
        # IdoB: That's not how breaking line works when using strings
        outputFileHandler.write('##FORMAT=<ID=GT,Number=1,Type=String,\
                                Description="Genotype">\n')
        outputFileHandler.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')

        nPop = len(hapData[0]) + len(hapData[1])
        for iPerson in range(nPop):
            outputFileHandler.write('\tPERSON' + str(iPerson + 1))
        outputFileHandler.write('\n')

        hapLen = len(hapData[0][0])
        for iSnp in range(hapLen):
            outputFileHandler.write('%s\t%s\t%s\tT\tC\t100\tPASS\t.\tGT' % \
                                    str(chrom), str(iSnp + 1), str(iSnp + 1))
            for pop in hapLen:
                for hapNum in range(len(pop)):
                    if hapNum % 2:
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
                if gen == '0':
                    outputFileHandler.write('\t0/0')
                if gen == '1':
                    outputFileHandler.write('\t0/1')
                if gen == '2':
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

def load_beagle_phased_data(chrom, windowList, inDir, inFile, outDir):
    '''
    Input:
    chrom - chromosome number
    windowList - list of windows dividing the chromosome to work by

    inDir - name of directory contains the admixture haplotype in VCF format
    inFile - name of file contains the admixture haplotype in VCF format
    outDir - name of directory to store the admixture haplotype

    Output:
    None - creates a file for each window containing the admixture haplotype data. Each row represents one haplotype. 
    '''
    # count windows
    winPath = winDir + '/' + str(chrom) + '/' + winFile
    windowFileHandler = open(winPath, 'r')
    numWindows = len(windowList)
    for iWindow in range(numWindows):
        if DEBUG:
                    print "--> --> load_beagle_phased_data: reading file for window %s..." % (iWindow)

    inPath = inDir + '/' + str(chrom) + '/' + inFile + str(iWindow + 1)
    inputFileHandler = open(inPath, 'r')
    # read data from beagle file
    hapData = []
    lineCounter = 0
    for line in inputFileHandler:
        if lineCounter >= 10:
        # IdoB: Indentation?
        splittedLine = line.split()
        snpData = []
        nPerson = len(splittedLine) - 9
        for iPerson in range(nPerson):
            # IdoB: Indentation?
                    personData = splittedLine[9 + iPerson]
            snpData.append(personData[0])
                snpData.append(personData[2])
            hapData.append(snpData)
            
        lineCounter = lineCounter + 1
    if DEBUG:
                    print "--> --> load_beagle_phased_data:creating file for window %s..." % (iWindow)

    outPath = outDir + '/chrom' + str(chrom) + '/pop0/' + str(iWindow)
    outputFileHandler = open(outPath, 'w')
    nHaplotypes = len(hapData[0])
    nSnps = len(hapData)
    for iHap in range(nHaplotypes):
        for iSNP in range(nSnps):
        outputFileHandler.write(hapData[iSNP][iHap])
        outputFileHandler.write('\n')

################################################################################
#                           compute_generation                                 #
################################################################################

def compute_generation(chrom, winDir):
    '''
    Input:
    chrom - chromosome number
    winDir - name of directory contains a list of LD windows

    Output:
    Number of generations passed since admixture.
    '''

    filename = winDir + '/' + str(chrom) + '/ld_windows'
    LDSnps = []
    inputFileHandler = open(filename, 'r')
    for line in inputFileHandler:
    splitLine = line.split()
    LDSnps.append(splitLine[0])
    
    hapData = get_hap_data(LDSnps, chrom, chromProcessedDirectory, windowListFile, phasedDirectory, phasedWindowFile)

    alleleFreq = compute_allele_frequencies(hapData)    
    genVec = []
    for iSnp in range(len(LDSnps)):
    for jSnp in range(iSnp + 1, len(LDSnps)):
        corr = compute_allele_correlation(hapData[2], iSnp, jSnp)
        snpCount = count_snps_in_chrom(snpDataDirectory)
        offsets = get_snp_offsets(snpDataDirectory, chrom, snpCount)
        d = offsets[LDSnps[jSnp] - 1] - offsets[LDSnps[iSnp] - 1]
        tmp = 4 * (corr - alleleFreq[2][iSnp] * alleleFreq[2][jSnp]) / ((alleleFreq[0][iSnp] - alleleFreq[1][iSnp]) * (alleleFreq[0][jSnp] - alleleFreq[1][jSnp]))
        n = math.log(tmp) / (-d)            
        genVec.append(n)
    return sum(x for x in genVec) / len(genVec)
            

################################################################################
#                                 get_hap_data                                         #
################################################################################

def get_hap_data(snpList, chrom, winDir, winFile, subWinDir, subWinFile):
    '''
    Input:
    snpList - list of relevant SNPs
    chrom - chromosome number
    winDir - name of directory contains a list of windows dividing the chromosome
    winFile - name of file contains a list of windows dividing the chromosome
    subWinDir - name of directory contains admixture population haplotype data stored by windows
    subWinFile - name of file contains admixture population haplotype data stored by windows

    Output:
    a list with haplotype data in the relevant SNPs in all populations.
    '''

    hapData = []
    for name in populationNames:
    filename = translatedRefDataDirecotry + '/' + name + '_' + str(chrom)
    fileHandle = open(filename, 'r')
        res = []
        for line in fileHandle:
            splitLine = line.strip()
        hapString = ''
        for iSnp in snpList:
        hapString = hapString + str(splitLine[0][iSnp - 1])
        res.append(hapString)
    
       fileHandle.close()
    hapData.append(res)

    #get data from simulator
    winPath = winDir + '/' + str(chrom) + '/' + winFile
    fileHandler = open(winPath, 'r')
    nLine = 1
    relSubWins = []
    for line in fileHandler:
    splitLine = line.strip()
    if splitLine[0] in snpList:
        relSubWins.append(nLine)
    nLine = nLine + 1
    
    res = []
    for subWin in relSubWins:
    filename = subWinDir + '/' + str(chrom) + '/' + subWinFile + str(subWin)
    fileHandler = open(filename, 'r')
        if len(res) == 0:
        for line in fileHandler:
        splited = line.split()[0]
        res.append(splited[0])
        else:
        nLine = 0
                for line in fileHandler:
            splited = line.split()[0]
            res[nLine] = res[nLine] + splited[0]
            nLine = nLine + 1
    fileHandler.close()
    hapData.append(res)
    return hapData
