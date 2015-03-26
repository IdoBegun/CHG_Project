#!/usr/bin/python

import os

from IBD.reference_data import *
from IBD.simulator_data import *
from IBD.translator import *
from IBD.common import *
from IBD.global_params import *
from IBD.phase import *

################################################################################
#                                    MAIN                                      #
################################################################################

if chromsToCompute == 0:
    chromsToCompute = numChrom

print "################################################################################"
print "Phase %s: Preprocessing reference data files" % phase
''
print "--> Loading SNP data..." ,
if os.path.exists(snpDataDirectory):
    print "Skipping"
else:
    simplify_snp_data(snpInfoFile, snpDataDirectory)
    print "Done"

snpCount = count_snps_in_chrom(snpDataDirectory)

print "--> Simplfying reference data...",
if os.path.exists(refDataDirectory):
    print "Skipping"
else:
    print "\n--> Simplfying European reference data..."
    simplify_ref_data(ceuFile, snpCount, refDataDirectory)
    if not DEBUG:
        print "\nDone"
    print "--> Simplfying African reference data..."
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
        create_translator(refDataDirectory, translationDirectory, snpCount, \
                          populationNames, i + 1)
        if DEBUG:
            print "Done."
    if not DEBUG:
        print "\nDone"

print "--> Translating reference data...", 
if os.path.exists(translatedRefDataDirecotry):
    print "Skipping"
else:
    for i in range(chromsToCompute):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpCount, \
                                            i + 1)
        if DEBUG:
            print "--> --> Translating reference data for chromosome "\
                    "%s/%s..." % (i+1, numChrom) ,
        for name in populationNames:
            translate_ref_data(refDataDirectory, translatedRefDataDirecotry, \
                               hapDict, name, i + 1)
        if DEBUG:
            print "Done."
    if not DEBUG:  
        print "\nDone"    


phase += 1   

print "################################################################################"
print "Phase %s: Preprocessing reference data files" % phase

if not os.path.exists(processedDataDirectory):
    os.makedirs(processedDataDirectory)


print "--> Creating LD windows..."
for chrom in range(chromsToCompute):
    chromProcessedDirectory = processedDataDirectory + '/' + str(chrom + 1)
    print "--> --> Creating windows for chromosome %s..." % (chrom + 1),
    if os.path.exists(chromProcessedDirectory):
        print "Skipping"
    else:
        print ""
        os.makedirs(chromProcessedDirectory)
        hapData = []
        for name in populationNames:
            filename = translatedRefDataDirecotry + '/' + name + '_' + \
                        str(chrom + 1)
            refHaps = read_translated_chrom_data(filename)
            hapData.append(refHaps)
            
        ld_windows = load_LD_windows(chromProcessedDirectory, hapData, \
                                     indEpsilon)
    
phase += 1   

print "################################################################################"
print "Phase %s: Preprocessing simulated data files" % phase

print "--> Preprocessing input data...",
if os.path.exists(inputDataDirectory):
    print "Skipping"
else:
    simplify_input_data(inputDataFile, snpCount, inputDataDirectory)
    if not DEBUG:
        print "\nDone."

print "--> Translating input data..." ,
if os.path.exists(translatedInputDataDirecotry):
    print "Skipping"
else:
    for i in range(chromsToCompute):
        hapDict = load_trans_dictionary_hap(translationDirectory, snpCount, \
                                            i + 1)
        genDict = load_trans_dictionary_gen(hapDict)
        translate_input_data(inputDataDirectory, \
                             translatedInputDataDirecotry, i + 1, genDict)
        
    print "Done"

personList = load_person_list(inputDataDirectory, personListFilename)
phase += 1

print "################################################################################"
print "Phase %s: Phasing simulated data" % phase

numGeneration = initNumGen
for chrom in range(chromsToCompute):
    #get hapData
    chromProcessedDirectory = processedDataDirectory + '/' + str(chrom + 1)
    if not os.path.exists(chromProcessedDirectory):
        print "Error: Directory %s does not exist" % (chromProcessedDirectory)
        sys.exit();

    hapData = []
    for name in populationNames:
        filename = translatedRefDataDirecotry + '/' + name + '_' + \
                    str(chrom + 1)
        refHaps = read_translated_chrom_data(filename)
        hapData.append(refHaps)
    # divide to windows
    if DEBUG:
        print "--> --> compute_windows: Started for chromosome %s..." \
                "" % (chrom + 1)

    windowList = compute_windows(chromProcessedDirectory, hapData, beagleEpsilon, \
                                 numGeneration, chrom + 1, \
                                 chromProcessedDirectory, windowListFile)
    
    load_ref_data_windows(hapData, chrom + 1, windowList, \
                          phasedDirectory)
    create_beagle_ref_data(hapData, chrom + 1, windowList, \
                           beaglePhaseDirectory, beagleRefFile)
    filename = translatedInputDataDirecotry + '/' + 'chrom_' + str(chrom + 1)
    genData = read_translated_chrom_data(filename)
    create_beagle_sim_data(genData, chrom + 1, windowList, \
                           beaglePhaseDirectory, beagleGenFile)

    nWindows = len(windowList)
    for iWindow in range(nWindows):
        winPath = beaglePhaseDirectory + '/' + str(chrom + 1) + '/' + beagleGenFile + str(iWindow) + '.vcf.gz'
        refPath = beaglePhaseDirectory + '/' + str(chrom + 1) + '/' + beagleRefFile + str(iWindow) + '.vcf'
        outputPath = phasedDirectory + '/chrom' + str(chrom + 1) + '/pop0/win' + str(iWindow)
        command = 'java -Xmx2000m -jar beagle.r1399.jar gt=' + winPath + ' ref=' + refPath + ' out=' + outputPath
    ###############  TODO ##########################
    # run command and gunzip the result

#################################################

    load_beagle_phased_data(chrom, windowList, inDir, inFile, \
                            phasedDirectory, phasedWindowFile, \
                            chromProcessedDirectory)

    if chrom == 0:
        numGeneration = compute_generation(chrom + 1, chromProcessedDirectory)


phase += 1

print "################################################################################"
print "Phase %s: Computing IBD" % phase

workingDir = os.getcwd() + "/" + beaglePhaseDirectory + "/"
personList = read_person_list_file(inputDataDirectory)
numHaps = len(personList) * 2  

for chrom in range(chromsToCompute):
    print "--> Computing IBD for chromosome %s..." % (chrom + 1)
    # Need a function to read the data saved instead of computing it again!
    windowList = compute_windows(chromProcessedDirectory, hapData, beagleEpsilon, \
                                 numGeneration, chrom + 1, \
                                 chromProcessedDirectory, windowListFile)
    
    numWindows = len(windowList)
    
    compute_ibd(ibdExe, workingDir, chrom, numWindows, numHaps, \
                ibdEpsilon, blockSize, ibdThreshold, maxDiff)

phase += 1

print "################################################################################"
print "Phase %s: Exporting results" % phase

resultsFileHandle = open(fileName, 'w')
header = "name1\tname2\tchr_start\tstart\tchr_end\tend\n"
resultsFileHandle.write(header)

for chrom in range(chromsToCompute):
    # Need a function to read the data saved instead of computing it again!
    windowList = compute_windows(chromProcessedDirectory, hapData, beagleEpsilon, \
                                 numGeneration, chrom + 1, \
                                 chromProcessedDirectory, windowListFile)
    
    results = read_ibd_results(beaglePhaseDirectory, chrom + 1, personList, \
                               windowList, blockSize)
    export_restuls(resultsFileHandle, chrom, personList, results)
    

resultsFileHandle.close()

phase += 1

print "################################################################################"
