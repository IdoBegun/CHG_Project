#!/usr/bin/python

import os, subprocess

import IBD.global_params
from IBD.reference_data import *
from IBD.simulator_data import *
from IBD.translator import *
from IBD.common import *
from IBD.phase import *
from IBD.results import *

################################################################################
#                                    MAIN                                      #
################################################################################

if len(sys.argv) != 6:
    print "Usage: %s <pedigree file> <SNP info file> <genotype data file> " \
            "<population1 data file> <population2 data file>\n"
    sys.exit()

global_params.snpInfoFile = sys.argv[2]
global_params.inputDataFile = sys.argv[3]
global_params.populationFilenames = sys.argv[4:6]
global_params.inputFileNumber = int(global_params.snpInfoFile.split(\
                                            global_params.snpInfoFilePrefix)[1])

global_params.populationNames = [x.split(global_params.populationFilenameSuffix)[0] \
                                 for x in global_params.populationFilenames]



if global_params.chromsToCompute == 0:
    global_params.chromsToCompute = global_params.numChrom

print "################################################################################"
print "Phase %s: Preprocessing reference data files" % phase
''
print "--> Loading SNP data..."
global_params.snpCount = count_snps_in_chrom()
print "--> Done\n"


print "--> Simplfying reference data..."
if os.path.exists(refDataDirectory):
    print "    --> Reference data already simplfied, skipping"
else:
    for popName in global_params.populationNames:
        print "    --> Simplfying reference data for population %s..." \
                % popName
        simplify_ref_data(popName)
        print "    --> Done"
print "--> Done\n"

print "--> Translating reference data..."
if os.path.exists(translatedRefDataDirecotry):
    print "Skipping"
else:
    for i in range(chromsToCompute):
        hapDict = load_trans_dictionary_hap(i + 1)
        for name in global_params.populationNames:
            translate_ref_data(hapDict, name, i + 1)    
print "--> Done"


if not os.path.exists(processedDataDirectory):
    os.makedirs(processedDataDirectory)

print "--> Creating LD independent windows..."
for chrom in range(chromsToCompute):
    chromProcessedDirectory = processedDataDirectory + '/' + str(chrom + 1)
    print "    --> Creating LD independent windows for chromosome %s..." \
            "" % (chrom + 1)
    if os.path.exists(chromProcessedDirectory):
        print "    --> LD independent windows already created, skipping"
    else:
        os.makedirs(chromProcessedDirectory)
        hapData = []
        for name in global_params.populationNames:
            filename = translatedRefDataDirecotry + '/' + name + '_' + \
                        str(chrom + 1)
            refHaps = read_translated_chrom_data(filename)
            hapData.append(refHaps)
            
        ld_windows = load_LD_ind_windows(chromProcessedDirectory, hapData)
print "--> Done"
    
global_params.phase += 1   

print "################################################################################"
print "Phase %s: Preprocessing simulated data files" % phase

print "--> Preprocessing input data..."
if os.path.exists(inputDataDirectory):
    print "    --> Input data already preprocessed, skipping"
else:
    simplify_input_data()

print "--> Done"

print "--> Translating input data..."
if os.path.exists(global_params.translatedInputDataDirecotry):
    print "--> Input data already translated, skipping"
else:
    for i in range(chromsToCompute):
        print "    --> Translating input data for chromosome %s..." \
                "" % str(chrom)
        hapDict = load_trans_dictionary_hap(i + 1)
        genDict = load_trans_dictionary_gen(hapDict)
        translate_input_data(i + 1, genDict)
        print "    --> Done"
print "--> Done"        

global_params.phase += 1

print "################################################################################"
print "Phase %s: Phasing simulated data" % phase

numGeneration = global_params.initNumGen
for chrom in range(chromsToCompute):
    #get hapData
    chromProcessedDirectory = global_params.processedDataDirectory + '/' + \
                                str(chrom + 1)
    if not os.path.exists(chromProcessedDirectory):
        print "Error: Directory %s does not exist" % (chromProcessedDirectory)
        sys.exit();

    hapData = []
    for name in global_params.populationNames:
        filename = translatedRefDataDirecotry + '/' + name + '_' + \
            str(chrom + 1)
        refHaps = read_translated_chrom_data(filename)
        hapData.append(refHaps)
    # divide to windows
    if DEBUG:
        print "--> --> compute_windows: Started for chromosome %s..." \
                "" % (chrom + 1)

    windowList = compute_windows(processedDataDirectory, hapData, \
                                 global_params.beagleEpsilon, numGeneration, \
                                 chrom + 1, global_params.minInd, \
                                 global_params.snpCount, \
                                 processedDataDirectory, \
                                 global_params.windowListFile)
    
    load_ref_data_windows(hapData, chrom + 1, windowList, \
                          global_params.phasedDirectory)
    create_beagle_ref_data(hapData, chrom + 1, windowList, \
                           global_params.beaglePhaseDirectory, \
                           global_params.beagleRefFile)
    filename = translatedInputDataDirecotry + '/' + \
                global_params.translatedDataPrefix + str(chrom + 1)
    genData = read_translated_chrom_data(filename)
    create_beagle_sim_data(genData, chrom + 1, windowList, \
                           global_params.beaglePhaseDirectory, \
                           global_params.beagleSimFile)

    outDir = global_params.phasedDirectory + '/chrom' + str(chrom + 1) + '/pop0'
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    nWindows = len(windowList)
    for iWindow in range(nWindows):
        winPath = global_params.beaglePhaseDirectory + '/' + str(chrom + 1) + \
                    '/' + global_params.beagleSimFile + str(iWindow) + \
                    '.vcf.gz'
        refPath = global_params.beaglePhaseDirectory + '/' + \
                    str(chrom + 1) + '/' + global_params.beagleRefFile + \
                    str(iWindow) + '.vcf'
        outputPath = global_params.phasedDirectory + '/chrom' + \
                        str(chrom + 1) + '/pop0/win' + str(iWindow)
        command = 'java -Xmx2000m -jar beagle.r1399.jar gt=' + winPath + \
                    ' ref=' + refPath + ' out=' + outputPath
        if not os.path.exists(outputPath + '.vcf'):
            output = subprocess.check_output(command, shell=True)
            if output:
                print output
            
            command = "gunzip " + outputPath + ".vcf.gz"
            
            output = subprocess.check_output(command, shell=True)
            if output:
                print output
    
    outputDir = global_params.phasedDirectory + '/chrom' + str(chrom + 1) + \
                '/pop0/'
    load_beagle_phased_data(chrom, windowList, outputDir, outputDir)

    if chrom == 0:
        # update the number of generations according to the first chromosome
        numGeneration = compute_generation(chrom + 1, \
                                           global_params.populationNames, \
                                           global_params.snpCount, \
                                           global_params.processedDataDirectory, \
                                           global_params.translatedRefDataDirecotry)


global_params.phase += 1

print "################################################################################"
print "Phase %s: Computing IBD" % phase

workingDir = os.getcwd() + "/" + beaglePhaseDirectory + "/"
personList = read_person_list_file()
numHaps = len(personList) * 2  

for chrom in range(chromsToCompute):
    print "--> Computing IBD for chromosome %s..." % (chrom + 1)
    chromProcessedDirectory = global_params.processedDataDirectory + '/' + \
                                str(chrom + 1)
    windowList = read_windows(chromProcessedDirectory, windowListFile)
    
    numWindows = len(windowList)
    
    compute_ibd(global_params.ibdExe, workingDir, chrom, numWindows, numHaps)

global_params.phase += 1

print "################################################################################"
print "Phase %s: Exporting results" % phase

resultsFileHandle = open(global_params.outputDataFilePrefix + \
                         global_params.inputFileNumber, 'w')
header = "name1\tname2\tchr_start\tstart\tchr_end\tend\n"
resultsFileHandle.write(header)

for chrom in range(chromsToCompute):
    chromProcessedDirectory = global_params.processedDataDirectory + '/' + \
                                str(chrom + 1)
    windowList = read_windows(chromProcessedDirectory, \
                              global_params.windowListFile)
    
    fileResults = read_ibd_results(beaglePhaseDirectory, chrom + 1)
    ibdResults = process_results(fileResults, chrom, personList, windowList)
    export_results(resultsFileHandle, chrom, personList, ibdResults)
    

resultsFileHandle.close()

global_params.phase += 1

print "################################################################################"
