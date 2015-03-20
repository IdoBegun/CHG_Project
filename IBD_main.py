#!/usr/bin/python

import os

import IBD.reference_data as RD
import IBD.simulator_data as SD
import IBD.translator as translator
import IBD.common as common
from IBD.global_params import *

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
    SD.simplify_snp_data(snpInfoFile, snpDataDirectory)
    print "Done"

snpCount = common.count_snps_in_chrom(snpDataDirectory)

print "--> Simplfying reference data...",
if os.path.exists(refDataDirectory):
    print "Skipping"
else:
    print "\n    --> Simplfying European reference data..." ,
    RD.simplify_ref_data(ceuFile, snpCount, refDataDirectory)
    if not DEBUG:
        print "\nDone"
    print "    --> Simplfying African reference data..." ,
    RD.simplify_ref_data(yriFile, snpCount, refDataDirectory)
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
        translator.create_translator(refDataDirectory, translationDirectory, \
                                     snpCount, populationNames, i + 1)
        if DEBUG:
            print "Done."
    if not DEBUG:
        print "\nDone"

print "--> Translating reference data...", 
if os.path.exists(translatedRefDataDirecotry):
    print "Skipping"
else:
    for i in range(chromsToCompute):
        hapDict = translator.load_trans_dictionary_hap(translationDirectory, \
                                                       snpCount, i + 1)
        if DEBUG:
            print "--> --> Translating reference data for chromosome "\
                    "%s/%s..." % (i+1, numChrom) ,
        for name in populationNames:
            translator.translate_ref_data(refDataDirectory, \
                                          translatedRefDataDirecotry, \
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
            refHaps = common.read_translated_chrom_data(filename)
            hapData.append(refHaps)
            
        ld_windows = RD.load_LD_windows(chromProcessedDirectory, hapData, \
                                        epsilon)
    
phase += 1   

print "################################################################################"
print "Phase %s: Preprocessing simulated data files" % phase

print "--> Preprocessing input data...",
if os.path.exists(inputDataDirectory):
    print "Skipping"
else:
    SD.simplify_input_data(inputDataFile, snpCount, inputDataDirectory)
    if not DEBUG:
        print "\nDone."

#personList = read_person_list_file(inputDataDirectory)

print "--> Translating input data..." ,
if os.path.exists(translatedInputDataDirecotry):
    print "Skipping"
else:
    for i in range(chromsToCompute):
        hapDict = translator.load_trans_dictionary_hap(translationDirectory, \
                                                       snpCount, i + 1)
        genDict = translator.load_trans_dictionary_gen(hapDict)
        translator.translate_input_data(inputDataDirectory, \
                                        translatedInputDataDirecotry, \
                                        i + 1, genDict)
        
    print "Done"

personList = SD.load_person_list(inputDataDirectory, personListFilename)
phase += 1

print "################################################################################"
print "Phase %s: Phasing simulated data" % phase

print "--> Not yet implemented"

phase += 1

print "################################################################################"
print "Phase %s: Creating windows" % phase

print "--> Not yet implemented"

phase += 1

print "################################################################################"
print "Phase %s: Computing IBD" % phase

print "--> Not yet implemented"

phase += 1

print "################################################################################"
print "Phase %s: Exporting results" % phase

print "--> Not yet implemented"

phase += 1

print "################################################################################"