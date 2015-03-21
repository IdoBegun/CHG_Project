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

# Location of the IBD executable file
ibdExe = ""

# Location of the BEAGLE executable file
beagleExe = ""

# When calculating dependencies between SNPs, this is the maximum value
# we consider as independent
indEpsilon = 0.00001

phase = 1

# Debug
DEBUG = True
chromsToCompute = 2