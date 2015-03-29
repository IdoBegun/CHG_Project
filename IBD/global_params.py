################################################################################
#                               Global Variable                                #
################################################################################

# Total number of chromosomes
numChrom = 22

# Suffix for input/output files (haplotypes, snpinfo, ibdresults)
inputFileNumber = "0"

# Filename of a person list from the input data
personListFilename = "personList"

# Reference data files
populationFilenameSuffix = ".haps.long"
ceuFile = "ceu.haps.long"
yriFile = "yri.haps.long"
populationNames = []

snpInfoFile = ""
inputDataFile = ""
populationFilenames = ""
inputFileNumber = 0

# Input data files
snpInfoFilePrefix = "snpinfo"
inputDataFilePrefix = "haplotypes"

# Output file
outputDataFilePrefix = "ibdresults1"

# Directories
# Directory of the simplified simulator data
inputDataDirectory = "input_data"
inputDataPrefix = "chrom_"
# Directory of the simplified reference data
refDataDirectory = "ref_data"
# Directory of the simplified SNP info
snpDataDirectory = "snp_data"
snpDataPrefix = "chrom_"
# Directory of the translator
translationDirectory = "translation"
translationaPrefix = "chrom_"
# Directory of the translated reference data
translatedRefDataDirecotry = "ref_data_trans"
# Directory of the translated simulator data
translatedInputDataDirecotry = "input_data_trans"
# Prefix for translated data files
translatedDataPrefix = "chrom_"
# Directory of the processesed data
processedDataDirectory = "processed_data"
# file prefix used for beagle phase input files
beagleRefFile = 'beagle_reference'
beagleSimFile = 'beagle_simulator'
beaglePhasedFile = 'beagle_phased_simulator'
# directory used for beagle phase
beaglePhaseDirectory = 'beagle_phase'

# directory of phased haplotypes by chromosome and window
phasedDirectory = 'phasedSimulatorData'
phasedWindowFile = 'phased_window'

# file prefix of window list
windowListFile = 'ld_windows'

# initial number of generations
initNumGen = 20

# Recombination factor
recombFact = 10

# max number of SNPs allowed in a window
maxWindowSNPs = 300

# Location of the IBD executable file
ibdExe = "exe"

# Location of the BEAGLE executable file
beagleExe = ""

# When calculating dependencies between SNPs, this is the maximum value
# we consider as independent
indEpsilon = 0

# Minimum number of SNPs in each LD independent window (don't set to less than 3!
minInd = 100

# Epsilon used for BEAGLE
beagleEpsilon = 0.00001

# Epsilon for the IBD executable
ibdEpsilon = 0.01

# Block size
blockSize = 20

# Threshold for accepting IBD
ibdThreshold = 1.5

# Maximum difference between IBD scores between populations
maxDiff = 0.3

phase = 1

# A list containing the number of SNPs in each chromosome
snpCount = []

# Debug
DEBUG = True
chromsToCompute = 2