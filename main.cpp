#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <assert.h>
#include "LDTree.h"

using namespace std;

void resizeVectors(unsigned int numOfWindows, unsigned int numOfHaplotypes,
                   vector< vector< vector<double> > >& windowScore/*out*/, 
                   vector< vector< vector<bool> > >& maxArgIBDWindow/*out*/, 
                   vector< vector< vector<double> > >& windowDiff/*out*/)
{
  windowScore.resize(numOfWindows);
  maxArgIBDWindow.resize(numOfWindows);
  windowDiff.resize(numOfWindows);
  for (unsigned int i = 0; i < numOfWindows; i++)
  {
    windowScore[i].resize(numOfHaplotypes);
    maxArgIBDWindow[i].resize(numOfHaplotypes);
    windowDiff[i].resize(numOfHaplotypes);
    for (unsigned int j = 0; j < numOfHaplotypes; ++j)
    {
      windowScore[i][j].resize(numOfHaplotypes);
      maxArgIBDWindow[i][j].resize(numOfHaplotypes);
      windowDiff[i][j].resize(numOfHaplotypes);
    }
  }
}

void getDataFromFile(string& fileName, vector<string>& population)
{
  string line;
  ifstream file(fileName.c_str());
  assert(file.is_open());
  while (getline(file,line) )
  {
    population.push_back(line);
  }
  file.close();
}

void getData(unsigned int chromozomeNumber, unsigned int window,
                       vector<string>& simulator/*out*/,
                       vector<string> populationA/*out*/,
                       vector<string> populationB/*out*/)
{
  // read the simulator data
  stringstream ssSimulator;
  ssSimulator << "chrom" << chromozomeNumber << "_win" << window << "_pop0.txt"; // population0 = simulator data
  string fileNameSimulator = ssSimulator.str();
  getDataFromFile(fileNameSimulator, simulator/*out*/);

  // read the reference data of population A
  stringstream ssPopulationA;
  ssPopulationA << "chrom" << chromozomeNumber << "_win" << window << "_pop1.txt"; // population1 = reference data of population A
  string fileNamePopulationA = ssPopulationA.str();
  getDataFromFile(fileNamePopulationA, populationA/*out*/);

  // read the reference data of population B
  stringstream ssPopulationB;
  ssPopulationB << "chrom" << chromozomeNumber << "_win" << window << "_pop2.txt"; // population2 = reference data of population B
  string fileNamePopulationB = ssPopulationB.str();
  getDataFromFile(fileNamePopulationB, populationB/*out*/);
}

void calculateWindows(unsigned int chromozomeNumber, unsigned int numOfWindows, unsigned int numOfHaplotypes, double epsilon,
                      vector< vector< vector<double> > >& windowScore/*out*/, 
                      vector< vector< vector<bool> > >& maxArgIBDWindow/*out*/, 
                      vector< vector< vector<double> > >& windowDiff/*out*/)
{
  resizeVectors(numOfWindows, numOfHaplotypes, windowScore/*out*/, maxArgIBDWindow/*out*/, windowDiff/*out*/);

  for (unsigned int window = 0 ; window < numOfWindows ; window++)
  {
    // read the data of the current window (simulator data & reference data from both populations
    vector<string> simulator;
    vector<string> populationA;
    vector<string> populationB;
    getData(chromozomeNumber, window, simulator/*out*/, populationA/*out*/, populationB/*out*/);

    // build the LD trees according to the populations' data
    LDTree treeA = LDTree(populationA);
    LDTree treeB = LDTree(populationB);

    // calculate the probability for each haplotype from the simulator data in each populatio
    vector<double> simulatorProbabilityA(simulator.size(), 0);
    vector<double> simulatorProbabilityB(simulator.size(), 0);
    for (unsigned int i = 0 ; i < simulator.size() ; i++)
    {
      simulatorProbabilityA[i] = treeA.getHaplotypeProbability(simulator[i]);
      simulatorProbabilityB[i] = treeB.getHaplotypeProbability(simulator[i]);
    }

    // for each two persons, calculate the IBD probability and the nonIBD probability
    // according to these values, fill in the output vectors
    for (unsigned int i = 0 ; i < simulator.size() ; i++)
    {
      for (unsigned int j = i+1 ; j < simulator.size() ; j++)
      {
        double origIA = simulatorProbabilityA[i];
        double origIB = simulatorProbabilityB[i];
        double origJA = simulatorProbabilityA[j];
        double origJB = simulatorProbabilityB[j];
        // origIA + origIB should be 1 ==> normalize the values
        double normIA = origIA / (origIA + origIB);
        double normIB = 1 - normIA;
        // origJA + origJB should be 1 ==> normalize the values
        double normJA = origJA / (origJA + origJB);
        double normJB = 1 - normJA;

        double ibdA = treeA.getIBDProbability(simulator[i], simulator[j], epsilon);
        double ibdB = treeB.getIBDProbability(simulator[i], simulator[j], epsilon);
        double totalIBD = 
          ibdA * normIA * normJA + // IBD(Hi,Hj|A) * Pr(Hi beolongs to A) * Pr(Hj beolongs to A)
          ibdB * normIB * normJB; // IBD(Hi,Hj|B) * Pr(Hi beolongs to B) * Pr(Hj beolongs to B)
        double totalNonIBD = 
          simulatorProbabilityA[i] * simulatorProbabilityA[j] * normIA * normJA +  // nonIBD(Hi,Hj|A) * Pr(Hi beolongs to A) * Pr(Hj beolongs to A)
          simulatorProbabilityB[i] * simulatorProbabilityB[j] * normIB * normJB + // nonIBD(Hi,Hj|B) * Pr(Hi beolongs to B) * Pr(Hj beolongs to B)
          simulatorProbabilityA[i] * simulatorProbabilityB[j] * normIA * normJB + // nonIBD(Hi,Hj|HiA, HjB) * Pr(Hi beolongs to A) * Pr(Hj beolongs to B)
          simulatorProbabilityB[i] * simulatorProbabilityA[j] * normIB * normJA; // nonIBD(Hi,Hj|HiB, HjA) * Pr(Hi beolongs to B) * Pr(Hj beolongs to A)
        windowScore[window][i][j] = totalIBD / totalNonIBD; // the window score for haplotypes i & j
        maxArgIBDWindow[window][i][j] = (ibdA > ibdB); // in which population there is a better chance for an IBD - true==>A, false==>B
        windowDiff[window][i][j] = abs(ibdA - ibdB); // what is the difference between the IBD probability in the populations
      }
    }
  }
}

// returns true if the windows population is the same as the population (the maximum IBD result is in this population)
// or if the difference between the IBD results of the populations (winPopulation) are not bigger than winDiff
bool isIBD(double maxDiff, bool winPopulation, double winDiff, bool population)
{
  return ((winPopulation==population) || (winDiff <= maxDiff));
}

// check whether there is an IBD from startWin to endWin (not include) of hap1 and hap2 in a specific population
// population: true = population A, false = population B
// returns true if there is an IBD and false otherwise
bool checkIBDPopulation(unsigned int startWin, unsigned int endWin, unsigned int hap1, unsigned int hap2, double maxDiff,
                        vector< vector< vector<bool> > >& maxArgIBDWindow, 
                        vector< vector< vector<double> > >& windowDiff,
                        bool population)
{
  for (unsigned int win = startWin ; win < endWin ; win++)
  {
    if (!isIBD(maxDiff, maxArgIBDWindow[win][hap1][hap2], windowDiff[win][hap1][hap2], population))
    {
      return false;
    }
  }
  return true;
}

// check whether there is an IBD from startWin to endWin (not include) of hap1 and hap2 in one of the populations
// returns true if there is an IBD and false otherwise
bool checkIBD(unsigned int startWin, unsigned int endWin, unsigned int hap1, unsigned int hap2, double maxDiff,
              vector< vector< vector<bool> > >& maxArgIBDWindow, 
              vector< vector< vector<double> > >& windowDiff)
{
  if (checkIBDPopulation(startWin, endWin, hap1, hap2, maxDiff, maxArgIBDWindow, windowDiff, true))
  {
    // check IBD in population A
    return true;
  }
  if (checkIBDPopulation(startWin, endWin, hap1, hap2, maxDiff, maxArgIBDWindow, windowDiff, false))
  {
    // check IBD in population B
    return true;
  }
  return false;
}

// IBDhaplotypes is a vector which size is as the number of blocks
// in each place, there is a vector of haplotypes pairs = pairs of IBDs in this block
void calculateBlocks(unsigned int numOfWindows, unsigned int numOfHaplotypes, unsigned int blockSize, double threshold, double maxDiff,
                     vector< vector< vector<double> > >& windowScore, 
                     vector< vector< vector<bool> > >& maxArgIBDWindow, 
                     vector< vector< vector<double> > >& windowDiff,
                     vector< vector < pair < unsigned int, unsigned int> > >& IBDHaplotypes/*out*/)
{
  unsigned int numOfBlocks = numOfWindows - blockSize + 1;
  IBDHaplotypes.resize(numOfBlocks);

  // go over all the blocks and over all the haplotypes and 
  // fill the output vector with pairs of haplotypes that are IBD in that block
  for (unsigned int blockNum = 0 ; blockNum < numOfBlocks ; blockNum++)
  {
    // the block is from startWin till endWin (not include)
    unsigned int startWin = blockNum;
    unsigned int endWin = blockNum + blockSize;

    for (unsigned int i = 0 ; i < numOfHaplotypes ; i++)
    {
      for (unsigned int j = i+1 ; j < numOfHaplotypes ; j++)
      {
        // the block score is a product of the windows' scores of these two haplotypes (i,j) 
        double blockScore = 1;
        for (unsigned int win = startWin ; win < endWin ; win++)
        {
          blockScore *= windowScore[blockNum][i][j];
        }

        // if the blockScore is big enough (bigger than the threshold) and there is an IBD series of the same population,
        // i and j are IBD in this block ==> they are added to IBDHaplotypes[blockNum]
        if ((blockScore >= threshold) && (checkIBD(startWin, endWin, i, j, maxDiff, maxArgIBDWindow, windowDiff)))
        {
          IBDHaplotypes[blockNum].push_back(make_pair(i, j));
        }
      }
    }
  }
}

void createResultFile(unsigned int chromozomeNumber,
                      vector< vector < pair < unsigned int, unsigned int> > >& IBDHaplotypes)
{
  stringstream ss;
  ss << "chrom" << chromozomeNumber << "_output.txt";
  string fileName = ss.str();

  ofstream file;
  file.open (fileName.c_str());
  for (vector< vector < pair < unsigned int, unsigned int> > >::const_iterator blockIter = IBDHaplotypes.begin(); blockIter != IBDHaplotypes.end(); blockIter ++)
  {
    for (vector < pair < unsigned int, unsigned int> >::const_iterator pairIter = (*blockIter).begin();
      pairIter != (*blockIter).end(); pairIter++)
    {
      file << (*pairIter).first << "-" << (*pairIter).second << "\t";
    }
    file << "\n";
  }
  file.close();
}

// the main parameters:
// 1. chromozome number
// 2. number of windows in the chromozome
// 3. number of haplotypes
// 4. epsilon for the IBD calculations in the LD trees
// 5. block size - number of windows in a block
// 6. threshold for the block score
// 7. the maximum difference between IBD in two populations
int main(int argc, char* argv[])
{
  // argv[0] is the path and name of the program itself
  unsigned int chromozomeNumber = atoi(argv[1]);
  unsigned int numOfWindows = atoi(argv[2]);
  unsigned int numOfHaplotypes = atoi(argv[3]);
  double epsilon = atof(argv[4]);
  unsigned int blockSize = atoi(argv[5]);
  double threshold = atof(argv[6]);
  double maxDiff = atof(argv[7]);

  // the external vector goes over all the windows
  // the interior vectors go over every two persons
  // for example: windowScore[a][b][c] is the window score of window a for persons b and c
  vector< vector< vector<double> > > windowScore; // the windowScore for each two persons
  vector< vector< vector<bool> > > maxArgIBDWindow; // the bigger IBD result for each two persons - true for population A, false for population B
  vector< vector< vector<double> > > windowDiff; // the difference between the IBD in population A and the IBD in population B for each two person
  
  // fill the above vectors for each window and for each two persons
  calculateWindows(chromozomeNumber, numOfWindows, numOfHaplotypes, epsilon,
    windowScore/*out*/, maxArgIBDWindow/*out*/, windowDiff/*out*/);

  // the external vector goes over all the blocks
  // the interior vector will be filled with pairs of haplotypes that considered as IBD in this block
  // for example: IBDHaplotypes[a] is a vector of pairs of haplotypes numbers <x,y> such that x-y are considered to have an IBD in block a
  vector< vector < pair < unsigned int, unsigned int> > > IBDHaplotypes;

  // fill the above vector according to the windows' calculations
  calculateBlocks(numOfWindows, numOfHaplotypes, blockSize, threshold, maxDiff, windowScore, maxArgIBDWindow, windowDiff, IBDHaplotypes/*out*/);

  // create an ouput file with the name "chromX_output" such that X is the chromozome number
  // the file is in the following format:
  // each line indicates pairs of IBD haplotypes in the coordinate block
  // each pair is written as "hap-hap" and there are tabs between the pairs
  createResultFile(chromozomeNumber, IBDHaplotypes);

  return 1;
}