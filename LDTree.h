#ifndef _LDTREE_H
#define _LDTREE_H

#include <set>
#include <vector>
#include <string>

using namespace std;

class LDTree
{
public:
  LDTree(const vector<string>& haplotypes, bool debug);
  ~LDTree();
  double getHaplotypeProbability(const string& haplotype);
  double getIBDProbability(const string& haplotype1, const string& haplotype2, double epsilon);
  void printLDTree();

private:
  class LDNode
  {
  public: 
    LDNode(LDNode* parent);
    ~LDNode();
    LDNode* childrenPtr[2];
    int childrenNum[2];
    set<string> possibleStrings;
    LDNode* originParent;
    LDNode* nextInDepth;
    LDNode* prevInDepth;
  };

  LDNode* root;
  vector<LDNode*> depthArray;

  // ctor - help functions
  // 1. build the tree
  void createChain(LDNode* node, string haplotype, int depth);
  void insertNodeToArray(LDNode* node, int depth);
  LDNode* createNode(LDNode* parent, int depth);
  void buildTree(const vector<string>& haplotypes);
  // 2. squeeze common parts of the tree
  void squeezeNumbers(LDNode* node, LDNode* temp);
  void squeezeNodes(LDNode* node, LDNode* temp);
  void squeezeLevel(int depth);
  void squeezeTree();

  // dtor - help functions
  void nullDepthChildren(LDNode* node, int currentChild);
  void updateDepthArray(LDNode* node, int depth);
  void freeNodes(LDNode* node, int depth);

  // getHaplotypeProbability + getIBDProbability - help function
  int getChildInfo(const LDNode* const parent, const char& SNP, LDNode** child); // fill the child ptr, return the child number

  // printLDTree - help functions
  void printRecursiveNodes(LDNode* node);
  void printNode(LDNode* node);
};

#endif