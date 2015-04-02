#include "LDTree.h"
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <iostream>

void LDTree::insertNodeToArray(LDNode* node, unsigned int depth)
{
  if (this->depthArray[depth]!=NULL)
  {
    LDNode* currentNode = this->depthArray[depth];
    currentNode->prevInDepth = node;
    node->nextInDepth = currentNode;
  }
  this->depthArray[depth] = node;
}

LDTree::LDNode::LDNode(LDNode* parent)
{
  this->childrenPtr[0] = this->childrenPtr[1] = NULL;
  this->childrenNum[0] = this->childrenNum[1] = 0;
  this->originParent = parent;
  this->nextInDepth = this->prevInDepth = NULL;
  // the possibleStrings set is empty
}

LDTree::LDNode::~LDNode()
{
  this->possibleStrings.clear(); // as the deafult dtor but more clear like this
}

LDTree::LDNode* LDTree::createNode(LDNode* parent, unsigned int depth)
{
  LDNode* node = new LDNode(parent);
  insertNodeToArray(node, depth);
  return node;
}

void LDTree::createChain(LDNode* node, string haplotype, unsigned int depth)
{
  node->possibleStrings.insert(haplotype);
  LDNode* child;
  bool isLastSNP = (haplotype.size()==1);
  string firstSNP = isLastSNP ? haplotype : haplotype.substr(0,1);
  if (firstSNP=="0")
  {
    if ((node->childrenPtr[0]==NULL) && (!isLastSNP))
    {
      node->childrenPtr[0] = createNode(node, depth+1);
    }
    child = node->childrenPtr[0];
    node->childrenNum[0]++;
  }
  else
  {
    assert(firstSNP=="1");
    if ((node->childrenPtr[1]==NULL) && (!isLastSNP))
    {
      node->childrenPtr[1] = createNode(node, depth+1);
    }
    child = node->childrenPtr[1];
    node->childrenNum[1]++;
  }
  if (!isLastSNP)
  {
    createChain(child, haplotype.substr(1), depth+1);
  }
}

void LDTree::buildTree(const vector<string>& haplotypes)
{
  unsigned int hapLength  = haplotypes[0].size();
  this->depthArray.resize(hapLength, NULL);
  this->root = createNode(NULL, 0);
  for (vector<string>::const_iterator iter = haplotypes.begin(); iter!=haplotypes.end(); iter++)
  {
    assert((*iter).size()==hapLength);
    createChain(this->root, *iter, 0);
  }
}

// squeeze the children numbers of temp into node
void LDTree::squeezeNumbers(LDNode* node, LDNode* temp)
{
  node->childrenNum[0] += temp->childrenNum[0];
  node->childrenNum[1] += temp->childrenNum[1];
  if (node->childrenPtr[0] != NULL)
  {
    assert(temp->childrenPtr[0]!=NULL);
    squeezeNumbers(node->childrenPtr[0], temp->childrenPtr[0]);
  }
  if (node->childrenPtr[1] != NULL)
  {
    assert(temp->childrenPtr[1]!=NULL);
    squeezeNumbers(node->childrenPtr[1], temp->childrenPtr[1]);
  }
}

// squeeze temp into node
void LDTree::squeezeNodes(LDNode* node, LDNode* temp)
{
  // update temp parent
  // after the squeeze, node has multiple parents (it keep saving it's origin one)
  LDNode* tempParent = temp->originParent; // can't be NULL since there temp isn't in depth==0
  if (tempParent->childrenPtr[0] == temp)
  {
    tempParent->childrenPtr[0] = node;
  }
  else
  {
    assert(tempParent->childrenPtr[1] == temp);
    tempParent->childrenPtr[1] = node;
  }

  // update the childrenNum of node according to temp children
  squeezeNumbers(node, temp);
}

void LDTree::squeezeLevel(unsigned int depth)
{
  LDNode* node = depthArray[depth];
  LDNode* next;
  while (node != NULL)
  {
    // for each node, go over all the nodes
    // and check whether they can be squeezed
    next = node->nextInDepth;
    while (next != NULL)
    {
      if (node->possibleStrings == next->possibleStrings)
      {
        LDNode* temp = next;
        LDNode* prev = temp->prevInDepth;
        next = temp->nextInDepth;
        // since we squeeze next into node, next->prevInDepth can't be NULL
        prev->nextInDepth = next;
        if (next!=NULL)
        {
          next->prevInDepth = prev;
        }
        squeezeNodes(node, temp); // squeeze temp into node
        freeNodes(temp, depth); // free the squeezed tree
      }
      else
      {
        next = next->nextInDepth;
      }
    }
    node = node->nextInDepth;
  }
}

void LDTree::squeezeTree()
{
  // starting from depth=1 since depth==0 contains only the root (no squeeze potential)
  for (unsigned int depth = 1; depth < this->depthArray.size(); depth++)
  {
    assert(this->depthArray[depth]!=NULL);
    squeezeLevel(depth);
  }
}

LDTree::LDTree(const vector<string>& haplotypes, bool debug)
{
  if (debug)
  {
    cout << "build the tree" << endl;
  }
  buildTree(haplotypes);

  if (debug)
  {
    cout << "squeeze the tree" << endl;
  }
  squeezeTree();
}

void LDTree::nullDepthChildren(LDNode* node, unsigned int currentChild)
{
  if ((currentChild==0) && (node->childrenPtr[0]==node->childrenPtr[1]))
  {
    // if checking the first child (0), need to check if its child is the same as the second child (1)
    node->childrenPtr[1] = NULL;
  }
  LDNode* child = node->childrenPtr[currentChild];
  LDNode* next = node->nextInDepth;
  while (next != NULL)
  {
    if (next->childrenPtr[0]==child)
    {
      next->childrenPtr[0] = NULL;
    }
    if (next->childrenPtr[1]==child)
    {
      next->childrenPtr[1] = NULL;
    }
    next = next->nextInDepth;
  }
  LDNode* prev = node->prevInDepth;
  while (prev != NULL)
  {
    if (prev->childrenPtr[0]==child)
    {
      prev->childrenPtr[0] = NULL;
    }
    if (prev->childrenPtr[1]==child)
    {
      prev->childrenPtr[1] = NULL;
    }
    prev = prev->prevInDepth;
  }
}

void LDTree::updateDepthArray(LDNode* node, unsigned int depth)
{
  LDNode* prev = node->prevInDepth;
  LDNode* next = node->nextInDepth;
  if (next!=NULL)
  {
    next->prevInDepth = prev;
  }
  if (prev!=NULL)
  {
    prev->nextInDepth = next;
  }
  if (this->depthArray[depth]==node)
  {
    assert(prev==NULL); // next can be NULL if it is been called from the dtor
    this->depthArray[depth] = next;
  }
}

void LDTree::freeNodes(LDNode* node, unsigned int depth)
{
  if (node->childrenPtr[0] != NULL)
  {
    nullDepthChildren(node, 0);
    freeNodes(node->childrenPtr[0], depth+1);
  }
  if (node->childrenPtr[1] != NULL)
  {
    nullDepthChildren(node, 1);
    freeNodes(node->childrenPtr[1], depth+1);
  }
  node->possibleStrings.clear();
  updateDepthArray(node, depth);
  delete node;
}

LDTree::~LDTree()
{
  freeNodes(this->root, 0);
  this->depthArray.clear();
}

unsigned int LDTree::getChildInfo(const LDNode* const parent, const char& SNP, LDNode** child)
{
  if (SNP=='0')
  {
    *child = parent->childrenPtr[0];
    return parent->childrenNum[0];
  }
  else
  {
    assert(SNP=='1');
    *child = parent->childrenPtr[1];
    return parent->childrenNum[1];
  }
}

double LDTree::getHaplotypeProbability(const string& haplotype, double forgivenessPercent)
{
  assert(this->depthArray.size()==haplotype.size());
  int maxWrongSNP = (int)floor(haplotype.size() * forgivenessPercent / 100);
  double res = 1;
  LDNode* node = this->root;
  LDNode* child;
  unsigned int childNum;
  int wrongCount = 0;
  for (unsigned int i = 0 ; i < haplotype.size(); i++)
  {
    char SNP = haplotype[i];
    childNum = getChildInfo(node, SNP, &child);
    if (childNum == 0)
    {
      if (wrongCount == maxWrongSNP)
      {
        // such haplotype doesn't exist according to the reference data
        // probably its probability is quite low
        return 0;
      }
      else
      {
        // assume the SNP is wrong --> take the other child
        SNP = (SNP == '1') ? '0' : '1';
        childNum = getChildInfo(node, SNP, &child);
        assert(childNum!=0); // it isn't possible that both children of node will be NULL
        wrongCount++;
      }
    }
    res *= (double)childNum / (double)(node->childrenNum[0] + node->childrenNum[1]);
    node = child;
  }
  return res;
}

double LDTree::getIBDProbability(const string& haplotype1, const string& haplotype2, double epsilon, double forgivenessPercent)
{
  assert((this->depthArray.size()==haplotype1.size()) && (haplotype1.size()==haplotype2.size()));
  int maxWrongSNP = (int)floor(haplotype1.size() * forgivenessPercent / 100);
  double res = 1;
  LDNode* node1 = this->root;
  LDNode* node2 = this->root;
  LDNode* child1;
  LDNode* child2;
  unsigned int child1Num, child2Num;
  double child1Prob, child2Prob;
  int wrongCount1 = 0;
  int wrongCount2 = 0;
  for (unsigned int i = 0; i < haplotype1.size(); i++)
  {
    char SNP1 = haplotype1[i];
    char SNP2 = haplotype2[i];
    child1Num = getChildInfo(node1, SNP1, &child1);
    child2Num = getChildInfo(node2, SNP2, &child2);
    if (child1Num==0)
    {
      if (wrongCount1 == maxWrongSNP)
      {
        // such haplotype doesn't exist according to the reference data
        // probably its probability is quite low
        return 0;
      }
      else
      {
        // assume the SNP is wrong --> take the other child
        SNP1 = (SNP1 == '1') ? '0' : '1';
        child1Num = getChildInfo(node1, SNP1, &child1);
        assert(child1Num!=0); // it isn't possible that both children of node will be NULL
        wrongCount1++;
      }
    }
    if (child2Num==0)
    {
      if (wrongCount2 == maxWrongSNP)
      {
        // such haplotype doesn't exist according to the reference data
        // probably its probability is quite low
        return 0;
      }
      else
      {
        // assume the SNP is wrong --> take the other child
        SNP2 = (SNP2 == '1') ? '0' : '1';
        child2Num = getChildInfo(node2, SNP2, &child2);
        assert(child2Num!=0); // it isn't possible that both children of node will be NULL
        wrongCount2++;
      }
    }
    child1Prob = (double)child1Num / (double)(node1->childrenNum[0] + node1->childrenNum[1]);
    child2Prob = (double)child2Num / (double)(node2->childrenNum[0] + node2->childrenNum[1]);
    res *= min(child1Prob, child2Prob);
    res *= ((SNP1==SNP2) ? (1-epsilon) : epsilon);
    node1 = child1;
    node2 = child2;
  }
  return res;
}

void LDTree::printNode(LDNode* node)
{
  cout << "printNode: " << node << endl;
  if (node->childrenPtr[0]==NULL)
  {
    cout << "first child - doesn't exist" << endl;
  }
  else
  {
    cout << "first child - " << node->childrenPtr[0] << endl;
  }
  if (node->childrenPtr[1]==NULL)
  {
    cout << "second child - doesn't exist" << endl;
  }
  else
  {
    cout << "second child - " << node->childrenPtr[1] << endl;
  }
  cout << "first child num = " << node->childrenNum[0]
    << " second child num = " << node->childrenNum[1] << endl;
  if (node->originParent==NULL)
  {
    cout << "origin parent - doesn't exist" << endl;
  }
  else
  {
    cout << "origin parent - " << node->originParent << endl;
  }
  if (node->nextInDepth==NULL)
  {
    cout << "next in depth - doesn't exist" << endl;
  }
  else
  {
    cout << "next in depth - " << node->nextInDepth << endl;
  }
  if (node->prevInDepth==NULL)
  {
    cout << "previous in depth - doesn't exist" << endl;
  }
  else
  {
    cout << "previous in depth - " << node->prevInDepth << endl;
  }
  cout << "possible strings:";
  for (set<string>::const_iterator iter = node->possibleStrings.begin();
    iter != node->possibleStrings.end(); iter++)
  {
    cout << "\t" << *iter;
  }
  cout << endl;
  cout << "done: " << node << endl;
}

void LDTree::printRecursiveNodes(LDNode* node)
{
  printNode(node);
  if (node->childrenPtr[0]!=NULL)
  {
    printRecursiveNodes(node->childrenPtr[0]);
  }
  if (node->childrenPtr[1]!=NULL)
  {
    printRecursiveNodes(node->childrenPtr[1]);
  }
}

void LDTree::printLDTree()
{
  printRecursiveNodes(this->root);
}