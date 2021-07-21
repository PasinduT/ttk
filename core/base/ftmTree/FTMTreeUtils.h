/// \ingroup base
/// \class FTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_H
#define _FTMTREEUTILS_H

#pragma once

#include <FTMTree_MT.h>
#include <stack>

using namespace ttk;
using namespace ftm;

// --------------------
// Utils
// --------------------

void myPause();

void printTreesStats(std::vector<ftm::FTMTree_MT *> &trees);

template <class dataType>
void printTree(MergeTree<dataType> &tree, bool doPrint = true) {
  tree.tree.printTree(doPrint);
}

template <class dataType>
void printTreeStats(MergeTree<dataType> &tree) {
  tree.tree.printTreeStats();
}

template <class dataType>
void printTreeScalars(MergeTree<dataType> &tree, bool printNodeAlone = true) {
  printTreeScalars<dataType>(&(tree.tree), printNodeAlone);
}

template <class dataType>
std::vector<ftm::FTMTree_MT *>
  mergeTreeToFTMTree(std::vector<MergeTree<dataType>> &trees) {
  std::vector<ftm::FTMTree_MT *> treesT;
  for(MergeTree<dataType> &t : trees)
    treesT.push_back(&(t.tree));
  return treesT;
}

template <class dataType>
MergeTree<double> mergeTreeTemplateToDouble(MergeTree<dataType> &mt) {
  std::vector<double> newScalarsValues;
  for(auto val : mt.scalarsValues)
    newScalarsValues.push_back(static_cast<double>(val));
  MergeTree<double> newMt
    = MergeTree<double>(mt.scalars, newScalarsValues, mt.params);
  newMt.tree.copyMergeTreeStructure(&(mt.tree));
  return newMt;
}

template <class dataType>
std::vector<MergeTree<double>>
  mergeTreesTemplateToDouble(std::vector<MergeTree<dataType>> &mts) {
  std::vector<MergeTree<double>> newMts;
  for(auto &mt : mts)
    newMts.push_back(mergeTreeTemplateToDouble<dataType>(mt));
  return newMts;
}

template <class dataType>
MergeTree<dataType> mergeTreeDoubleToTemplate(MergeTree<double> &mt) {
  std::vector<dataType> newScalarsValues;
  for(auto val : mt.scalarsValues)
    newScalarsValues.push_back(static_cast<dataType>(val));
  MergeTree<dataType> newMt
    = MergeTree<dataType>(mt.scalars, newScalarsValues, mt.params);
  newMt.tree.copyMergeTreeStructure(&(mt.tree));
  return newMt;
}

template <class dataType>
std::vector<MergeTree<dataType>>
  mergeTreesDoubleToTemplate(std::vector<MergeTree<double>> &mts) {
  std::vector<MergeTree<dataType>> newMts;
  for(auto &mt : mts)
    newMts.push_back(mergeTreeDoubleToTemplate<dataType>(mt));
  return newMts;
}

// --------------------
// MergeTree
// --------------------

template <class dataType>
MergeTree<dataType> createEmptyMergeTree(int scalarSize) {
  // Init Scalars
  ftm::Scalars scalars;
  scalars.size = scalarSize;
  dataType *scalarsValues = nullptr;
  scalars.values = (void *)scalarsValues;

  // Init Params
  ftm::Params params;
  params.treeType = ftm::Join_Split;

  // Init tree
  MergeTree<dataType> mergeTree(scalars, params);

  return mergeTree;
}

template <class dataType>
void setTreeScalars(MergeTree<dataType> &mergeTree,
                    std::vector<dataType> &scalarsVector) {
  mergeTree.scalarsValues = scalarsVector;
  mergeTree.scalars.values = (void *)(mergeTree.scalarsValues.data());
  mergeTree.scalars.size = mergeTree.scalarsValues.size();
}

template <class dataType>
std::vector<dataType> getTreeScalars(ftm::FTMTree_MT *tree) {
  std::vector<dataType> scalarsVector;
  for(unsigned int i = 0; i < tree->getNumberOfNodes(); ++i)
    scalarsVector.push_back(tree->getValue<dataType>(i));
  return scalarsVector;
}

template <class dataType>
std::vector<dataType> getTreeScalars(MergeTree<dataType> &mergeTree) {
  return getTreeScalars<dataType>(&(mergeTree.tree));
}

template <class dataType>
MergeTree<dataType> copyMergeTree(ftm::FTMTree_MT *tree,
                                  bool doSplitMultiPersPairs = false) {
  std::vector<dataType> scalarsVector = getTreeScalars<dataType>(tree);

  // Get multi persistence pairs
  std::vector<ftm::idNode> multiPersOrigins;
  if(doSplitMultiPersPairs) {
    multiPersOrigins = tree->getMultiPersOrigins<dataType>(true);
    for(ftm::idNode nodeOrigin : multiPersOrigins) {
      scalarsVector[nodeOrigin]
        = tree->getValue<dataType>(tree->getNode(nodeOrigin)->getOrigin());
      scalarsVector.push_back(tree->getValue<dataType>(nodeOrigin));
    }
  }

  // Create new tree
  MergeTree<dataType> mergeTree
    = createEmptyMergeTree<dataType>(scalarsVector.size());
  ftm::FTMTree_MT *treeNew = &(mergeTree.tree);
  setTreeScalars<dataType>(mergeTree, scalarsVector);

  // Copy tree structure
  treeNew->copyMergeTreeStructure(tree);

  // Add multi persistence nodes origins
  if(doSplitMultiPersPairs) {
    for(ftm::idNode nodeOrigin : multiPersOrigins) {
      int nodeCpt = treeNew->getNumberOfNodes();
      treeNew->makeNode(nodeCpt);
      treeNew->getNode(nodeCpt)->setOrigin(nodeOrigin);
      treeNew->getNode(nodeOrigin)->setOrigin(nodeCpt);
    }
  }

  return mergeTree;
}

template <class dataType>
MergeTree<dataType> copyMergeTree(MergeTree<dataType> &mergeTree,
                                  bool doSplitMultiPersPairs = false) {
  return copyMergeTree<dataType>(&(mergeTree.tree), doSplitMultiPersPairs);
}

// Remove unused nodes
template <class dataType>
MergeTree<dataType> cleanMergeTree(ftm::FTMTree_MT *tree,
                                   std::vector<int> &nodeCorr,
                                   bool useBD = true) {
  int verboseT = 0;

  // Create new tree
  if(verboseT > 0)
    std::cout << "// Create new tree" << std::endl;
  int newNoNodes = tree->getRealNumberOfNodes() * 2;
  MergeTree<dataType> mTreeNew = createEmptyMergeTree<dataType>(newNoNodes);
  ftm::FTMTree_MT *treeNew = &(mTreeNew.tree);
  std::vector<dataType> newScalarsVector(newNoNodes, 0);

  // Copy the old tree structure
  if(verboseT > 0)
    std::cout << "// Copy the old tree structure" << std::endl;
  std::vector<int> nodeDone(tree->getNumberOfNodes(), 0);
  nodeCorr = std::vector<int>(tree->getNumberOfNodes(), -1);
  std::vector<std::vector<ftm::idNode>> treeMultiPers;
  if(not useBD)
    treeMultiPers = tree->getMultiPersOriginsVectorFromTree();
  std::queue<ftm::idNode> queue;
  for(auto leaf : tree->getLeavesFromTree())
    queue.push(leaf);
  while(!queue.empty()) {
    ftm::idNode node = queue.front();
    queue.pop();
    ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();

    if(useBD) {
      int nodeOriginIndex = treeNew->getNumberOfNodes();
      if(nodeCorr[nodeOrigin] == -1)
        treeNew->makeNode(nodeOriginIndex);
      else
        nodeOriginIndex = nodeCorr[nodeOrigin];
      int nodeIndex = treeNew->getNumberOfNodes();
      if(nodeCorr[node] == -1)
        treeNew->makeNode(nodeIndex);
      else
        nodeIndex = nodeCorr[node];

      if(nodeCorr[nodeOrigin] == -1)
        treeNew->getNode(nodeOriginIndex)->setOrigin(nodeIndex);
      treeNew->getNode(nodeIndex)->setOrigin(nodeOriginIndex);

      newScalarsVector[nodeOriginIndex] = tree->getValue<dataType>(nodeOrigin);
      newScalarsVector[nodeIndex] = tree->getValue<dataType>(node);
      nodeCorr[nodeOrigin] = nodeOriginIndex;
      nodeCorr[node] = nodeIndex;
    } else {
      int nodeCpt = treeNew->getNumberOfNodes();
      treeNew->makeNode(nodeCpt);
      if(!tree->isLeaf(node)) {
        treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
        if(not(tree->isRoot(node) and node == nodeOrigin))
          treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
        for(auto nodeMultiPers : treeMultiPers[node])
          treeNew->getNode(nodeCorr[nodeMultiPers])->setOrigin(nodeCpt);
      } else if(tree->isNodeAlone(nodeOrigin)) { // saddle merged
        treeNew->makeNode(nodeCpt + 1);
        newScalarsVector[nodeCpt + 1] = tree->getValue<dataType>(nodeOrigin);
        nodeCorr[nodeOrigin] = nodeCpt + 1;
        treeNew->getNode(nodeCpt)->setOrigin(nodeCorr[nodeOrigin]);
        treeNew->getNode(nodeCorr[nodeOrigin])->setOrigin(nodeCpt);
      }
      newScalarsVector[nodeCpt] = tree->getValue<dataType>(node);
      nodeCorr[node] = nodeCpt;
    }

    for(auto child : tree->getChildren(node))
      treeNew->makeSuperArc(nodeCorr[child], nodeCorr[node]);

    if(!tree->isRoot(node)) {
      ftm::idNode parent = tree->getParentSafe(node);
      nodeDone[parent] += 1;
      if(nodeDone[parent] == tree->getNumberOfChildren(parent))
        queue.push(parent);
    }
  }

  // Set new scalars
  if(verboseT > 0)
    std::cout << "// Set new scalars" << std::endl;
  setTreeScalars<dataType>(mTreeNew, newScalarsVector);

  // Manage full merge
  if(verboseT > 0)
    std::cout << "// Manage full merge" << std::endl;
  auto treeRoot = tree->getRoot();
  if(tree->getNode(treeRoot)->getOrigin() == (int)treeRoot) {
    auto treeNewRoot = treeNew->getRoot();
    // auto treeNewRootOrigin = tree->getNode(treeNewRoot)->getOrigin();
    // treeNew->getNode(treeNewRootOrigin)->setOrigin(treeNewRootOrigin);
    treeNew->getNode(treeNewRoot)->setOrigin(treeNewRoot);
  }

  // Return new tree
  if(verboseT > 0)
    std::cout << "// Return new tree" << std::endl;
  return mTreeNew;
}

template <class dataType>
void cleanMergeTree(MergeTree<dataType> &mTree,
                    std::vector<int> &nodeCorr,
                    bool useBD = true) {
  mTree = cleanMergeTree<dataType>(&(mTree.tree), nodeCorr, useBD);
}

template <class dataType>
void cleanMergeTree(MergeTree<dataType> &mTree, bool useBD = true) {
  std::vector<int> nodeCorr;
  cleanMergeTree<dataType>(mTree, nodeCorr, useBD);
}

#endif
