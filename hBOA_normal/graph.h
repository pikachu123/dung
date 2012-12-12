#ifndef _graph_h_
#define _graph_h_

#include <stdio.h>

#define CONNECTED      1
#define NOT_CONNECTED  0

#define SUCCESS        1
#define FAIL           0

//#define REDUCTION_STRATEGY_1
#define REDUCTION_STRATEGY_2

//===============================================================

class OrientedGraph {

 private:

  int  N;              // the number of vertices
  long N2;             // N square
  char **coincidence;  // the coincidence matrix
  char **path;         // the matrix for maintanance of paths
  int  *numIn;         // the number of incoming vertices
  int  *numOut;        // the number of outcoming vertices
  int  *mark;          // the array for vertex-marks
  int  **parentList;    // the list of parents for each node

 public:

  OrientedGraph(int n);
  ~OrientedGraph();

  int size();

  int addEdge(int i, int j);
  int removeEdge(int i, int j);
  int reverseEdge(int i, int j);
  int removeAllEdges();
  int setMark(int i, int val);
  int setAllMarks(int val);
  int removeMark(int i);
  int removeAllMarks();

  int getNumberOfVertices();
  int connected(int i, int j);
  int notConnected(int i, int j);
  int existsPath(int i, int j);
  int getNumIn(int i);
  int getNumOut(int i);
  int getMark(int i);
  int *getParentList(int i);
  int canAddEdge(int i, int j);
  int canReverseEdge(int i, int j);

  int printCoincidenceMatrix(FILE *out);
  int printPathMatrix(FILE *out);
  int printNumInArray(FILE *out);
  int printNumOutArray(FILE *out);
};

//===============================================================

class AcyclicOrientedGraph: public OrientedGraph {

 public:

  AcyclicOrientedGraph(int n);
  ~AcyclicOrientedGraph();

  int addEdge(int i, int j);
  int reverseEdge(int i, int j);
  int canAddEdge(int i, int j);
  int canReverseEdge(int i, int j);
};

class LoudAcyclicOrientedGraph:public AcyclicOrientedGraph {

 public:

  LoudAcyclicOrientedGraph(int n);
  ~LoudAcyclicOrientedGraph();

  int addEdge(int i, int j);
  int removeEdge(int i, int j);
  int reverseEdge(int i, int j);
  int removeAllEdges();

  int connected(int i, int j);
  int notConnected(int i, int j);
  int existsPath(int i, int j);
};

//===============================================================

class LeastIncommingEdgesAcyclicOrientedGraph:public AcyclicOrientedGraph
{

 public:

  LeastIncommingEdgesAcyclicOrientedGraph(int n);
  ~LeastIncommingEdgesAcyclicOrientedGraph();

  int addUnorientedEdge(int i, int j);
};

//===============================================================

class BoundedIncommingEdgesAcyclicOrientedGraph:public AcyclicOrientedGraph
{
 private:

  int K;

  int reduce(int i);

#ifdef REDUCTION_STRATEGY_2
  int findBestReduction(int i, int *reducedTo);
  int applyBestReduction(int i);
#endif

 public:
  
  BoundedIncommingEdgesAcyclicOrientedGraph(int n, int k);
  ~BoundedIncommingEdgesAcyclicOrientedGraph();

  int addUnorientedEdge(int i, int j);
};


#endif
