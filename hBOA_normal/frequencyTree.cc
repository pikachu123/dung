#include "frequencyTree.h"

#include "memalloc.h"
#include "population.h"

#define SPLIT 0
#define NODE  1

#define allocateNewNode(x) { (x) = (TreeNode*) Malloc(sizeof(TreeNode)); (x)->which = NODE; (x)->left = (x)->right = NULL; (x)->value = 0; }
#define allocateNewSplit(x) { (x) = (TreeSplit*) Malloc(sizeof(TreeSplit)); ((TreeSplit*) x)->which = SPLIT; ((TreeSplit*) x)->left = ((TreeSplit*) x)->right = NULL; ((TreeSplit*) x)->value = 0; ((TreeSplit*) x)->numSubtrees = 0; ((TreeSplit*) x)->subtree = NULL; numLeaves++; }
//#define allocateNewLeaf(x) { (x) = (TreeLeaf*) Malloc(sizeof(TreeLeaf)); ((TreeLeaf*) x)->left = ((TreeLeaf*) x)->right = NULL; ((TreeLeaf*) x)->list = NULL; (x)->value = 0; numLeaves++; }
#define isLeaf(x)          ((x->left==NULL)&&(x->right==NULL))
#define createEmptyTree    { allocateNewNode(root); numLeaves = 0; listSize = 0; }
#define deleteSubtrees(x)  { deleteSubtree(x->left); deleteSubtree(x->right); x->left = x->right = NULL; }

// ====================================================================

FrequencyTree::FrequencyTree()
{
    // the tree is empty 

    createEmptyTree;
};

// ------------------------------------------------

FrequencyTree::FrequencyTree(const FrequencyTree& frequencyTree)
{
    // copy the variables first

    numLeaves = frequencyTree.numLeaves;
    n         = frequencyTree.n;

    // copy the tree starting in the root of the other tree

    root = copySubtree(frequencyTree.root);
};

// ------------------------------------------------

FrequencyTree::FrequencyTree(const FrequencyTree& frequencyTree, int *index, int indexLength)
{
    // create an empty tree

    createEmptyTree;

  // number of leaves (excluding empty root) is 0

    numLeaves = 0;

    // copy the tree starting in the root of the other tree

    reduceByIndex(frequencyTree.root,index,indexLength);
};

// ------------------------------------------------

TreeNode *FrequencyTree::copySubtree(TreeNode *node)
{
    fprintf(stderr,"ERROR: copy subtree doesn't really work.");

    if (node)
	{
	    TreeNode *currentNode;

	    if ((isLeaf(node))&&(node!=root))
		{
		    if (listSize==0)
			{
			    allocateNewNode(currentNode); 
			    currentNode->value = node->value;
			}
		    else
			{
			    int i;

			    allocateNewSplit(currentNode);
	      
			    ((TreeSplit*)currentNode)->value       = ((TreeSplit*)node)->value;
			    ((TreeSplit*)currentNode)->numSubtrees = ((TreeSplit*)node)->numSubtrees;
			    ((TreeSplit*)currentNode)->subtree     = (FrequencyTree**) Calloc(((TreeSplit*)currentNode)->numSubtrees,sizeof(FrequencyTree*));
			    for (i=0; i<((TreeSplit*)currentNode)->numSubtrees; i++)
				((TreeSplit*)currentNode)->subtree[i] = new FrequencyTree(*((TreeSplit*)node)->subtree[i]);
			}
		}
	    else
		{
		    allocateNewNode(currentNode);
	  
		    currentNode->value = node->value;
		    currentNode->left  = copySubtree(node->left);
		    currentNode->right = copySubtree(node->right);
		};

	    return currentNode;
	}
    else
	return NULL;
};

// ------------------------------------------------

FrequencyTree::~FrequencyTree()
{
    // delete the tree

    deleteSubtree(root);

};

// ------------------------------------------------

int FrequencyTree::reset()
{
    // is the tree non-empty? if yes, delete the tree starting in the root

    deleteSubtree(root);

  // the tree is empty
  
    createEmptyTree;

    // no leaves anymore

    numLeaves = 0;

    // get back

    return 0;
};

// ------------------------------------------------

int FrequencyTree::deleteSubtree(TreeNode *node)
{
    if (node)
	{
	    if (isLeaf(node))
		{
		    if (listSize>0)
			{
			    int i;

			    for (i=0; i<((TreeSplit*)node)->numSubtrees; i++)
				delete ((TreeSplit*) node)->subtree[i];

			    Free(((TreeSplit*)node)->subtree);
			}
		}
	    else
		{
		    // is there a left subtree? if yes, delete it
	  
		    deleteSubtree(node->left);
	  
		    // is there a right subtree? if yes, delete it
	  
		    deleteSubtree(node->right);
		}

	    // delete the current node
      
	    Free(node);
	};
  
    // get back

    return 0;
};

// ------------------------------------------------

int FrequencyTree::computeFrequencies(Population *p)
{
    register long  i;
    register char *x;
    register long  N;

    // set the local variable for size of strings

    n = p->n;
    N = p->N;

    // delete the tree

    reset();

    // go through data and add each instance into the tree

    for (i=0; i<N; i++)
	{
	  x=p->individual[i].chromosome;
// 	    x = getChromosome(&(p->individual[i]),d);
	  addInstance(x);
	};

    // divide the leaves by the number of points (to get the frequencies)

    divideNodes(double(N));

    // get back

    return 0;
};


// -------------------------------------------------

TreeNode *FrequencyTree::addInstance(char *x)
{
    register int       i;
    register TreeNode *currentNode;

    // set some variables
  
    currentNode = root;

    // trace the x using index

    for (i=0; i<n; i++)
	{
	    // increase the number of instances corresponding to this node

	    currentNode->value++;

	    // is the bit a set?

	    if (x[i])
		{
		    if (!currentNode->right)
			allocateNewNode(currentNode->right);
	      
		    currentNode = currentNode->right;
		}
	    else
		{
		    if (!currentNode->left)
			allocateNewNode(currentNode->left);
	  
		    currentNode = currentNode->left;
		}
	};
 
    // new instance?

    if (!currentNode->value)
	numLeaves++;

    // increase the number of instances for this node (leaf)

    currentNode->value++;

    // get back

    return currentNode;
};

// -------------------------------------------------

int FrequencyTree::computeIndexedFrequencies(Population *p, int *index, int indexLength)
{
    register long  i;
    register char *x;
    register long  N;

    // set the local variable for size of strings

    n = indexLength;
    N = p->N;
  
    // delete the tree
  
    reset();

    // go through data and add each instance into the tree
  
    for (i=0; i<N; i++)
	{
	  x = p->individual[i].chromosome;
// 	    x = getChromosome(&(p->individual[i]),d);
	  addIndexedInstance(x,index);
	};

    // divide the leaves by the number of points (to get the frequencies)

    divideNodes(double(N));

    // get back
 
    return 0;
};
// -------------------------------------------------

TreeNode *FrequencyTree::addIndexedInstance(char *x, int *index)
{
    register int       i;
    register TreeNode *currentNode;
  
    // set some variables
  
    currentNode = root;
  
    // trace the x using index
  
    for (i=0; i<n; i++, index++)
	{
	    // increase the number of instances corresponding to this node

	    currentNode->value++;

	    // is the bit a set?

	    if (x[*index])
		{
		    if (!currentNode->right)
			allocateNewNode(currentNode->right);

		    currentNode = currentNode->right;
		}
	    else
		{
		    if (!currentNode->left)
			allocateNewNode(currentNode->left);
	  
		    currentNode = currentNode->left;
		}
	};

    // new instance?

    if (!currentNode->value)
	numLeaves++;

    // increase the number of instances for this node (leaf)

    currentNode->value++;

    // get back

    return currentNode;
};

// -------------------------------------------------

int FrequencyTree::computeIndexedFrequenciesWithSplit(Population *p, int *index, int indexLength, int **list, int *listLength, int listSize)
{
    register long  i;
    register char *x;
    register long  N;

    // delete the tree
  
    reset();

    // set the local variable for size of strings

    this->listSize   = listSize;
    n                = indexLength;
    N                = p->N;
  
    // go through data and add each instance into the tree
  
    for (i=0; i<N; i++)
	{
	  x=p->individual[i].chromosome;
// 	    x = getChromosome(&(p->individual[i]),d);
	    addIndexedInstanceWithSplit(x,index,list,listLength);
	};
  
    // divide the leaves by the number of points (to get the frequencies)

    divideNodes(double(N));

    // get back
 
    return 0;
};

// -------------------------------------------------

TreeNode *FrequencyTree::addIndexedInstanceWithSplit(char *x, int *index, int **list, int *listLength)
{
    register int       i;
    register TreeNode *currentNode;

    // set some variables
  
    currentNode = root;

    // trace the x using index

    for (i=0; i<n; i++, index++)
	{
	    // increase the number of instances corresponding to this node

	    currentNode->value++;

	    // is the bit a set?

	    if (x[*index])
		if (currentNode->right)
		    currentNode = currentNode->right;
		else
		    {
			if (i<n-1)
			    allocateNewNode(currentNode->right)
				else
				    allocateNewSplit(currentNode->right);

			currentNode = currentNode->right;
		    }
	    else
		if (currentNode->left)
		    currentNode = currentNode->left;
		else
		    {
			if (i<n-1)
			    allocateNewNode(currentNode->left)
				else
				    allocateNewSplit(currentNode->left);

			currentNode = currentNode->left;
		    };
	};

    // split empty?

    if (((TreeSplit*)currentNode)->subtree==NULL)
	{
	    ((TreeSplit*)currentNode)->numSubtrees = listSize;
	    ((TreeSplit*)currentNode)->subtree     = (FrequencyTree**) Calloc(listSize,sizeof(FrequencyTree*));
	    for (i=0; i<listSize; i++)
		{
		    ((TreeSplit*)currentNode)->subtree[i] = new FrequencyTree();
		    ((TreeSplit*)currentNode)->subtree[i]->n = listLength[i];
		};
	};
 
    // new instance?

    if (!currentNode->value)
	numLeaves++;

    // increase the number of instances for this node (leaf)

    currentNode->value++;

    // add the rest of the instance where it belongs

    for (i=0; i<listSize; i++)
	(((TreeSplit*)currentNode)->subtree[i])->addIndexedInstance(x,list[i]);
  
    // get back

    return currentNode;
};

// -------------------------------------------------

int FrequencyTree::divideNodes(double x)
{
    // divide the leaves in a tree starting at the root

    globalDivisor = x;

    divideSubtreeNodes(root);

    // get back

    return 0;
};

// -------------------------------------------------

int FrequencyTree::divideSubtreeNodes(TreeNode *node)
{
    // is the tree non-empty?

    if (node)
	{
	    // divide the node by the divisor (set globally)

	    node->value /= globalDivisor;

	    if (isLeaf(node))
		{
		    if (listSize)
			{
			    int i;
			    for (i=0; i<listSize; i++)
				((TreeSplit*)node)->subtree[i]->divideNodes(globalDivisor);
			};
		}
	    else
		{
		    // divide the subtrees
	  
		    divideSubtreeNodes(node->left);
		    divideSubtreeNodes(node->right);
		};
	};

    // get back

    return 0;
};

// -------------------------------------------------

long FrequencyTree::getNumInstances()
{
    // return the number of leaves

    return numLeaves;
};

// -------------------------------------------------

long FrequencyTree::getNumInstancesWithSplit(int l)
{
    // is the tree non-empty?

    return getSubtreeNumInstancesWithSplit(root,l);
};

// -------------------------------------------------

long FrequencyTree::getSubtreeNumInstancesWithSplit(TreeNode *node, int l)
{
    int result;

    result = 0;

    // check it out

    if (node->which==SPLIT)
	result += ((TreeSplit*) node)->subtree[l]->getNumInstances();
    else
	if (node->left)
	    {
		if (node->right)
		    result+=getSubtreeNumInstancesWithSplit(node->right,l);

		result+=getSubtreeNumInstancesWithSplit(node->left,l);
	    }
	else
	    {
		if (node->right)
		    result+=getSubtreeNumInstancesWithSplit(node->right,l);
	    };

    // get back with result
  
    return result;
};

// -------------------------------------------------

long FrequencyTree::getInstances(char **x)
{
    // get the instances from the tree beginning in the root

    if (root)
	return getSubtreeInstances(root,x,0,0);
    else
	return 0;
};

// -------------------------------------------------

long FrequencyTree::getSubtreeInstances(TreeNode *node, char **x, int depth, long i)
{
    long      newInstances; 

    // just do it

    if (node->left)
	if (node->right)
	    {
		newInstances  = (globalK = getSubtreeInstances(node->left,x,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth]=0;

		i += newInstances;

		newInstances += (globalK = getSubtreeInstances(node->right,x,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth]=1;
	    }
	else
	    {
		newInstances  = (globalK = getSubtreeInstances(node->left,x,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth]=0;
	    }
    else
	if (node->right)
	    {
		newInstances  = (globalK = getSubtreeInstances(node->right,x,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth]=1;
	    }
	else
	    newInstances    = 1;

    return newInstances;
};

// -------------------------------------------------

double FrequencyTree::getFrequency(char *x, int l)
{
    TreeNode *node, *next;

    node = root;

    for (int i=0; i<l; i++)
	{
	    if (x[i]==0)
		next = node->left;
	    else
		next = node->right;

	    if (next==NULL)
		return 0;
	    else
		node = next;
	};

    return node->value;
};

// -------------------------------------------------

long FrequencyTree::getFrequencyList(double *x)
{
    // get the instances from the tree beginning in the root

    if (root)
	return getSubtreeFrequencyList(root,x,0,0);
    else
	return 0;
};

// -------------------------------------------------

long FrequencyTree::getSubtreeFrequencyList(TreeNode *node, double *x, int depth, long i)
{
    long      newInstances; 

    if (node->left)
	if (node->right)
	    {
		newInstances  = getSubtreeFrequencyList(node->left,x,depth+1,i);
		i            += newInstances;
		newInstances += getSubtreeFrequencyList(node->right,x,depth+1,i);
	    }
	else
	    newInstances    = getSubtreeFrequencyList(node->left,x,depth+1,i);
    else
	if (node->right)
	    newInstances    = getSubtreeFrequencyList(node->right,x,depth+1,i);
	else
	    {
		x[i]          = node->value;
		newInstances  = 1;
	    };
  
    return newInstances;
};

// ------------------------------------------------

long FrequencyTree::getInstancesAndFrequencies(char **x, double *f)
{
    // get the instances from the tree beginning in the root

    if (root)
	return getSubtreeInstancesAndFrequencies(root,x,f,0,0);
    else
	return 0;  
};

// ------------------------------------------------

long FrequencyTree::getSubtreeInstancesAndFrequencies(TreeNode *node, char **x, double *f, int depth, long i)
{
    long      newInstances; 

    // just do it

    if (node->left)
	if (node->right)
	    {
		x[i][depth]   = 0;
		newInstances  = (globalK = getSubtreeInstancesAndFrequencies(node->left,x,f,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth] = 0;

		i            += newInstances;
	
		newInstances += (globalK = getSubtreeInstancesAndFrequencies(node->right,x,f,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth] = 1;
	    }
	else
	    {
		x[i][depth]   = 0;
		newInstances  = (globalK = getSubtreeInstancesAndFrequencies(node->left,x,f,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth] = 0;
	    }
    else
	if (node->right)
	    {
		x[i][depth]   = 1;
		newInstances  = (globalK = getSubtreeInstancesAndFrequencies(node->right,x,f,depth+1,i));

		for (globalK--; globalK>=0; globalK--)
		    x[i+globalK][depth] = 1;
	    }
	else
	    {
		f[i]            = node->value;
		newInstances    = 1;
	    };
  
    return newInstances;
};

// ------------------------------------------------

long FrequencyTree::getShiftedInstancesAndFrequencies(char **x, double *f, int shift)
{
    // get the instances from the tree beginning in the root

    if (root)
	return getSubtreeInstancesAndFrequencies(root,x,f,shift,0);
    else
	return 0;  
};

// ------------------------------------------------

long FrequencyTree::getInstancesAndFrequenciesWithSplit(char **x, double *f, int l)
{
    // get the instances from the tree beginning in the root

    return getSubtreeInstancesAndFrequenciesWithSplit(root,x,f,l,0,0);
};

// ------------------------------------------------

long FrequencyTree::getSubtreeInstancesAndFrequenciesWithSplit(TreeNode *node, char **x, double *f, int l, int depth, long i)
{
    long      newInstances;

    // just do it

    if (depth<n)
	{
	    if (node->left)
		if (node->right)
		    {
			newInstances  = (globalK = getSubtreeInstancesAndFrequenciesWithSplit(node->left,x,f,l,depth+1,i));
	    
			for (globalK--; globalK>=0; globalK--)
			    x[i+globalK][depth] = 0;
	    
			i            += newInstances;
	    
			newInstances += (globalK = getSubtreeInstancesAndFrequenciesWithSplit(node->right,x,f,l,depth+1,i));
	    
			for (globalK--; globalK>=0; globalK--)
			    x[i+globalK][depth] = 1;
		    }
		else
		    {
			newInstances  = (globalK = getSubtreeInstancesAndFrequenciesWithSplit(node->left,x,f,l,depth+1,i));
	    
			for (globalK--; globalK>=0; globalK--)
			    x[i+globalK][depth] = 0;
		    }
	    else
		if (node->right)
		    {
			newInstances  = (globalK = getSubtreeInstancesAndFrequenciesWithSplit(node->right,x,f,l,depth+1,i));
	    
			for (globalK--; globalK>=0; globalK--)
			    x[i+globalK][depth] = 1;
		    }
		else
		    {
			fprintf(stderr,"ERROR: Something weird is going on in FrequencyTree::getSubtreeInstancesAndFrequenciesWithSplit!\n");
			exit(-1);
		    }
	}
    else
	{      
	    newInstances = (globalK = ((TreeSplit*)node)->subtree[l]->getShiftedInstancesAndFrequencies(x+i,f+i,depth));
	};

    return newInstances;
};

// ------------------------------------------------

int FrequencyTree::reduceByIndex(TreeNode *source, int *index, int indexLength)
{
    // the depth of the tree is equal to the length of index (if there is something)

    n = indexLength;

  // pack the subtree starting in the root

    reduceSubtreeByIndex(&root,source,index,indexLength,0);

    // get back

    return 0;
};

// ------------------------------------------------

int FrequencyTree::reduceSubtreeByIndex(TreeNode **currentNode, TreeNode *sourceTree, int *index, int indexLength, int depth)
{
    // update for lists!!!!!!!!!!!!!!!!!

    if (sourceTree)
	{
	    if (*index==depth)
		if (indexLength==1)
		    {
			(*currentNode)->value = sourceTree->value;

			if (sourceTree->left)
			    {
				allocateNewNode((*currentNode)->left);
				(*currentNode)->left->value = sourceTree->left->value;

				numLeaves++;
			    }
			else
			    (*currentNode)->left = NULL;
	   
			if (sourceTree->right)
			    {
				allocateNewNode((*currentNode)->right);
				(*currentNode)->right->value = sourceTree->right->value;

				numLeaves++;
			    }
			else
			    (*currentNode)->right = NULL;
		    }
		else
		    {
			index++;
			indexLength--;
			depth++;

			(*currentNode)->value = sourceTree->value;

			allocateNewNode((*currentNode)->left);
			allocateNewNode((*currentNode)->right);

			reduceSubtreeByIndex(&((*currentNode)->left),sourceTree->left,index,indexLength,depth);
			reduceSubtreeByIndex(&((*currentNode)->right),sourceTree->right,index,indexLength,depth);
		    }
	    else
		{
		    TreeNode *newNode;

		    depth++;

		    allocateNewNode(newNode);

		    reduceSubtreeByIndex(currentNode,sourceTree->left,index,indexLength,depth);
		    reduceSubtreeByIndex(&newNode,sourceTree->right,index,indexLength,depth);

		    addTree(currentNode,newNode);

		    deleteSubtree(newNode);
		}
	}
    else
	{
	    if (*currentNode)
		Free(*currentNode);
	    *currentNode = NULL;
	};

    return 0;
};

// ------------------------------------------------

int FrequencyTree::addTree(TreeNode **destination, TreeNode *source)
{
    // update for lists !!!!!!!!!!!!!!!!!!

    if (source)
	{
	    if (*destination==NULL)
		*destination = copySubtree(source);
	    else
		{
		    (*destination)->value += source->value;

		    if (isLeaf(source))
			numLeaves--;
		    else
			{
			    addTree(&((*destination)->left),source->left);
			    addTree(&((*destination)->right),source->right);
			};
		};
	}

    return 0;
}

// ------------------------------------------------

long FrequencyTree::recomputeNumLeaves(TreeNode *node)
{
    if (node)
	if (isLeaf(node))
	    return 1;
	else
	    return recomputeNumLeaves(node->left) +
		recomputeNumLeaves(node->right);
    else
	return 0;
};
