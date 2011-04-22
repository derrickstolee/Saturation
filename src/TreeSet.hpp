/*
 * TreeSet.hpp
 *
 *  Created on: Feb 3, 2010
 *      Author: derrickstolee
 *
 * Oct 14, 2010: Refactored Set.hpp into TreeSet.hpp.
 */

#ifndef TREESET_HPP_
#define TREESET_HPP_

#include "Set.hpp"

#define BTREE_T 3
#define NUM_B_KEYS 10
#define NUM_B_CHILDREN 10

typedef struct bnode_struct* bnode_ptr;
typedef struct bnode_struct
{
	/* number of keys in this node */
	int key_count;

	/* whether a leaf or not */
	int leaf;

	/* the keys */
	int* keys;

	/* the child nodes */
	bnode_ptr* children;

} bnode;

/**
 * binit - initialize the root of a tree
 */
bnode_ptr binit();

/**
 * bsearch - recursively search the tree for an element.
 *
 * @param elt the element
 * @return elt if it is there, -1 otherwise
 */
int bsearch(bnode_ptr node, int elt);

/**
 * binsert - insert an element into the btree.
 *
 * @return the new root of the tree
 */
bnode_ptr binsert(bnode_ptr node, int elt);

/**
 * bsplit - split a node into parts at position 'pos'
 */
bnode_ptr bsplit(bnode_ptr node, int pos);

/**
 * bremove - remove an element 'elt' from the node
 *
 * @return the new root, if needed
 */
bnode_ptr bremove(bnode_ptr node, int elt);

/**
 * bfree - free up the memory from a node, including recursing below.
 */
void bfree(bnode_ptr node);

/**
 * bnext - the element in the set strictly larger than this element.
 *
 * @param node - the node to recurse over
 * @param elt - the element to be strictly larger
 * @return the next element in the list.  -1 if no more.
 */
int bnext(bnode_ptr node, int elt);

/**
 * bmax - get the maximum element
 */
int bmax(bnode_ptr node);

/**
 * bmin - get the minimum element
 */
int bmin(bnode_ptr node);

/**
 * bprint - output everything about this node and its children
 */
void bprint(bnode_ptr node);

/**
 * A set stores a set of numerical values in a binary tree.
 *
 * The specific storage mechanism is a B-tree.
 */
class TreeSet : public Set
{
protected:
	/* the mininum degree of the B-tree */
	int t;

	/* the root of the tree */
	bnode_ptr root;

	/* the number of elements in the tree */
	int count;

	/* the current iterator element */
	int cur_elt;

	/* the next iterator element */
	int next_elt;

	/* the array */
	int* array;

	/* the array is up-to-date? */
	bool array_up_to_date;

public:
	/**
	 * Constructor.
	 */
	TreeSet();

	/**
	 * Constructor.
	 *
	 * @param T - the minimum degree
	 */
	TreeSet(int T);

	/**
	 * Copy Constructor
	 */
	TreeSet(TreeSet& s);

	/**
	 * Destructor
	 */
	virtual ~TreeSet();

	/**
	 * size - return the number of elements in the set.
	 *
	 * @return the number of elements
	 */
	virtual int size();

	/**
	 * contains - checks if an element is in the set
	 *
	 * @param elt the element tho find.
	 * @return 1 if elt is in the Set. 0 otherwise.
	 */
	virtual int contains(int elt);

	/**
	 * remove - removes an element from the set
	 *
	 * @param elt the element to find
	 * @return elt if the element was found and removed. -1 otherwise.
	 */
	virtual int remove(int elt);

	/**
	 * add - adds an element to the set
	 *
	 * @param elt the element to add
	 * @return elt
	 */
	virtual int add(int elt);

	/**
	 * join - form a new Set as the union of this and another Set.
	 *
	 * @param set the other set
	 * @return a new Set filled with elements from both sets
	 */
	virtual Set* join(Set* set);

	/**
	 * intersect - form a new Set as the intersection of this and another Set.
	 *
	 * @param set the other set
	 * @return a new Set filled with elements in both sets
	 */
	virtual Set* intersect(Set* set);

	/**
	 * resetIterator - reset the iteration through the set.
	 * That is, we want to enumerate the elements in the set, but
	 * need to start from the beginning after this call.
	 */
	virtual void resetIterator();

	/**
	 * next - get the next element in the Set.
	 *
	 * @return the next element in the traversal
	 */
	virtual int next();

	/**
	 * hasNext - there is another element in the Set that we have not
	 * enumerated yet.
	 *
	 * @return 1 if there is, 0 if there is not.
	 */
	virtual int hasNext();

	/**
	 * clear -- remove all elements from the set.
	 */
	virtual void clear();

	/**
	 * min -- get the minimum key.
	 */
	virtual int min();

	/**
	 * max -- get the maximum key.
	 */
	virtual int max();

	/**
	 * toArray -- get the set as an array that the set will clean up.
	 */
	virtual int* toArray();

	/**
	 * print
	 */
	virtual void print();

};

#endif /* TREESET_HPP_ */
