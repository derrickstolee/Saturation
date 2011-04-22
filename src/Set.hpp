/*
 * Set.hpp
 *
 *  Created on: Feb 3, 2010
 *      Author: derrickstolee
 *
 * Oct, 13, 2010: Refactored to make Set an interface.
 */

#ifndef SET_HPP_
#define SET_HPP_


/**
 * A set stores a set of numerical values in a binary tree.
 *
 * The specific storage mechanism is a B-tree.
 */
class Set
{

public:

	virtual ~Set();

	/**
	 * size - return the number of elements in the set.
	 *
	 * @return the number of elements
	 */
	virtual int size() = 0;

	/**
	 * contains - checks if an element is in the set
	 *
	 * @param elt the element tho find.
	 * @return 1 if elt is in the Set. 0 otherwise.
	 */
	virtual int contains(int elt) = 0;

	/**
	 * remove - removes an element from the set
	 *
	 * @param elt the element to find
	 * @return elt if the element was found and removed. -1 otherwise.
	 */
	virtual int remove(int elt) = 0;

	/**
	 * add - adds an element to the set
	 *
	 * @param elt the element to add
	 * @return elt
	 */
	virtual int add(int elt) = 0;

	/**
	 * join - form a new Set as the union of this and another Set.
	 *
	 * @param set the other set
	 * @return a new Set filled with elements from both sets
	 */
	virtual Set* join(Set* set) = 0;

	/**
	 * intersect - form a new Set as the intersection of this and another Set.
	 *
	 * @param set the other set
	 * @return a new Set filled with elements in both sets
	 */
	virtual Set* intersect(Set* set) = 0;

	/**
	 * resetIterator - reset the iteration through the set.
	 * That is, we want to enumerate the elements in the set, but
	 * need to start from the beginning after this call.
	 */
	virtual void resetIterator() = 0;

	/**
	 * next - get the next element in the Set.
	 *
	 * @return the next element in the traversal
	 */
	virtual int next() = 0;

	/**
	 * hasNext - there is another element in the Set that we have not
	 * enumerated yet.
	 *
	 * @return 1 if there is, 0 if there is not.
	 */
	virtual int hasNext() = 0;

	/**
	 * clear -- remove all elements from the set.
	 */
	virtual void clear() = 0;

	/**
	 * min -- get the minimum key.
	 */
	virtual int min() = 0;

	/**
	 * max -- get the maximum key.
	 */
	virtual int max() = 0;

	/**
	 * toArray -- get the set as an array that the set will clean up.
	 */
	virtual int* toArray() = 0;

	/**
	 * print
	 */
	virtual void print() = 0;

};

#endif /* SET_HPP_ */
