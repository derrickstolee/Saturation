/*
 * SaturationAugmenter.hpp
 *
 *  Created on: Apr 17, 2011
 *      Author: stolee
 */

#ifndef SATURATIONAUGMENTER_HPP_
#define SATURATIONAUGMENTER_HPP_

#include "SaturationData.hpp"
#include "SearchManager.hpp"

#define _MODE_CANONICAL_ 1
#define _MODE_ORBITAL_ 2

class SaturationNode: public SearchNode
{
public:
	SaturationNode(LONG_T label);
	virtual ~SaturationNode();

	bool initialized;
	int numOrbits;
	int curOrbit;
	int* orbitReps;
	int curNumCompleters;
	int* curOrbitCompleterReps;
	int chosenOrbit;
};

/**
 * SaturationAugmenter performs augmentations
 * for the uniquely K_r-saturated graph search.
 *
 * There are two modes:
 * 	1: Canonical -- Use canonical deletion to be
 * 		isomorph-free.
 *  2: Orbital -- Use orbital branching with constraint propagation
 *  	to search for solutions (with possible repetition)
 */
class SaturationAugmenter: public SearchManager
{
protected:
	typedef SearchManager super;

	/**
	 * mode -- Which mode is the augmentation using?
	 *
	 * Either _MODE_CANONICAL_ or
	 *  _MODE_ORBITAL_
	 */
	int mode;

	int r;
	int maxN;

	/**
	 * All data is stored in a SaturationData object
	 */
	SaturationData* satdata;

	/**
	 * augmentation_succeeded
	 */
	bool augmentation_succeeded;

	/**
	 * initializeNode
	 *
	 * Run the symmetry methods on this node
	 */
	void initializeNode(SaturationNode* node, int& orbit, int& completer);

	/**
	 * getInitialOrbitCompleter
	 *
	 * Pick the orbit and completer to START the augmentations.
	 */
	void getInitialOrbitCompleter(SaturationNode* node, int& orbit, int& completer);

	/**
	 * getNextOrbitCompleter
	 *
	 * Given the current node, orbit, and completer,
	 * 	use the mode to decide the next orbit and completer.
	 */
	void getNextOrbitCompleter(SaturationNode* node, int& orbit, int& completer);

public:
	/**
	 * Constructor
	 */
	SaturationAugmenter(int r, int maxN, int mode);

	/**
	 * Destructor
	 */
	virtual ~SaturationAugmenter();

	/**
	 * pushNext -- deepen the search to the next child
	 * 	of the current node.
	 *
	 * @return the label for the new node. -1 if none.
	 */
	virtual LONG_T pushNext();

	/**
	 * prune -- Perform a check to see if this nod
	 * pushTo -- deepen the search to the specified child
	 * 	of the current node.
	 *
	 * @param child the specified label for the new node
	 * @return the label for the new node. -1 if none, or failed.
	 */
	virtual LONG_T pushTo(LONG_T child);

	/**
	 * pop -- remove the current node and move up the tree.
	 *
	 * @return the label of the node after the pop.
	 * 		This return value is used for validation purposes
	 * 		to check proper implementation of push*() and pop().
	 */
	virtual LONG_T pop();

	/**
	 * prune -- Perform a check to see if this node should be pruned.
	 *
	 * @return 0 if no prune should occur, 1 if prune should occur.
	 */
	virtual int prune();

	/**
	 * isSolution -- Perform a check to see if a solution exists
	 * 		at this point.
	 *
	 * @return 0 if no solution is found, 1 if a solution is found.
	 */
	virtual int isSolution();

	/**
	 * writeSolution -- create a buffer that contains a
	 * 	description of the solution.
	 *
	 * Solution strings start with S.
	 *
	 * @return a string of data allocated with malloc().
	 * 	It will be deleted using free().
	 */
	virtual char* writeSolution();

	/**
	 * writeStatistics -- create a buffer that contains a
	 * 	description of the solution.
	 *
	 * Statistics take the following format in each line:
	 *
	 * T [TYPE] [ID] [VALUE]
	 *
	 * @return a string of data allocated with malloc().
	 * 	It will be deleted using free().
	 */
	virtual char* writeStatistics();

};

#endif /* SATURATIONAUGMENTER_HPP_ */
