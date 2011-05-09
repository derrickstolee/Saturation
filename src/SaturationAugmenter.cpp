/*
 * SaturationAugmenter.cpp
 *
 *  Created on: Apr 17, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "translation.hpp"
#include "SaturationData.hpp"
#include "SaturationSymmetry.hpp"
#include "SaturationAugmenter.hpp"

/**
 * Used for nauty statistics
 */
extern double nauty_time;
extern double nauty_layer_time;
extern double nauty_zero_time;
extern double nauty_stab_time;
extern double orbit_layer_time;
extern double orbit_zero_time;
extern double orbit_stab_time;
extern long long int nauty_calls;
extern long long int nauty_layer_calls;
extern long long int nauty_stab_calls;
extern long long int nauty_zero_calls;

SaturationNode::SaturationNode( LONG_T label ) :
	SearchNode(label)
{
	this->curChild = -1;
	this->initialized = false;
	this->numOrbits = -1;
	this->curOrbit = -1;
	this->orbitReps = 0;
	this->curNumCompleters = -1;
	this->curOrbitCompleterReps = 0;
	this->chosenOrbit = -1;
}

SaturationNode::~SaturationNode()
{
	if ( this->orbitReps != 0 )
	{
		free(this->orbitReps);
		this->orbitReps = 0;
	}

	if ( this->curOrbitCompleterReps != 0 )
	{
		free(this->curOrbitCompleterReps);
		this->curOrbitCompleterReps = 0;
	}
}

/**
 * initializeNode
 *
 * Run the symmetry methods on this node
 */
void SaturationAugmenter::initializeNode( SaturationNode* node, int& orbit, int& completer )
{
	if ( node->numOrbits < 0 )
	{

#ifndef _AUGMENT_ALL_VERTICES_FIRST_
		if ( mode == _MODE_ORBITAL_ && this->satdata->getNum2Edges() == 0 )
		{
			/* add a vertex in this case */
			bool success = this->satdata->augmentVertex();
			if ( !success )
			{
				return;
			}
		}
#else
		if ( mode == _MODE_ORBITAL_ && this->satdata->getNum2Edges() == 0 )
		{
			node->numOrbits = 0;
			return;
		}
#endif

		node->orbitReps = this->satdata->get2Orbits(node->numOrbits);
	}

	if ( node->numOrbits > 0 )
	{
		if ( orbit < 0 )
		{
			orbit = 0;
		}

		if ( node->curOrbit != orbit )
		{
			node->curOrbit = orbit;

			if ( node->curOrbitCompleterReps != 0 )
			{
				free(node->curOrbitCompleterReps);
				node->curOrbitCompleterReps = 0;
			}

			if ( node->curOrbit < node->numOrbits )
			{
				node->curOrbitCompleterReps = this->satdata->getCompletionOrbits(node->orbitReps[node->curOrbit],
				                                                                 node->curNumCompleters);
			}
		}
	}

	node->initialized = true;
}

/**
 * getInitialOrbitCompleter
 *
 * Pick the orbit and completer to START the augmentations.
 */
void SaturationAugmenter::getInitialOrbitCompleter( SaturationNode* node, int& orbit, int& completer )
{
	/* we need to compute orbits for this data */
	this->initializeNode(node, orbit, completer);

	if ( this->mode == _MODE_CANONICAL_ )
	{
		orbit = 0;
		completer = 0;

		if ( node->numOrbits <= 0 )
		{
			/* there are no orbits? */
			completer = 0xFFFF;
		}
		else
		{
			/* we need to compute completion orbits for fixing this orbit */
			this->initializeNode(node, orbit, completer);

			while ( orbit < node->numOrbits && node->curNumCompleters == 0 )
			{
				orbit++;

				if ( orbit < node->numOrbits )
				{
					this->initializeNode(node, orbit, completer);
				}
			}
		}
	}
	else if ( this->mode == _MODE_ORBITAL_ )
	{
		/* We need to select an orbit by largest size */
		orbit = 0;
		completer = 0;

		if ( node->numOrbits > 1 )
		{
			int maxOrbitSize = 0;
			double min0degree = 0;
			double min1degree = this->maxN * this->maxN * 2;
			int nChoose2 = nChooseK(this->satdata->getN(), 2);
			int* orbit_sizes = (int*) malloc(nChoose2 * sizeof(int));
			bzero(orbit_sizes, nChoose2 * sizeof(int));

			int orbitRep = 0;
			for ( int j = 0; j < nChoose2; j++ )
			{
				if ( this->satdata->getAdjacency(j) == 2 )
				{
					int orbitIndex = this->satdata->get2OrbitIndexByRepresentative(j);
					int osize = orbit_sizes[orbitIndex] + 1;
					orbit_sizes[orbitIndex] = osize;
				}
			}

			/* Select the "best" orbit */
			/* TODO: Insert branching rules here! */
			for ( int j = 0; j < node->numOrbits; j++ )
			{
				int ej = node->orbitReps[j];
				if ( this->satdata->getAdjacency(ej) == 2 )
				{
					int osize = orbit_sizes[j];
					int vi, vj;
					indexToPair(ej, vi, vj);
					double deg0sum = osize / (double) (1.0 + this->satdata->getDegree(vi, 0)
					        + this->satdata->getDegree(vj, 0));
					double deg1sum = osize / (double) (1.0 + this->satdata->getDegree(vi, 1)
					        + this->satdata->getDegree(vj, 1));

#ifdef _DEG_SUM_ORBIT_CHOICE
					if ( deg0sum > min0degree )
					{
						orbit = j;
						maxOrbitSize = osize;
						min0degree = deg0sum;
						min1degree = deg1sum;
					}
					else if ( deg0sum == min0degree && deg1sum < min1degree )
					{
						orbit = j;
						maxOrbitSize = osize;
						min1degree = deg1sum;
					}
					else if ( deg0sum == min0degree && deg1sum == min1degree && osize > maxOrbitSize )
					{
						orbit = j;
						maxOrbitSize = osize;
					}
#else
					if ( osize > maxOrbitSize )
					{
						orbit = j;
						maxOrbitSize = osize;
					}
#endif
				}
			}
			free(orbit_sizes);
			node->chosenOrbit = orbit;
			//			printf("Selected orbit %d of size %d (numOrbits %d).\n", orbit, maxOrbitSize, node->numOrbits);
		}
		else if ( node->numOrbits == 1 )
		{
			//			printf("Only one orbit.\n");
			node->chosenOrbit = 0;
		}

		if ( node->numOrbits <= 0 )
		{
			/* we can't do anything here */
			node->chosenOrbit = -1;
			orbit = node->numOrbits + 3;
			completer = 0xFFFF;
		}
		else
		{
			/* we need to compute completion orbits for fixing this orbit */
			this->initializeNode(node, orbit, completer);

			if ( node->curNumCompleters <= 0 )
			{
				/* just jump to filling */
				orbit = node->numOrbits;
				completer = 0xFFFF;
			}
		}
	}
}

/**
 * getNextOrbitCompleter
 *
 * Given the current node, orbit, and completer,
 * 	use the mode to decide the next orbit and completer.
 */
void SaturationAugmenter::getNextOrbitCompleter( SaturationNode* node, int& orbit, int& completer )
{
	if ( orbit < node->numOrbits && completer < node->curNumCompleters - 1 )
	{
		/* just use the next completer */
		completer++;
		return;
	}

	if ( node->numOrbits <= 0 && this->mode == _MODE_ORBITAL_ )
	{
		node->chosenOrbit = -1;
		/* just jump to adding a vertex */
		orbit = node->numOrbits + 3;
		completer = 0xFFFF;
	}

	if ( orbit >= node->numOrbits || this->mode == _MODE_CANONICAL_ )
	{
		/* try all orbits */
		orbit++;
	}
	else
	{
		/* ORBITAL and below numOrbits */
		if ( node->numOrbits == 0 )
		{
			orbit = node->numOrbits + 1;
		}
		else
		{
			/* just jump to filling */
			orbit = node->numOrbits;
		}
	}

	if ( orbit >= node->numOrbits + 1 && mode == _MODE_ORBITAL_ )
	{
		orbit = -1;
		completer = -1;
		return;
	}

	if ( orbit >= node->numOrbits )
	{
		completer = 0xFFFF;
	}
	else
	{

		completer = 0;

		/* COMPUTE NEW curNumCompleters */
		this->initializeNode(node, orbit, completer);

		while ( orbit < node->numOrbits && node->curNumCompleters == 0 )
		{
			orbit++;

			if ( orbit < node->numOrbits )
			{
				this->initializeNode(node, orbit, completer);
			}
		}
	}

	if ( orbit >= node->numOrbits + 2 )
	{
		/* we have surpassed the possible augmentations */
		orbit = -1;
		completer = -1;
	}
}

/**
 * Constructor
 */
SaturationAugmenter::SaturationAugmenter( int r, int maxN, int mode )
{
	this->haltAtSolutions = false;

	this->r = r;
	this->maxN = maxN;
	this->mode = mode;

	this->satdata = new SaturationData(r, maxN);

#ifdef _AUGMENT_ALL_VERTICES_FIRST_
	while ( this->satdata->getN() < maxN )
	{
		this->satdata->augmentVertex();
	}
#endif

	if ( this->root != 0 )
	{
		delete this->root;
		this->root = new SaturationNode(0);

		this->root->curChild = -1;
	}
}

/**
 * Destructor
 */
SaturationAugmenter::~SaturationAugmenter()
{
	delete this->satdata;
}

/**
 * pushNext -- deepen the search to the next child
 * 	of the current node.
 *
 * @return the label for the new node. -1 if none.
 */
LONG_T SaturationAugmenter::pushNext()
{
	SaturationNode* node = 0;
	if ( this->stack.size() == 0 )
	{
		/* we need to look at the root */
		node = (SaturationNode*) this->root;
	}
	else
	{
		node = (SaturationNode*) this->stack.back();
	}

	LONG_T child = node->curChild;

	int orbit = child >> 16;
	int completer = child & 0xFFFF;

	if ( child < 0 )
	{
		this->getInitialOrbitCompleter(node, orbit, completer);
	}
	else
	{
		this->getNextOrbitCompleter(node, orbit, completer);
	}

#ifdef _AUGMENT_ALL_VERTICES_FIRST_
	if ( node->numOrbits <= 0 )
	{
		return -1;
	}
#endif

	if ( orbit < 0 )
	{
		return -1;
	}

	child = ((LONG_T) orbit << 16) | (LONG_T) (completer & 0xFFFF);

	return this->pushTo(child);
}

/**
 * pushTo -- deepen the search to the specified child
 * 	of the current node.
 *
 * @param child the specified label for the new node
 * @return the label for the new node. -1 if none, or failed.
 */
LONG_T SaturationAugmenter::pushTo( LONG_T child )
{
	SaturationNode* node = 0;

	if ( this->stack.size() == 0 )
	{
		node = (SaturationNode*) this->root;
	}
	else
	{
		node = (SaturationNode*) this->stack.back();
	}

	node->curChild = child;

	int orbit = child >> 16;
	int completer = child & 0xFFFF;

	/* if Node is not initialized with orbit information, do so! */
	if ( node->initialized == false || node->curOrbit != orbit )
	{
		this->initializeNode(node, orbit, completer);
	}

#ifdef _AUGMENT_ALL_VERTICES_FIRST_
	if ( node->numOrbits <= 0 )
	{
		return -1;
	}
#endif

	/********************************** TAKE SNAPSHOT **************/
	this->satdata->snapshot();

	if ( orbit < node->numOrbits )
	{
		/* Augment with a completer */
		int i, j;
		indexToPair(node->orbitReps[orbit], i, j);
		this->augmentation_succeeded = this->satdata->augmentCompletion(i, j, node->curOrbitCompleterReps[completer]);
	}
	else if ( orbit == node->numOrbits )
	{
		/* Fill */
		if ( this->mode == _MODE_CANONICAL_ )
		{
			this->augmentation_succeeded = this->satdata->augmentFill(-1);
		}
		else
		{
			this->augmentation_succeeded = this->satdata->augmentFill(node->chosenOrbit);
		}
	}
	else if ( orbit == node->numOrbits + 1 && mode == _MODE_CANONICAL_ )
	{
		/* Add a vertex */
		this->augmentation_succeeded = this->satdata->augmentVertex();
	}
	else
	{
		return -1;
	}

	this->stack.push_back(new SaturationNode(child));

	return child;
}

/**
 * pop -- remove the current node and move up the tree.
 *
 * @return the label of the node after the pop.
 * 		This return value is used for validation purposes
 * 		to check proper implementation of push*() and pop().
 */
LONG_T SaturationAugmenter::pop()
{
	if ( this->stack.size() > 0 )
	{
		SearchNode* node = this->stack.back();
		this->stack.pop_back();

		LONG_T label = node->label;
		delete node;

		this->satdata->rollback();
		return label;
	}

	return this->root->label;
}

/**
 * prune -- Perform a check to see if this node should be pruned.
 *
 * @return 0 if no prune should occur, 1 if prune should occur.
 */
int SaturationAugmenter::prune()
{
	if ( this->augmentation_succeeded == false )
	{
		this->augmentation_succeeded = true;
		return 1;
	}

	if ( this->mode == _MODE_CANONICAL_ )
	{
		//		this->satdata->printEdgeStack();
		//		printf("%s", this->satdata->getOneGraphString());

		if ( this->satdata->isCanonical() == false )
		{
			//			printf("-- NOT CANONICAL\n\n");
			return 1;
		}
	}

	if ( this->mode == _MODE_ORBITAL_ )
	{
		/* for all 2-edges with completions, complete them! */
		if ( this->satdata->has2TypeWithCompletion() )
		{
			bool success = this->satdata->augmentAutoCompletions();

			if ( !success )
			{
				/* Did this ever fail? */
				return 1;
			}
		}
	}

	return 0;
}

/**
 * isSolution -- Perform a check to see if a solution exists
 * 		at this point.
 *
 * @return 0 if no solution is found, 1 if a solution is found.
 */
int SaturationAugmenter::isSolution()
{
	if ( this->augmentation_succeeded )
	{
		if ( this->satdata->isUniquelyKrSaturated() )
		{
			return 1;
		}
	}

	return 0;
}

/**
 * writeSolution -- create a buffer that contains a
 * 	description of the solution.
 *
 * Solution strings start with S.
 *
 * @return a string of data allocated with malloc().
 * 	It will be deleted using free().
 */
char* SaturationAugmenter::writeSolution()
{
	//	this->satdata->printEdgeStack();
	char* onestr = this->satdata->getOneGraphString();
	char* os_copy = (char*) malloc(strlen(onestr) + 5);
	strcpy(os_copy, onestr);
	return os_copy;
}

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
char* SaturationAugmenter::writeStatistics()
{
	char* buffer = super::writeStatistics();

	/* THIS IS TOTALLY CHEATING */
	printf("T SUM NUM_NAUTY_CALLS %lld\n", nauty_calls);
	printf("T SUM NUM_LAYER_CALLS %lld\n", nauty_layer_calls);
	printf("T SUM NUM_ZERO_CALLS %lld\n", nauty_zero_calls);
	printf("T SUM NUM_STAB_CALLS %lld\n", nauty_stab_calls);
	printf("T SUM NAUTY_TIME %lf\n", nauty_time);
	printf("T SUM NAUTY_LAYER_TIME %lf\n", nauty_layer_time);
	printf("T SUM NAUTY_ZERO_TIME %lf\n", nauty_zero_time);
	printf("T SUM NAUTY_STAB_TIME %lf\n", nauty_stab_time);
	printf("T SUM ORBIT_LAYER_TIME %lf\n", orbit_layer_time);
	printf("T SUM ORBIT_ZERO_TIME %lf\n", orbit_zero_time);
	printf("T SUM ORBIT_STAB_TIME %lf\n", orbit_stab_time);

	return buffer;
}
