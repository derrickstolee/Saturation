/*
 * SaturationManager.cpp
 *
 *  Created on: Apr 6, 2011
 *      Author: stolee
 */

#include "SaturationGraph.hpp"
#include "SaturationManager.hpp"
#include <stdlib.h>
#include <string.h>

void printString( SaturationGraph* satgraph )
{
	char* str = satgraph->getString();
	printf("%s", str);
	free(str);
}

/**
 * Constructor
 *
 * @param r the 'saturation parameter' which determines
 * 	which K_r is used for unique saturation.
 * @param N the maximum number of vertices to use.
 */
SaturationManager::SaturationManager( int r, int N )
{
	this->r = r;
	this->N = N;

	this->satgraph = new SaturationGraph(r, N);
<<<<<<< HEAD

=======
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
	this->haltAtSolutions = false;
}

/**
 * Destructor
 */
SaturationManager::~SaturationManager()
{
	if ( this->satgraph != 0 )
	{
		delete this->satgraph;
		this->satgraph = 0;
	}
}

/**
 * pushNext -- deepen the search to the next child
 * 	of the current node.
 *
 * @return the label for the new node. -1 if none.
 */
LONG_T SaturationManager::pushNext()
{
	if ( this->stack.size() == 0 )
	{
		/* this is the first augmentation */

		/* (r) types of first augmentations:
		 *
		 * Exactly two non-edge types:
		 * 	(a) Incident to the non-edge and the extra vertex
		 *  (b) Incident to the completion and the extra vertex
		 *
		 * (r-1) types of completers,
		 * 	use i vertices from the current completer,
		 * 	for each i in {r-2,...,1,0}.
		 *
		 *   */
		LONG_T label = this->root->curChild;

		label++;
		int orbit = (int) (label >> 24);
		int completer = (int) (label & 0xFFF);

		if ( orbit == 0 && completer >= r - 1 )
		{
			/* Wrap-around */
			completer = 0;
			orbit++;
		}
		else if ( orbit >= 1 && orbit <= 2 && completer >= r - 2 )
		{
			completer = 0;
			orbit++;
		}

		if ( orbit >= 3 )
		{
			return -1;
		}

		label = (((LONG_T) orbit) << 24) | (LONG_T) completer;
		LONG_T result = this->pushTo(label);
		return result;
	}

	/* Get the label from the current top of the stack */
	SearchNode* node = this->stack.back();

	LONG_T label = node->curChild;

	int orbit = (int) (label >> 24);
	int completer = (int) (label & 0xFFF);

	if ( orbit < 0 )
	{
		orbit = 0;
		completer = 0;
	}
	else if ( completer < 0 )
	{
		completer = 0;
	}
	else
	{
		completer++;
	}

	/* modify to the next orbit/completer pair */
	int noni, nonj;

	/* this may trigger a nauty call */
	int num_orbits = this->satgraph->numOrbits();
	if ( orbit < num_orbits )
	{
		int aug_size = this->satgraph->getDepth();
		int my_size = this->stack.size();

		int* orbit_array = this->satgraph->getOrbit(orbit);
		int orbit_rep_pair = orbit_array[0];

		this->satgraph->indexToPair(orbit_rep_pair, noni, nonj);

		if ( completer >= this->satgraph->numStabilizedOrbits(noni, nonj) )
		{
			orbit++;
			completer = 0;
		}
	}
	else
	{
		orbit++;
		completer = 0;

		if ( orbit > this->satgraph->numOrbits() + 1 )
		{
			/* out of new orbit-completer pairs */
			return -1;
		}
	}

	label = (((LONG_T) orbit) << 24) | (LONG_T) completer;
	LONG_T result = this->pushTo(label);
	return result;
}

/**
 * pushTo -- deepen the search to the specified child
 * 	of the current node.
 *
 * The LABEL for this augmentation is:
 * (a) the orbit-index for the non-edge selection.
 * (b) the stabilized-orbit-index for the completion.
 *
 * @param child the specified label for the new node
 * @return the label for the new node. -1 if none, or failed.
 */
LONG_T SaturationManager::pushTo( LONG_T child )
{
	int orbitIndex = (int) (child >> 24);
	int completionIndex = (int) (child & 0xFFF);

	if ( this->stack.size() == 0 )
	{
		/* this is the first augmentation */

		/* (r) types of first augmentations:
		 *
		 * Exactly two non-edge types:
		 * 	(a) Incident to the non-edge and the extra vertex
		 *  (b) Incident to the completion and the extra vertex
		 *
		 * (r-1) types of completers,
		 * 	use i vertices from the current completer,
		 * 	for each i in {r-2,...,1,0}.
		 *
		 *   */
		int orbit = orbitIndex;
		int completer = completionIndex;

		if ( orbit == 0 )
		{
			/* first vertex is 0: on the non-edge */
			/* second vertex is r: the new vertex */
			/* completion is given by 'completer': *
			 * (completer current verts)+(r-2-completer other verts)
			 */
			int* completion = (int*) malloc((r - 2) * sizeof(int));

			/* there are r-2 possible completions */
			for ( int i = 0; i < completer; i++ )
			{
				completion[i] = 2 + i;
			}
			for ( int i = completer; i < r - 2; i++ )
			{
				/* the r-2-completer vertices beyond the first r vertices */
				completion[i] = r + (i - completer) + 1;
			}

			this->augmentation_suceeded = this->satgraph->augment(0, r, completion);

			if ( this->augmentation_suceeded == false )
			{
				/* completion was copied by the augment */
				free(completion);
			}
		}
		else if ( orbit == 1 )
		{
			/* first vertex is 2: off the non-edge */
			/* second vertex is r: the new vertex */
			/* completion is given by 'completer', which
			 * does NOT use 0 or 1: *
			 * 	{2,3,...,r-1-completer} + { r+1,r+2,...,r-2-completer }
			 */
			int* completion = (int*) malloc((r - 2) * sizeof(int));

			/* there are r-3 possible completions */
			for ( int i = 0; i < completer; i++ )
			{
				completion[i] = 3 + i;
			}
			for ( int i = completer; i < r - 2; i++ )
			{
				/* the r-2-completer vertices beyond the first r vertices */
				completion[i] = r + (i - completer) + 1;
			}

			this->augmentation_suceeded = this->satgraph->augment(2, r, completion);

			if ( this->augmentation_suceeded == false )
			{
				/* completion was copied by the augment */
				free(completion);
			}
		}
		else if ( orbit == 2 )
		{
			/* first vertex is 2: off the non-edge */
			/* second vertex is r: the new vertex */
			/* completion is given by 'completer', which
			 * does uses 0: *i
			 * 	{0,2,3,...,r-2-completer} + { r+1,r+2,...,r-2-completer }
			 */
			int* completion = (int*) malloc((r - 2) * sizeof(int));

			completion[0] = 0;
			for ( int i = 1; i < completer; i++ )
			{
				completion[i] = 2 + i; /* really: 3+(i-1) */
			}
			for ( int i = completer; i < r - 2; i++ )
			{
				/* the r-2-completer vertices beyond the first r vertices */
				completion[i] = r + (i - completer) + 1;
			}

			this->augmentation_suceeded = this->satgraph->augment(2, r, completion);

			if ( this->augmentation_suceeded == false )
			{
				/* completion was copied by the augment */
				free(completion);
			}
		}
		else
		{
			return -1;
		}

		LONG_T label = (((LONG_T) orbit) << 24) | (LONG_T) completer;

		this->stack.push_back(new SearchNode(label));
		this->root->curChild = label;

		return label;
	}

	SearchNode* node = this->stack.back();

	/* modify to the next orbit/completer pair */
	int noni, nonj;

	/* this may trigger a nauty call */
	if ( orbitIndex > this->satgraph->numOrbits() + 1 )
	{
		return -1;
	}
	else if ( orbitIndex == this->satgraph->numOrbits() )
	{
		/* Special Augmentation: New Vertex */
		this->augmentation_suceeded = this->satgraph->addVertex();

		LONG_T label = (((LONG_T) orbitIndex) << 24) | (LONG_T) completionIndex;
		node->curChild = label;

		this->stack.push_back(new SearchNode(label));

		return label;
	}
	else if ( orbitIndex == this->satgraph->numOrbits() + 1 )
	{
		/* Special Augmentation: Fill */
		this->augmentation_suceeded = this->satgraph->fill();

		LONG_T label = (((LONG_T) orbitIndex) << 24) | (LONG_T) completionIndex;
		node->curChild = label;

		this->stack.push_back(new SearchNode(label));

		return label;
	}
	else if ( orbitIndex < this->satgraph->numOrbits() )
	{
		int* orbit_array = this->satgraph->getOrbit(orbitIndex);
		int orbit_rep_pair = orbit_array[0];

		this->satgraph->indexToPair(orbit_rep_pair, noni, nonj);

		if ( completionIndex >= this->satgraph->numStabilizedOrbits(noni, nonj) )
		{
			/* the completionIndex is too high */
			return -1;
		}

		/* get the completion representative */
		int* completion = this->satgraph->getStabilizedCompletion(noni, nonj, completionIndex);

		int* completecopy = (int*) malloc((this->r - 2) * sizeof(int));
		bcopy(completion, completecopy, (this->r - 2) * sizeof(int));

		this->augmentation_suceeded = this->satgraph->augment(noni, nonj, completecopy);

		LONG_T label = (((LONG_T) orbitIndex) << 24) | (LONG_T) completionIndex;
		node->curChild = label;

		this->stack.push_back(new SearchNode(label));

		return label;
	}

	printf(">> at the end of pushTo when I shouldn't be...\n");
	return -1;
}

/**
 * pop -- remove the current node and move up the tree.
 *
 * @return the label of the node after the pop.
 * 		This return value is used for validation purposes
 * 		to check proper implementation of push*() and pop().
 */
LONG_T SaturationManager::pop()
{
	if ( this->stack.size() == 0 )
	{
		return this->root->label;
	}

	SearchNode* node = this->stack.back();

	this->satgraph->pop();

	this->stack.pop_back();

	LONG_T label = node->label;

	delete node;

	return label;
}

/**
 * prune -- Perform a check to see if this node should be pruned.
 *
 * @return 0 if no prune should occur, 1 if prune should occur.
 */
int SaturationManager::prune()
{
	if ( this->augmentation_suceeded == false )
	{
		/* failed during the augmentation call */
		return 1;
	}

	if ( this->satgraph->isFeasible() == false )
	{
		/* Did not fail during the augmentation call,
		 * 	do an extra check on the resulting graph.
		 */
		return 1;
	}

	/* Now check that this augmentation is canonical */
	if ( this->satgraph->isCanonical() == false )
	{
		return 1;
	}
//
//	if ( this->searchDepth >= 1 && this->searchDepth <= 4 )
//	{
//		this->searchDepth = this->searchDepth + 1;
//		writeJob(stdout);
//		printf("-- isCanonical: ");
//		printString(this->satgraph);
//		this->searchDepth = this->searchDepth - 1;
//	}

	return 0;
}

/**
 * isSolution -- Perform a check to see if a solution exists
 * 		at this point.
 *
 * @return 0 if no solution is found, 1 if a solution is found.
 */
int SaturationManager::isSolution()
{
	return (this->satgraph->isFeasible() && this->satgraph->isUniquelySaturated()) ? 1 : 0;
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
char* SaturationManager::writeSolution()
{
	return this->satgraph->getString();
}

