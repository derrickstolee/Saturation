/*
 * Augmentation.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Augmentation.hpp"

/**
 * Constructor
 */
Augmentation::Augmentation( AUGMENT_TYPE type, int i, int j, int* completion, int r, int oldN )
{
	//	printf("\t\t\t--CREATING AUGMENTATION OF ");
	//
	//	switch ( type )
	//	{
	//		case NONEDGE:
	//			printf("NONEDGE");
	//			break;
	//
	//		case FILL:
	//			printf("FILL");
	//			break;
	//
	//		case ADDVERTEX:
	//			printf("ADDVERTEX");
	//			break;
	//	}
	//	printf(" TYPE.\n");

	this->type = type;

	this->i = i;
	this->j = j;

	/* copy! */
	if ( completion != 0 )
	{
		if ( type == NONEDGE )
		{
			this->completion = (int*) malloc((r - 2) * sizeof(int));
			bcopy(completion, this->completion, (r - 2) * sizeof(int));
		}
		else if ( type == FILL )
		{
			this->completion = (int*) malloc(i * sizeof(int));
			bcopy(completion, this->completion, i * sizeof(int));
		}
		else
		{
			this->completion = 0;
		}
	}
	else
	{
		this->completion = 0;
	}

	this->numOrbits = -1;
	this->orbits = 0;
	this->orbitLabels = 0;
	this->stabilizer = -1;
	this->numStabilizedOrbits = 0;
	this->completionOrbitReps = 0;
	this->canonicalLabels = 0;
	this->zeroOrbitLabels = 0;
<<<<<<< HEAD

	this->zeroEdges = 0;
	this->numZeroEdges = -1;

	this->numZeroGraphOrbits = -1;
	this->zeroGraphOrbitLabels = 0;
	this->zeroGraphCanonLabels = 0;
=======
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

	this->oldN = oldN;
}

/**
 * Destructor
 */
Augmentation::~Augmentation()
{
	//	printf("\t\t\t--DELETING AUGMENTATION OF ");
	//
	//	switch ( type )
	//	{
	//		case NONEDGE:
	//			printf("NONEDGE");
	//			break;
	//
	//		case FILL:
	//			printf("FILL");
	//			break;
	//
	//		case ADDVERTEX:
	//			printf("ADDVERTEX");
	//			break;
	//	}
	//	printf(" TYPE.\n");

	if ( this->completion != 0 )
	{
		free(this->completion);
		this->completion = 0;
	}

	if ( this->orbits != 0 )
	{
		for ( int i = 0; i < this->numOrbits; i++ )
		{
			if ( this->orbits[i] != 0 )
			{
				free(this->orbits[i]);
				this->orbits[i] = 0;
			}
		}

		free(this->orbits);
		this->orbits = 0;
	}

	if ( this->orbitLabels != 0 )
	{
		free(this->orbitLabels);
		this->orbitLabels = 0;
	}

	if ( this->canonicalLabels != 0 )
	{
		free(this->canonicalLabels);
		this->canonicalLabels = 0;
	}

	if ( this->zeroOrbitLabels != 0 )
	{
		free(this->zeroOrbitLabels);
		this->zeroOrbitLabels = 0;
	}

	if ( this->completionOrbitReps != 0 )
	{
		for ( int i = 0; i < this->numStabilizedOrbits; i++ )
		{
			if ( this->completionOrbitReps[i] != 0 )
			{
				free(this->completionOrbitReps[i]);
				this->completionOrbitReps[i] = 0;
			}
		}

		free(this->completionOrbitReps);
		this->completionOrbitReps = 0;
	}

	if ( this->zeroEdges != 0 )
	{
		free(this->zeroEdges);
		this->zeroEdges = 0;
		this->numZeroEdges = -1;
	}

	if ( this->zeroGraphOrbitLabels != 0 )
	{
		free(this->zeroGraphOrbitLabels);
		this->zeroGraphOrbitLabels = 0;
	}

	if ( this->zeroGraphCanonLabels != 0 )
	{
		free(this->zeroGraphCanonLabels);
		this->zeroGraphCanonLabels = 0;
	}

}
