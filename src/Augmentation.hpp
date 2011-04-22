/*
 * Augmentation.hpp
 *
 *  Created on: Apr 7, 2011
 *      Author: stolee
 */

#ifndef AUGMENTATION_HPP_
#define AUGMENTATION_HPP_

#include "nausparse.h"

typedef enum augment_type
{
	/**
	 * a NONEDGE augmentation
	 * takes a 2-edge and makes it 0-edge,
	 * and also fills in a completion.
	 */
	NONEDGE,

	/**
	 * a FILL augmentation
	 * takes all remaining 2-edges and makes them 1-edges.
	 * This is only a FINAL type of augmentation.
	 */
	FILL,

	/**
	 * an ADDVERTEX augmentation
	 * adds a vertex which has only 2-edges.
	 */
	ADDVERTEX
} AUGMENT_TYPE;

/**
 * Augmentation
 *
 * The Augmentation class stores data corresponding to
 * 	the augmentations for the SaturationGraph.
 * It contains enough information to augment as well as
 * 	delete.
 */
class Augmentation
{
	public:
		/**
		 * Constructor
		 */
		Augmentation( AUGMENT_TYPE type, int i, int j, int* completion, int r, int oldN );

		/**
		 * Destructor
		 */
		virtual ~Augmentation();

		/**
		 * type The Type of augmentation
		 */
		AUGMENT_TYPE type;

		/**
		 * i the first vertex in the non-edge pair.
		 *
		 * OR for FILL type, stores the number of filled edges.
		 */
		int i;

		/**
		 * j the second vertex in the non-edge pair.
		 */
		int j;

		/**
		 * completion The set of vertices in the completion.
		 *
		 * OR for FILL type, stores the indices of the filled edges.
		 */
		int* completion;

		/****** FOR LAYERED GRAPH SYMMETRY **********/

		/**
		 * numOrbits the number of 2-edge pair orbits.
		 */
		int numOrbits;

		/**
		 * orbits The list of -1-terminated orbits.
		 */
		int** orbits;

		/**
		 * orbitLabels Given a pair index, return the orbit index for that pair.
		 */
		int* orbitLabels;

		/**
		 * oneOrbitLabels Count orbits for 1-type edges.
		 */
		int* zeroOrbitLabels;

		/**
		 * canonicalLabels The canonical labeling of the vertices.
		 */
		int* canonicalLabels;

		/****** FOR LAYERED GRAPH STABILIZED SYMMETRY **********/

		/**
		 * stabilizer The pair index for the most recent stabilizer.
		 *
		 * -1 means it has not been set.
		 */
		int stabilizer;

		/**
		 * numStabilizedOrbits The number of completion orbits in the stabilizer.
		 */
		int numStabilizedOrbits;

		/**
		 * completionOrbitReps Representative completions for each orbit in the
		 * 	stabilizer.
		 *
		 * Each sub-array has (r-2) entries.
		 */
		int** completionOrbitReps;

		/**
		 * oldN The previous value of n, in order to undo vertex additions.
		 */
		int oldN;

		/**
		 * small_g the graph associated with this leve.
		 */
		sparsegraph* small_g;

		/**
		 * zeroEdges the set of implied zeroEdges
		 */
		int* zeroEdges;

		int numZeroEdges;

		/****** FOR ZERO GRAPH SYMMETRY **********/

		int numZeroGraphOrbits;

		int* zeroGraphOrbitLabels;

		int* zeroGraphCanonLabels;

};

#endif /* AUGMENTATION_HPP_ */
