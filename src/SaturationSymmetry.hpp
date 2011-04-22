/*
 * SaturationSymmetry.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: stolee
 */

#ifndef SATURATIONSYMMETRY_HPP_
#define SATURATIONSYMMETRY_HPP_

#include "SaturationGraph.hpp"
#include "Augmentation.hpp"
#include "nausparse.h"
#include "gtools.h"


/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeUnassignedOrbits(SaturationGraph* satgraph, Augmentation* augment);

/**
 * Given a graph and an unassigned edge, compute all completion orbits
 * 	for that edge.
 *
 * BONUS: the completion orbits may be filtered by removing completion sets
 * 	where there is a 0-edge which would need to be a 1-edge.
 */
void computeStabilizedOrbits(int r, SaturationGraph* satgraph, int pairindex, Augmentation* augment);

#endif /* SATURATIONSYMMETRY_HPP_ */
