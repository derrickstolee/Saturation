/*
 * SaturationDataSymmetry.hpp
 *
 *  Created on: Apr 18, 2011
 *      Author: stolee
 */

#ifndef SATURATIONDATASYMMETRY_HPP_
#define SATURATIONDATASYMMETRY_HPP_

#include "SaturationData.hpp"
#include "Augmentation.hpp"
#include "nausparse.h"
#include "gtools.h"

/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeUnassignedOrbits( SaturationData* satdata, Augmentation* augment );

/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeZeroGraphSymmetry( SaturationData* satdata, Augmentation* augment );

/**
 * Given a graph and an unassigned edge, compute all completion orbits
 * 	for that edge.
 *
 * BONUS: the completion orbits may be filtered by removing completion sets
 * 	where there is a 0-edge which would need to be a 1-edge.
 */
void computeStabilizedOrbits( int r, SaturationData* satdata, int pairindex, Augmentation* augment );

#endif /* SATURATIONDATASYMMETRY_HPP_ */
