/*
 * graphcompactor.h
 *
 *  Created on: Apr 18, 2011
 *      Author: stolee
 */

#ifndef GRAPHCOMPACTOR_H_
#define GRAPHCOMPACTOR_H_

#include "nausparse.h"

/**
 * compactsg
 *
 * Given a sparsegraph which may have gaps in its representation (the e array),
 * 	compact it to have no gaps.
 */
sparsegraph* compactsg( sparsegraph* g );

#endif /* GRAPHCOMPACTOR_H_ */
