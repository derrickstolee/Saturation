/*
 * connectivity.hpp
 *
 *  Created on: Dec 18, 2010
 *      Author: derrickstolee
 */

#ifndef CONNECTIVITY_HPP_
#define CONNECTIVITY_HPP_

#include "nausparse.h"

/**
 * is2connected
 *
 * Is the graph still 2-connected without the specified ear?
 *
 * @param g an undirected sparsegraph
 * @param v1 an index of a vertex
 * @param v2 an index of another vertex
 * @param without_ear an array of vertices, -1 terminated.
 */
bool is2connected( sparsegraph* g, int v1, int v2, char* without_ear );

/**
 * isCutEdge
 *
 * Is the given edge a cut edge?
 */
bool isCutEdge( sparsegraph* g, int v1, int v2 );

#endif /* CONNECTIVITY_HPP_ */
