/*
 * SaturationData.cpp
 *
 *  Created on: Apr 14, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Augmentation.hpp"
#include "nausparse.h"
#include "gtools.h"
#include "connectivity.hpp"
#include "translation.hpp"
#include "TreeSet.hpp"
#include "SaturationData.hpp"
#include "SaturationDataSymmetry.hpp"

#define SD_DEBUG_MODE
#define SD_DEBUG_PRIORITY 3

#ifdef SD_DEBUG_MODE
#define SD_DEBUG1(PRIORITY,STRING,ARG) if ( PRIORITY <= SD_DEBUG_PRIORITY) printf(STRING,ARG)
#define SD_DEBUG2(PRIORITY,STRING,ARG1,ARG2) if ( PRIORITY <= SD_DEBUG_PRIORITY) printf(STRING,ARG1,ARG2)
#define SD_DEBUG3(PRIORITY,STRING,ARG1,ARG2,ARG3) if ( PRIORITY <= SD_DEBUG_PRIORITY) printf(STRING,ARG1,ARG2,ARG3)
#else
#define SD_DEBUG1(PRIORITY,STRING,ARG)
#define SD_DEBUG2(PRIORITY,STRING,ARG1,ARG2)
#define SD_DEBUG3(PRIORITY,STRING,ARG1,ARG2,ARG3)
#endif

/**
 * snapshots -- The number of edge augmentations at the current depth.
 */
int SaturationData::initSnapshots()
{
	/* this is just a guess */
	this->depth = 0;
	this->maxDepth = (this->maxN * this->maxN) / 8;
	int snapshots_size = this->maxDepth * sizeof(int);
	this->snapshots = (int*) malloc(snapshots_size);
	this->snapshotorder = (int*) malloc(snapshots_size);
	this->augmentorder = (int*) malloc(snapshots_size);

	bzero(this->snapshots, snapshots_size);
	bzero(this->snapshotorder, snapshots_size);
	bzero(this->augmentorder, snapshots_size);

	int backup_size = 0;
#ifdef _USE_BACKUPS_
	backup_size += 5 * this->maxDepth * sizeof(char*);
	backup_size += this->maxDepth * sizeof(Set**);
	backup_size += this->maxDepth * sizeof(int*);
	backup_size += this->maxDepth * sizeof(Set*);
	this->bu_adjmat = (char**) malloc(this->maxDepth * sizeof(char*));
	this->bu_degseq = (char**) malloc(this->maxDepth * sizeof(char*));
	this->bu_adjlist = (char**) malloc(this->maxDepth * sizeof(char*));
	this->bu_set_common_neighs = (Set***) malloc(this->maxDepth * sizeof(Set**));
	this->bu_edges_in_set = (char**) malloc(this->maxDepth * sizeof(char*));
	this->bu_edges_to_set = (char**) malloc(this->maxDepth * sizeof(char*));
	this->bu_completion_sets = (int**) malloc(this->maxDepth * sizeof(int*));
	this->bu_two_edges = (Set**) malloc(this->maxDepth * sizeof(Set*));

	bzero(this->bu_adjmat, this->maxDepth * sizeof(char*));
	bzero(this->bu_degseq, this->maxDepth * sizeof(char*));
	bzero(this->bu_adjlist, this->maxDepth * sizeof(char*));
	bzero(this->bu_set_common_neighs, this->maxDepth * sizeof(Set**));
	bzero(this->bu_edges_in_set, this->maxDepth * sizeof(char*));
	bzero(this->bu_edges_to_set, this->maxDepth * sizeof(char*));
	bzero(this->bu_completion_sets, this->maxDepth * sizeof(int*));
	bzero(this->bu_two_edges, this->maxDepth * sizeof(Set*));
#endif

	return 3 * snapshots_size + backup_size;
}

void SaturationData::testSnapshotSize()
{
	if ( this->depth >= this->maxDepth )
	{
		this->maxDepth = 2 * this->maxDepth;
		int snapshots_size = this->maxDepth * sizeof(int);
		this->snapshots = (int*) realloc(this->snapshots, snapshots_size);
		this->snapshotorder = (int*) realloc(this->snapshotorder, snapshots_size);
		this->augmentorder = (int*) realloc(this->augmentorder, snapshots_size);

		for ( int i = this->depth; i < this->maxDepth; i++ )
		{
			this->snapshots[i] = 0;
			this->snapshotorder[i] = 0;
			this->augmentorder[i] = 0;
		}

#ifdef _USE_BACKUPS_
		this->bu_adjmat = (char**) realloc(this->bu_adjmat, this->maxDepth * sizeof(char*));
		this->bu_degseq = (char**) realloc(this->bu_degseq, this->maxDepth * sizeof(char*));
		this->bu_adjlist = (char**) realloc(this->bu_adjlist, this->maxDepth * sizeof(char*));
		this->bu_set_common_neighs = (Set***) realloc(this->bu_set_common_neighs, this->maxDepth * sizeof(Set**));
		this->bu_edges_in_set = (char**) realloc(this->bu_edges_in_set, this->maxDepth * sizeof(char*));
		this->bu_edges_to_set = (char**) realloc(this->bu_edges_to_set, this->maxDepth * sizeof(char*));
		this->bu_completion_sets = (int**) realloc(this->bu_completion_sets, this->maxDepth * sizeof(int*));
		this->bu_two_edges = (Set**) realloc(this->bu_two_edges, this->maxDepth * sizeof(Set*));

		for ( int i = this->depth; i < this->maxDepth; i++ )
		{
			this->bu_adjmat[i] = 0;
			this->bu_degseq[i] = 0;
			this->bu_adjlist[i] = 0;
			this->bu_set_common_neighs[i] = 0;
			this->bu_edges_in_set[i] = 0;
			this->bu_edges_to_set[i] = 0;
			this->bu_completion_sets[i] = 0;
			this->bu_two_edges[i] = 0;
		}
#endif
	}
}

void SaturationData::cleanSnapshots()
{
	if ( this->snapshots != 0 )
	{
		free(this->snapshots);
		this->snapshots = 0;
	}

	if ( this->snapshotorder != 0 )
	{
		free(this->snapshotorder);
		this->snapshotorder = 0;
	}

	if ( this->augmentorder != 0 )
	{
		free(this->augmentorder);
		this->augmentorder = 0;
	}
}

/**
 * augmentations -- An array of augmentations for
 * 	storing orbis.
 */
int SaturationData::initAugmentations()
{
	this->numAugmentations = 0;
	this->augmentationsSize = this->maxDepth;
	this->augmentations = (Augmentation**) malloc(this->maxDepth * sizeof(Augmentation*));
	bzero(this->augmentations, this->maxDepth * sizeof(Augmentation*));

	return this->maxDepth * sizeof(Augmentation*);
}

void SaturationData::testAugmentationSize()
{
	if ( this->numAugmentations >= this->augmentationsSize )
	{
		this->augmentationsSize = this->augmentationsSize * 2;
		this->augmentations = (Augmentation**) realloc(this->augmentations, this->augmentationsSize
				* sizeof(Augmentation*));

		for ( int i = this->numAugmentations; i < this->augmentationsSize; i++ )
		{
			this->augmentations[i] = 0;
		}
	}
}

void SaturationData::cleanAugmentations()
{
	if ( this->augmentations != 0 )
	{
		for ( int i = 0; i < this->numAugmentations; i++ )
		{
			if ( this->augmentations[i] != 0 )
			{
				delete this->augmentations[i];
				this->augmentations[i] = 0;
			}
		}

		free(this->augmentations);
		this->augmentations = 0;
	}
}

void SaturationData::addAugmentation(Augmentation* augment)
{
	this->testAugmentationSize();

	this->augmentations[this->numAugmentations] = augment;
	(this->numAugmentations)++;
}

void SaturationData::removeAugmentation()
{
	(this->numAugmentations)--;
	delete this->augmentations[this->numAugmentations];
	this->augmentations[this->numAugmentations] = 0;
}

/**
 * layer_graph is the two-layered graph which places edges in each layer
 * 	for type 0 and type 1 edges.
 * It is built to grow as edges and vertices are added.
 * The layers alternate: normal vertex 0 maps to 0 and 1, vert 1 is 2 and 3, etc.
 *
 * The graph is modified by the following operations:
 *  addEdgeToLayerGraph( i,  j,  type)
 *  removeEdgeFromLayerGraph(i, j, type)
 *  addVertexToLayerGraph()
 *  removeVertexFromLayerGraph()
 */
int SaturationData::initLayerGraph()
{
	int layer_graph_size = sizeof(sparsegraph);
	this->layer_graph = (sparsegraph*) malloc(layer_graph_size);

	SG_INIT(*(this->layer_graph));

	this->layer_graph->nv = 2 * this->n;
	this->layer_graph->vlen = 2 * this->maxN;
	this->layer_graph->dlen = 2 * this->maxN;

	int vert_array_size = 2 * this->maxN * sizeof(int);
	this->layer_graph->v = (int*) malloc(vert_array_size);
	this->layer_graph->d = (int*) malloc(vert_array_size);

	for ( int i = 0; i < 2 * this->maxN; i++ )
	{
		this->layer_graph->v[i] = this->maxN * i;
		this->layer_graph->d[i] = 1;
	}

	/* max degree is maxN over 2maxN vertices */
	this->layer_graph->elen = 2 * this->maxN * this->maxN;
	int edge_array_size = 2 * this->maxN * this->maxN * sizeof(int);
	this->layer_graph->e = (int*) malloc(edge_array_size);

	for ( int i = 0; i < this->layer_graph->elen; i++ )
	{
		this->layer_graph->e[i] = -1;
	}

	/* add crossing edges */
	for ( int i = 0; i < this->maxN; i++ )
	{
		this->layer_graph->e[this->layer_graph->v[2 * i]] = 2 * i + 1;
		this->layer_graph->e[this->layer_graph->v[2 * i + 1]] = 2 * i;
	}

	/* TODO: fill in layer_graph with current graph information */

	return layer_graph_size + 2 * vert_array_size + edge_array_size;
}

void SaturationData::cleanLayerGraph()
{
	if ( this->layer_graph != 0 )
	{
		SG_FREE(*(this->layer_graph));
		free(this->layer_graph);
		this->layer_graph = 0;
	}
}

void SaturationData::addEdgeToLayerGraph(int i, int j, char type)
{
	/* This will only add the edge i->j */
	/* For j->i, it will be called again */
	int offset = 0;
	if ( type == 1 )
	{
		offset = 1;
	}

	SD_DEBUG3(10, "--[SaturationData::addEdgeToLayerGraph(%d,%d,%d)]\n", i,j,type);

	int vi = this->layer_graph->v[2 * i + offset];
	int di = this->layer_graph->d[2 * i + offset];

	/* add the edge */
	this->layer_graph->e[vi + di] = 2 * j + offset;

	/* increase the degree */
	this->layer_graph->d[2 * i + offset] = di + 1;
}

void SaturationData::removeEdgeFromLayerGraph(int i, int j, char type)
{
	/* This will only remove the edge i->j */
	/* For j->i, it will be called again */
	int offset = 0;
	if ( type == 1 )
	{
		offset = 1;
	}

	int vi = this->layer_graph->v[2 * i + offset];
	int di = this->layer_graph->d[2 * i + offset];

	/* decrease the degree */
	this->layer_graph->d[2 * i + offset] = di - 1;

	/* remove the edge */
	if ( this->layer_graph->e[vi + di - 1] != 2 * j + offset )
	{
		SD_DEBUG3(3,"--[SaturationData::removeEdgeFromLayerGraph(%d,%d,%d)] This neighbor does not agree!\n", i,j,type);
	}

	this->layer_graph->e[vi + di - 1] = -1;
}

void SaturationData::addVertexToLayerGraph()
{
	/* add the vertex! */
	this->layer_graph->nv = this->layer_graph->nv + 2;
}

void SaturationData::removeVertexFromLayerGraph()
{
	/*  remove the vertex! */
	this->layer_graph->nv = this->layer_graph->nv - 2;
}

/**
 * one_graph is the standard graph of just the 1-type edges.
 *
 * The graph is built to grow with vertex and edge additions.
 *
 * The graph is modified by the methods
 * 	addEdgeToOneGraph(i,j)
 *  removeEdgeFromOneGraph(i,j)
 */
int SaturationData::initOneGraph()
{
	int one_graph_size = sizeof(sparsegraph);
	this->one_graph = (sparsegraph*) malloc(one_graph_size);

	SG_INIT(*(this->one_graph));

	this->one_graph->nv = this->n;
	this->one_graph->vlen = this->maxN;
	this->one_graph->dlen = this->maxN;

	int vert_array_size = this->maxN * sizeof(int);
	this->one_graph->v = (int*) malloc(vert_array_size);
	this->one_graph->d = (int*) malloc(vert_array_size);

	for ( int i = 0; i < this->maxN; i++ )
	{
		this->one_graph->v[i] = i * (this->maxN - 1);
		this->one_graph->d[i] = 0;
	}

	/* max degree is maxN-1 over maxN vertices */
	this->one_graph->elen = this->maxN * (this->maxN - 1);
	int edge_array_size = this->maxN * (this->maxN - 1) * sizeof(int);
	this->one_graph->e = (int*) malloc(edge_array_size);

	for ( int i = 0; i < this->one_graph->elen; i++ )
	{
		this->one_graph->e[i] = -1;
	}

	/* TODO: fill in one_graph with current graph information */

	return one_graph_size + 2 * vert_array_size + edge_array_size;
}

void SaturationData::cleanOneGraph()
{
	if ( this->one_graph != 0 )
	{
		SG_FREE(*(this->one_graph));
		free(this->one_graph);
		this->one_graph = 0;
	}
}

void SaturationData::addEdgeToOneGraph(int i, int j)
{
	int vi = this->one_graph->v[i];
	int di = this->one_graph->d[i];

	this->one_graph->e[vi + di] = j;
	this->one_graph->d[i] = di + 1;
}

void SaturationData::removeEdgeFromOneGraph(int i, int j)
{
	int vi = this->one_graph->v[i];
	int di = this->one_graph->d[i];

	if ( this->one_graph->e[vi + di - 1] != j )
	{
		SD_DEBUG2 (3,"--[SaturationData::removeEdgeFromOneGraph(%d,%d)] This is not the right neighbor!\n",i,j);
	}
	this->one_graph->e[vi + di - 1] = -1;
	this->one_graph->d[i] = di - 1;
}

void SaturationData::addVertexToOneGraph()
{
	this->one_graph->nv = this->one_graph->nv + 1;
}

void SaturationData::removeVertexFromOneGraph()
{
	this->one_graph->nv = this->one_graph->nv - 1;
}

/**
 * zero_graph is the standard graph of just the 0-type edges.
 *
 * The graph is built to grow with vertex and edge additions.
 *
 * The graph is modified by the methods
 * 	addEdgeToZeroGraph(i,j)
 *  removeEdgeFromZeroGraph(i,j)
 */
int SaturationData::initZeroGraph()
{
	int zero_graph_size = sizeof(sparsegraph);
	this->zero_graph = (sparsegraph*) malloc(zero_graph_size);

	SG_INIT(*(this->zero_graph));

	this->zero_graph->nv = this->n;
	this->zero_graph->vlen = this->maxN;
	this->zero_graph->dlen = this->maxN;

	int vert_array_size = this->maxN * sizeof(int);
	this->zero_graph->v = (int*) malloc(vert_array_size);
	this->zero_graph->d = (int*) malloc(vert_array_size);

	for ( int i = 0; i < this->maxN; i++ )
	{
		this->zero_graph->v[i] = i * (this->maxN - 1);
		this->zero_graph->d[i] = 0;
	}

	/* max degree is maxN-1 over maxN vertices */
	this->zero_graph->elen = this->maxN * (this->maxN - 1);
	int edge_array_size = this->maxN * (this->maxN - 1) * sizeof(int);
	this->zero_graph->e = (int*) malloc(edge_array_size);

	for ( int i = 0; i < this->zero_graph->elen; i++ )
	{
		this->zero_graph->e[i] = -1;
	}

	/* TODO: fill in zero_graph with current graph information */

	return zero_graph_size + 2 * vert_array_size + edge_array_size;
}

void SaturationData::cleanZeroGraph()
{
	if ( this->zero_graph != 0 )
	{
		SG_FREE(*(this->zero_graph));
		free(this->zero_graph);
		this->zero_graph = 0;
	}
}

void SaturationData::addEdgeToZeroGraph(int i, int j)
{
	int vi = this->zero_graph->v[i];
	int di = this->zero_graph->d[i];

	this->zero_edges->add(indexOf(i, j));

	this->zero_graph->e[vi + di] = j;
	this->zero_graph->d[i] = di + 1;
}

void SaturationData::removeEdgeFromZeroGraph(int i, int j)
{
	int vi = this->zero_graph->v[i];
	int di = this->zero_graph->d[i];

	this->zero_edges->remove(indexOf(i, j));

	if ( this->zero_graph->e[vi + di - 1] != j )
	{
		SD_DEBUG2 (3,"--[SaturationData::removeEdgeFromZeroGraph(%d,%d)] This is not the right neighbor!\n",i,j);
	}
	this->zero_graph->e[vi + di - 1] = -1;
	this->zero_graph->d[i] = di - 1;
}

void SaturationData::addVertexToZeroGraph()
{
	this->zero_graph->nv = this->zero_graph->nv + 1;
}

void SaturationData::removeVertexFromZeroGraph()
{
	this->zero_graph->nv = this->zero_graph->nv - 1;
}

/**
 * The adjacency matrix stores edge types of 0, 1, 2, and 3.
 *
 * 0: A known non-edge. It MUST have a completion.
 * 1: A known edge.
 * 2: An unknown edge.
 * 3: An edge pair which uses a vertex that is not currently in the graph.
 */
int SaturationData::initAdjacencyMatrix()
{
	int adjmat_size = nChooseK(this->maxN, 2);

	this->adjmat = (char*) malloc(adjmat_size);

	for ( int i = 0; i < adjmat_size; i++ )
	{
		this->adjmat[i] = 3;
	}

	return adjmat_size;
}

void SaturationData::cleanAdjacencyMatrix()
{
	if ( this->adjmat != 0 )
	{
		free(this->adjmat);
		this->adjmat = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyAdjacencyMatrix()
{
	if ( this->bu_adjmat[this->depth] == 0 )
	{
		this->bu_adjmat[this->depth] = (char*) malloc(nChooseK(this->maxN, 2));
	}

	int nC2 = nChooseK(this->maxN, 2);
	for ( int i = 0; i < nC2; i++ )
	{
		this->bu_adjmat[this->depth][i] = this->adjmat[i];
	}
}

bool SaturationData::compareAdjacencyMatrix()
{
	bool result = true;
	int num_diff = 0;

	int nC2 = nChooseK(this->maxN, 2);
	for ( int i = 0; i < nC2; i++ )
	{
		if ( this->bu_adjmat[this->depth][i] != this->adjmat[i] )
		{
			result = false;
			num_diff++;
		}
	}

	if ( !result )
	{
		SD_DEBUG2(1,"--[SaturationData::compareAdjacencyMatrix()] %d different entries at depth %d.\n", num_diff, this->depth);
		return false;
	}

	return true;
}
#endif

/**
 * The degree sequence stores two types of degrees: 0 and 1.
 *
 * It is a 2-dimensional array, maxN x 2, packed into a single character list.
 *
 * The data type is char since we will never reach close to even 50 vertices.
 *
 * The getDegree(i,type) and setDegree(i,type) methods modify this sequence.
 */
int SaturationData::initDegreeSequence()
{
	int degseq_size = 2 * this->maxN;

	this->degseq = (char*) malloc(degseq_size);
	bzero(this->degseq, degseq_size);

	return degseq_size;
}

void SaturationData::cleanDegreeSequence()
{
	if ( this->degseq != 0 )
	{
		free(this->degseq);
		this->degseq = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyDegreeSequences()
{
	if ( this->bu_degseq[this->depth] == 0 )
	{
		this->bu_degseq[this->depth] = (char*) malloc(this->maxN * 2);
	}

	for ( int i = 0; i < this->maxN * 2; i++ )
	{
		this->bu_degseq[this->depth][i] = this->degseq[i];
	}
}

bool SaturationData::compareDegreeSequences()
{
	bool result = true;
	int num_diff = 0;

	for ( int i = 0; i < this->maxN * 2; i++ )
	{
		if ( this->bu_degseq[this->depth][i] != this->degseq[i] )
		{
			result = false;
			num_diff++;
		}
	}

	if ( !result )
	{
		SD_DEBUG2(1,"--[SaturationData::compareDegreeSequences()] %d different at depth %d.\n", num_diff, this->depth);
	}

	return result;
}
#endif

char SaturationData::getDegree(int i, int type)
{
	if ( type < 0 || type > 1 || i < 0 || i >= this->n )
	{
		return 0;
	}

	return this->degseq[2 * i + type];
}

void SaturationData::setDegree(int i, int type, char value)
{
	if ( type < 0 || type > 1 || i < 0 || i >= this->n )
	{
		return;
	}

	this->degseq[2 * i + type] = value;
}

/**
 * the adjacency lists store the (0/1)-neighborhoods of each vertex,
 * 	in order of augmentations.
 *
 * It is a 3-dimensional array, maxN x 2 x maxN, packed into a single list
 *
 * The methods addNeighbor(i,j,type) and removeNeighbor(i,j,type) modify
 * 	the lists, the degrees, AND the adjacency matrix.
 */
int SaturationData::initAdjacencyList()
{
	int adjlist_size = 2 * this->maxN * this->maxN;
	this->adjlist = (char*) malloc(adjlist_size);
	bzero(this->adjlist, adjlist_size);

	return adjlist_size;
}

void SaturationData::cleanAdjacencyList()
{
	if ( this->adjlist != 0 )
	{
		free(this->adjlist);
		this->adjlist = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyAdjacencyLists()
{
	int adjlist_size = 2 * this->maxN * this->maxN;
	if ( this->bu_adjlist[this->depth] == 0 )
	{
		this->bu_adjlist[this->depth] = (char*) malloc(adjlist_size);
	}

	for ( int i = 0; i < adjlist_size; i++ )
	{
		this->bu_adjlist[this->depth][i] = this->adjlist[i];
	}
}

bool SaturationData::compareAdjacencyLists()
{
	int num_diff = 0;

	int adjlist_size = 2 * this->maxN * this->maxN;

	for ( int i = 0; i < adjlist_size; i++ )
	{
		if ( this->bu_adjlist[this->depth][i] != this->adjlist[i] )
		{
			num_diff++;
		}
	}

	if ( num_diff > 0 )
	{
		SD_DEBUG2(1,"--[SaturationDatat::compareAdjacencyLists()] %d different at depth %d.\n", num_diff, this->depth);
		return false;
	}
	return true;
}
#endif

int SaturationData::getBasePosition(int i, char type)
{
	if ( type < 2 )
	{
		return 2 * i * this->maxN + type * this->maxN;
	}
	else
	{
		/* a bad type! */
		return -1;
	}
}

bool SaturationData::addNeighbor(int i, int j, char type)
{
	if ( type >= 2 )
	{
		return false;
	}

	int ideg = this->getDegree(i, type);

	this->adjlist[this->getBasePosition(i, type) + ideg] = j;
	this->setDegree(i, type, ideg + 1);

	/* also modify the layered graph */
	this->addEdgeToLayerGraph(i, j, type);
	/* these modifications will happen no matter what, and will be undone */
	if ( type == 0 )
	{
		this->addEdgeToZeroGraph(i, j);
	}
	else if ( type == 1 )
	{
		this->addEdgeToOneGraph(i, j);
	}

	return true;
}

void SaturationData::removeNeighbor(int i, int j, char type)
{
	int ideg = this->getDegree(i, type);

	if ( ideg == 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::removeNeighbor(%d,%d,%d)] Trying to remove neighbor when degree is 0.\n", i, j,
				(int) type);
	}
	else if ( this->adjlist[this->getBasePosition(i, type) + ideg - 1] == j )
	{
		this->adjlist[this->getBasePosition(i, type) + ideg - 1] = 0;
		this->setDegree(i, type, ideg - 1);

		/* also modify the layered graph */
		this->removeEdgeFromLayerGraph(i, j, type);
		/* these modifications will happen no matter what, and will be undone */
		if ( type == 0 )
		{
			this->removeEdgeFromZeroGraph(i, j);
		}
		else if ( type == 1 )
		{
			this->removeEdgeFromOneGraph(i, j);
		}

	}
	else
	{
		SD_DEBUG3(3,"--[SaturationData::removeNeighbor(%d,%d,%d)] Trying to remove wrong neighbor..\n", i, j,
				(int)type);
	}

}

/**
 * set_common_neighs stores sets of common neighbors for each set
 * of order s (s from 2 to r-2)
 *
 * set_common_base_index stores the base index for each size.
 *
 * So, the array is (N choose 2) + (N choose 3) + ... + (N choose r-2)
 */
int SaturationData::initSetCommonNeighs()
{
	int base_index_size = (this->getMaxCommonNeighsSetSize()) * sizeof(int);
	this->set_common_base_index = (int*) malloc(base_index_size);

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxCommonNeighsSetSize(); i++ )
	{
		this->set_common_base_index[i - 2] = base_index;
		base_index += nChooseK(this->maxN, i);
	}

	int set_common_size = base_index * sizeof(Set*);
	this->set_common_neighs = (Set**) malloc(set_common_size);

	for ( int i = 0; i < base_index; i++ )
	{
		this->set_common_neighs[i] = new TreeSet();
	}

	return base_index_size + set_common_size + base_index * sizeof(TreeSet);
}

void SaturationData::cleanSetCommonNeighs()
{
	if ( this->set_common_neighs != 0 )
	{
		for ( int i = 2; i <= this->getMaxCommonNeighsSetSize(); i++ )
		{
			int base = this->set_common_base_index[i - 2];
			int iCount = nChooseK(this->maxN, i);
			for ( int j = 0; j < iCount; j++ )
			{
				if ( this->set_common_neighs[base + j] != 0 )
				{
					delete this->set_common_neighs[base + j];
					this->set_common_neighs[base + j] = 0;
				}
			}
		}

		free(this->set_common_neighs);
		this->set_common_neighs = 0;

		free(this->set_common_base_index);
		this->set_common_base_index = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyCommonNeighbors()
{
	int base_index = 0;
	for ( int i = 2; i <= this->getMaxCommonNeighsSetSize(); i++ )
	{
		base_index += nChooseK(this->maxN, i);
	}

	if ( this->bu_set_common_neighs[this->depth] == 0 )
	{
		int set_common_size = base_index * sizeof(Set*);
		this->bu_set_common_neighs[this->depth] = (Set**) malloc(set_common_size);
		bzero(this->bu_set_common_neighs[this->depth], set_common_size);
	}

	for ( int i = 0; i < base_index; i++ )
	{
		if ( this->bu_set_common_neighs[this->depth][i] == 0 )
		{
			this->bu_set_common_neighs[this->depth][i] = new TreeSet();
		}

		this->bu_set_common_neighs[this->depth][i]->clear();

		this->set_common_neighs[i]->resetIterator();
		while ( this->set_common_neighs[i]->hasNext() )
		{
			this->bu_set_common_neighs[this->depth][i]->add(this->set_common_neighs[i]->next());
		}
	}
}

bool SaturationData::compareCommonNeighbors()
{
	int num_diff = 0;

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxCommonNeighsSetSize(); i++ )
	{
		base_index += nChooseK(this->maxN, i);
	}

	for ( int i = 0; i < base_index; i++ )
	{
		Set* set_inplace = this->set_common_neighs[i];
		Set* set_copy = this->bu_set_common_neighs[this->depth][i];

		if ( set_inplace->size() != set_copy->size() )
		{
			SD_DEBUG3(1,"--[SaturationData::compareCommonNeighbors()] Depth %d size difference IN PLACE %d != COPY %d.\n", this->depth, set_inplace->size(), set_copy->size());
			num_diff++;
		}

		set_inplace->resetIterator();
		while ( set_inplace->hasNext() )
		{
			int next = set_inplace->next();
			if ( set_copy->contains(next) == 0 )
			{
				num_diff++;
			}
		}
	}

	if ( num_diff > 0 )
	{
		SD_DEBUG2(1,"--[SaturationData::compareCommonNeighbors()] %d diff at depth %d.\n", num_diff, this->depth);
		return false;
	}

	return true;
}
#endif

int SaturationData::getMaxCommonNeighsSetSize()
{
	return this->r - 1;
}

Set* SaturationData::getSetCommonNeighs(int size, int index)
{
	if ( size < 2 )
	{
		return 0;
	}
	else if ( size > this->getMaxCommonNeighsSetSize() )
	{
		return 0;
	}

	/* use this->n to bound, as we should not ask outside of the current graph */
	if ( index < 0 || index >= nChooseK(this->n, size) )
	{
		SD_DEBUG2(3,"--[SaturationData::getSetCommonNeighs(%d,%d)] Index out of bounds.\n", size, index);
		return 0;
	}

	int base = this->set_common_base_index[size - 2];

	return this->set_common_neighs[base + index];
}

bool SaturationData::addSetCommonNeigh(int size, int set_index, int i)
{
	if ( size < 2 )
	{
		return false;
	}
	else if ( size > this->getMaxCommonNeighsSetSize() )
	{
		return false;
	}

	if ( i < 0 || i >= this->n )
	{
		return false;
	}

	/* use this->n to bound, as we should not ask outside of the current graph */
	if ( set_index < 0 || set_index >= nChooseK(this->n, size) )
	{
		SD_DEBUG3(3,"--[SaturationData::addSetCommonNeigh(%d,%d,%d)] Index out of bounds.\n", size, set_index, i);
		return false;
	}

	int base = this->set_common_base_index[size - 2];

	Set* neighs = this->set_common_neighs[base + set_index];

	if ( neighs->contains(i) == 1 )
	{
		SD_DEBUG3(3,"--[SaturationData::addSetCommonNeigh(%d,%d,%d)] Already contains this vertex!\n", size, set_index, i);
		return false;
	}

	neighs->add(i);

	return true;
}

void SaturationData::removeSetCommonNeigh(int size, int set_index, int i)
{
	if ( size < 2 )
	{
		return;
	}
	else if ( size > this->getMaxCommonNeighsSetSize() )
	{
		return;
	}

	if ( i < 0 || i >= this->n )
	{
		return;
	}

	/* use this->n to bound, as we should not ask outside of the current graph */
	if ( set_index < 0 || set_index >= nChooseK(this->n, size) )
	{
		SD_DEBUG3(3,"--[SaturationData::removeSetCommonNeigh(%d,%d,%d)] Index out of bounds.\n", size, set_index, i);
		return;
	}

	int base = this->set_common_base_index[size - 2];

	Set* neighs = this->set_common_neighs[base + set_index];

	if ( neighs->contains(i) == 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::addSetCommonNeigh(%d,%d,%d)] Does not contain this vertex!\n", size, set_index, i);
		return;
	}

	neighs->remove(i);
}

/**
 * edgesinset stores the number of edges within a set of order 2...(r)
 *
 * This is modified by the methods addEdgeInSet(size, index) and
 * 	removeEdgeInSet(size, index).
 * Returns FALSE if a K_r is created.
 */

int SaturationData::initEdgesInSet()
{
	int base_index_size = (this->getMaxEdgesInSetSize()) * sizeof(int);
	this->edges_in_set_base_index = (int*) malloc(base_index_size);

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesInSetSize(); i++ )
	{
		this->edges_in_set_base_index[i - 2] = base_index;
		base_index += nChooseK(this->maxN, i);
	}

	int edges_in_size = base_index;
	this->edges_in_set = (char*) malloc(edges_in_size);

	for ( int i = 0; i < base_index; i++ )
	{
		this->edges_in_set[i] = 0;
	}

	return base_index_size + edges_in_size;
}

void SaturationData::cleanEdgesInSet()
{
	if ( this->edges_in_set != 0 )
	{
		free(this->edges_in_set);
		this->edges_in_set = 0;

		free(this->edges_in_set_base_index);
		this->edges_in_set_base_index = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyEdgesInSet()
{
	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesInSetSize(); i++ )
	{
		base_index += nChooseK(this->maxN, i);
	}

	if ( this->bu_edges_in_set[this->depth] == 0 )
	{
		this->bu_edges_in_set[this->depth] = (char*) malloc(base_index);
	}

	for ( int i = 0; i < base_index; i++ )
	{
		this->bu_edges_in_set[this->depth][i] = this->edges_in_set[i];
	}
}

bool SaturationData::compareEdgesInSet()
{
	int num_diff = 0;

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesInSetSize(); i++ )
	{
		base_index += nChooseK(this->maxN, i);
	}

	for ( int i = 0; i < base_index; i++ )
	{
		if ( this->bu_edges_in_set[this->depth][i] != this->edges_in_set[i] )
		{
			num_diff++;
		}
	}

	if ( num_diff > 0 )
	{
		SD_DEBUG2(1,"--[SaturationData::compareEdgesInSet()] %d diff at depth %d.\n", num_diff, this->depth);
		return false;
	}

	return true;
}
#endif

int SaturationData::getMaxEdgesInSetSize()
{
	return this->r;
}

bool SaturationData::addEdgeInSet(int size, int index)
{
	if ( size < 2 || size > this->getMaxEdgesInSetSize() )
	{
		SD_DEBUG2(3,"--[SaturationData::addEdgeInSet(%d,%d)] Size out of bounds.\n", size, index);
		return false;
	}

	if ( index < 0 || index >= nChooseK(this->maxN, size) )
	{
		SD_DEBUG2(3,"--[SaturationData::addEdgeInSet(%d,%d)] Index out of bounds.\n", size, index);
		return false;
	}

	int base = this->edges_in_set_base_index[size - 2];

	this->edges_in_set[base + index] = this->edges_in_set[base + index] + 1;

	return true;
}

void SaturationData::removeEdgeInSet(int size, int index)
{
	if ( size < 2 || size > this->getMaxEdgesInSetSize() )
	{
		SD_DEBUG2(3,"--[SaturationData::removeEdgeInSet(%d,%d)] Size out of bounds.\n", size, index);
		return;
	}

	if ( index < 0 || index >= nChooseK(this->maxN, size) )
	{
		SD_DEBUG2(3,"--[SaturationData::removeEdgeInSet(%d,%d)] Index out of bounds.\n", size, index);
		return;
	}

	int base = this->edges_in_set_base_index[size - 2];

	if ( this->edges_in_set[base + index] <= 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::removeEdgeInSet(%d,%d)] Count (%d) too low to remove.\n", size, index,this->edges_in_set[base+index] );
		return;
	}

	this->edges_in_set[base + index] = this->edges_in_set[base + index] - 1;

}

bool SaturationData::isSetFull(int size, int index)
{
	if ( size < 2 || size > this->getMaxEdgesInSetSize() )
	{
		SD_DEBUG2(3,"--[SaturationData::isSetFull(%d,%d)] Size out of bounds.\n", size, index);
		return false;
	}

	if ( index < 0 || index >= nChooseK(this->maxN, size) )
	{
		SD_DEBUG2(3,"--[SaturationData::isSetFull(%d,%d)] Index out of bounds.\n", size, index);
		return false;
	}

	int base = this->edges_in_set_base_index[size - 2];

	if ( this->edges_in_set[base + index] >= nChooseK(size, 2) )
	{
		return true;
	}

	return false;
}

/**
 * edges_to_set stores the number of edges between a set of order 2...(r-1)
 * 	and another vertex i
 */
int SaturationData::initEdgesToSet()
{
	int edges_base_size = (this->getMaxEdgesToSetSize() - 1) * sizeof(int);
	this->edges_to_set_base_index = (int*) malloc(edges_base_size);

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesToSetSize(); i++ )
	{
		this->edges_to_set_base_index[i - 2] = base_index;
		base_index += nChooseK(this->maxN, i) * this->maxN;
	}

	int edges_to_set_size = base_index * sizeof(char);
	this->edges_to_set = (char*) malloc(edges_to_set_size);
	bzero(this->edges_to_set, edges_to_set_size);

	return edges_base_size + edges_to_set_size;
}

void SaturationData::cleanEdgesToSet()
{
	if ( this->edges_to_set != 0 )
	{
		free(this->edges_to_set);
		this->edges_to_set = 0;
		free(this->edges_to_set_base_index);
		this->edges_to_set_base_index = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyEdgesToSet()
{
	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesToSetSize(); i++ )
	{
		this->edges_to_set_base_index[i - 2] = base_index;
		base_index += nChooseK(this->maxN, i) * this->maxN;
	}

	if ( this->bu_edges_to_set[this->depth] == 0 )
	{
		this->bu_edges_to_set[this->depth] = (char*) malloc(base_index);
	}

	for ( int i = 0; i < base_index; i++ )
	{
		this->bu_edges_to_set[this->depth][i] = this->edges_to_set[i];
	}
}

bool SaturationData::compareEdgesToSet()
{
	int num_diff = 0;

	int base_index = 0;
	for ( int i = 2; i <= this->getMaxEdgesToSetSize(); i++ )
	{
		this->edges_to_set_base_index[i - 2] = base_index;
		base_index += nChooseK(this->maxN, i) * this->maxN;
	}

	for ( int i = 0; i < base_index; i++ )
	{
		if ( this->bu_edges_to_set[this->depth][i] != this->edges_to_set[i] )
		{
			num_diff++;
		}
	}

	if ( num_diff > 0 )
	{
		SD_DEBUG2(1,"--[SaturationData::compareEdgesToSet()] %d different at depth %d.\n", num_diff, this->depth);
		return false;
	}

	return true;
}
#endif

int SaturationData::getMaxEdgesToSetSize()
{
	return this->r - 1;
}

int SaturationData::getEdgesToSetPosition(int size, int index, int i)
{
	if ( size < 2 || size > this->getMaxEdgesToSetSize() || index < 0 || index >= nChooseK(this->maxN, size) || i < 0
			|| i >= this->n )
	{
		SD_DEBUG3(3,"--[SaturationData::getEdgesToSetPosition(%d,%d,%d)] A bad set of input.\n", size, index, i);
		return -1;
	}

	return this->edges_to_set_base_index[size - 2] + this->maxN * index + i;
}

bool SaturationData::addEdgeToSet(int size, int index, int i)
{
	int pos = this->getEdgesToSetPosition(size, index, i);

	if ( pos < 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::addEdgeToSet(%d,%d,%d)] A bad set of input.\n", size, index, i);
		return false;
	}

	this->edges_to_set[pos] = this->edges_to_set[pos] + 1;

	return true;
}

void SaturationData::removeEdgeToSet(int size, int index, int i)
{
	int pos = this->getEdgesToSetPosition(size, index, i);

	if ( pos < 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::removeEdgeToSet(%d,%d,%d)] A bad set of input.\n", size, index, i);
		return;
	}

	if ( this->edges_to_set[pos] <= 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::removeEdgeToSet(%d,%d,%d)] Trying to lower below 0.\n", size, index, i);
		return;
	}

	this->edges_to_set[pos] = this->edges_to_set[pos] - 1;
}

bool SaturationData::doesSetDominate(int size, int index, int i)
{
	int pos = this->getEdgesToSetPosition(size, index, i);

	if ( pos < 0 )
	{
		SD_DEBUG3(3,"--[SaturationData::doesSetDominate(%d,%d,%d)] A bad set of input.\n", size, index, i);
		return false;
	}

	if ( this->edges_to_set[pos] >= size )
	{
		return true;
	}
	else
	{
		return false;
	}
}

/**
 * completion_sets stores the index of the (r-2)-set which
 * 	completes the given edge. Value is -1 if there is no completion.
 */
int SaturationData::initCompletionSets()
{
	int num_completions = nChooseK(this->maxN, 2);
	int completion_size = num_completions * sizeof(int);

	this->completion_sets = (int*) malloc(completion_size);

	for ( int i = 0; i < num_completions; i++ )
	{
		this->completion_sets[i] = -1;
	}

	return completion_size;
}

void SaturationData::cleanCompletionSets()
{
	if ( this->completion_sets != 0 )
	{
		free(this->completion_sets);
		this->completion_sets = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyCompletionSets()
{
	int num_completions = nChooseK(this->maxN, 2);

	if ( this->bu_completion_sets[this->depth] == 0 )
	{
		this->bu_completion_sets[this->depth] = (int*) malloc(num_completions * sizeof(int));
	}

	for ( int i = 0; i < num_completions; i++ )
	{
		this->bu_completion_sets[this->depth][i] = this->completion_sets[i];
	}
}

bool SaturationData::compareCompletionSets()
{
	int num_diff = 0;

	int num_completions = nChooseK(this->maxN, 2);

	for ( int i = 0; i < num_completions; i++ )
	{
		if ( this->bu_completion_sets[this->depth][i] != this->completion_sets[i] )
		{
			num_diff++;
		}
	}

	if ( num_diff > 0 )
	{
		SD_DEBUG2(1,"--[SaturationData::compareCompletionSets()] %d different at depth %d.\n", num_diff, this->depth);
		return false;
	}

	return true;
}
#endif

int SaturationData::getCompletionSet(int i, int j)
{
	if ( i < 0 || i >= this->n || j < 0 || j >= this->n )
	{
		SD_DEBUG2(3,"--[SaturationData::getCompletionSet(%d,%d)] Inputs out of bounds.\n", i, j);
		return -1;
	}

	return this->getCompletionSet(indexOf(i, j));
}

int SaturationData::getCompletionSet(int index)
{
	if ( index < 0 || index >= nChooseK(this->n, 2) )
	{
		SD_DEBUG1(3,"--[SaturationData::getCompletionSet(%d)] Inputs out of bounds.\n", index);
		return -1;
	}

	return this->completion_sets[index];
}

bool SaturationData::setCompletionSet(int i, int j, int completionIndex)
{
	if ( i < 0 || i >= this->n || j < 0 || j >= this->n )
	{
		SD_DEBUG3(3,"--[SaturationData::setCompletionSet(%d,%d,%d)] Inputs out of bounds.\n", i, j, completionIndex);
		return false;
	}

	int index = indexOf(i, j);

	if ( this->completion_sets[index] >= 0 )
	{
		SD_DEBUG3(5,"--[SaturationData::setCompletionSet(%d,%d,%d)] Entry not -1.\n", i, j, completionIndex);
		return false;
	}

	this->completion_sets[index] = completionIndex;

	return true;
}

void SaturationData::removeCompletionSet(int i, int j, int completionIndex)
{
	if ( i < 0 || i >= this->n || j < 0 || j >= this->n )
	{
		SD_DEBUG3(3,"--[SaturationData::removeCompletionSet(%d,%d,%d)] Inputs out of bounds.\n", i, j, completionIndex);
		return;
	}

	int index = indexOf(i, j);

	if ( this->completion_sets[index] < 0 )
	{
		/* this behavior may occur... we'll see */
		SD_DEBUG3(5,"--[SaturationData::removeCompletionSet(%d,%d,%d)] Entry is -1.\n", i, j, completionIndex);
		return;
	}
	else if ( this->completion_sets[index] != completionIndex )
	{
		/* this behavior occurs during a failed augmetation */
		SD_DEBUG3(5,"--[SaturationData::removeCompletionSet(%d,%d,%d)] Entry is not given completionIndex.\n", i, j, completionIndex);
		return;
	}

	this->completion_sets[index] = -1;
}

/**
 * completion_mult stores the number of times a 1-edge is in
 * 	a completion.
 */
int SaturationData::initCompletionMult()
{
	int completion_mult_size = nChooseK(this->maxN, 2);

	this->completion_mult = (char*) malloc(completion_mult_size);

	bzero(this->completion_mult, completion_mult_size);

	return completion_mult_size;
}

void SaturationData::cleanCompletionMult()
{
	if ( this->completion_mult != 0 )
	{
		free(this->completion_mult);
		this->completion_mult = 0;
	}
}

char SaturationData::getCompletionMult(int index)
{
	if ( index < 0 || index >= nChooseK(this->n, 2) )
	{
		SD_DEBUG1(3,"--[SaturationData::getCompletionMult(%d)] Bad input.\n", index);
		return 0;
	}

	return this->completion_mult[index];
}

bool SaturationData::addCompletionMult(int index)
{
	if ( index < 0 || index >= nChooseK(this->n, 2) )
	{
		SD_DEBUG1(3,"--[SaturationData::addCompletionMult(%d)] Bad input.\n", index);
		return false;
	}
	this->completion_mult[index] = this->completion_mult[index] + 1;
	return true;
}

void SaturationData::removeCompletionMult(int index)
{
	if ( index < 0 || index >= nChooseK(this->n, 2) )
	{
		SD_DEBUG1(3,"--[SaturationData::removeCompletionMult(%d)] Bad input.\n", index);
		return;
	}
	this->completion_mult[index] = this->completion_mult[index] - 1;
}

int SaturationData::initTwoEdges()
{
	this->two_edges = new TreeSet();
	this->zero_edges = new TreeSet();
	return 2 * sizeof(TreeSet);
}

void SaturationData::cleanTwoEdges()
{
	if ( this->two_edges != 0 )
	{
		delete this->two_edges;
		this->two_edges = 0;
	}
	if ( this->zero_edges != 0 )
	{
		delete this->zero_edges;
		this->zero_edges = 0;
	}
}

#ifdef _USE_BACKUPS_
void SaturationData::copyTwoEdges()
{
	this->bu_two_edges[this->depth] = new TreeSet();

	for ( this->two_edges->resetIterator(); this->two_edges->hasNext(); )
	{
		this->bu_two_edges[this->depth]->add(this->two_edges->next());
	}
}

bool SaturationData::compareTwoEdges()
{
	bool result = true;
	if ( this->bu_two_edges[this->depth]->size() != this->two_edges->size() )
	{
		SD_DEBUG2(1,"--[SaturationData::compareTwoEdges()] Copy has size %d while in-place has size %d.\n", this->bu_two_edges[this->depth]->size(),this->two_edges->size());
		result = false;
	}

	for ( this->two_edges->resetIterator(); this->two_edges->hasNext(); )
	{
		int x = this->two_edges->next();

		if ( this->bu_two_edges[this->depth]->contains(x) == 0 )
		{
			SD_DEBUG1(1,"--[SaturationData::compareTwoEdges()] in-place has entry %d but copy does not.\n",x);
			result = false;
		}
	}

	delete this->bu_two_edges[this->depth];
	this->bu_two_edges[this->depth] = 0;

	return result;
}
#endif

bool SaturationData::addTwoEdge(int index)
{
	if ( this->two_edges->contains(index) == 1 )
	{
		SD_DEBUG1(3,"--[SaturationData::addTwoEdge(%d)] Already contains this edge!\n", index);
		return false;
	}

	this->two_edges->add(index);

	return true;
}

void SaturationData::removeTwoEdge(int index)
{
	if ( this->two_edges->contains(index) == 0 )
	{
		SD_DEBUG1(3,"--[SaturationData::removeTwoEdge(%d)] Does not contain this edge!\n", index);
		return;
	}

	this->two_edges->remove(index);
}

/**
 * Constructor
 */
SaturationData::SaturationData(int r, int maxN)
{
	this->r = r;
	this->maxN = maxN;

	/* this will be the size of the initial graph */
	this->n = 2;

	int data_size = 0;

	data_size += this->initSnapshots();
	data_size += this->initAugmentations();
	data_size += this->initLayerGraph();
	data_size += this->initOneGraph();
	data_size += this->initZeroGraph();
	data_size += this->initAdjacencyMatrix();
	this->adjmat[0] = 2;

	data_size += this->initDegreeSequence();
	data_size += this->initAdjacencyList();
	data_size += this->initSetCommonNeighs();
	data_size += this->initEdgesInSet();
	data_size += this->initEdgesToSet();
	data_size += this->initCompletionSets();
	data_size += this->initCompletionMult();
	data_size += this->initTwoEdges();
	this->two_edges->add(0);

	SD_DEBUG3(5,"--[SaturationData::SaturationData(%d,%d)] Allocated %d bytes.\n", r, maxN, data_size);
}

/**
 * Destructor
 */
SaturationData::~SaturationData()
{
	this->cleanSnapshots();
	this->cleanAugmentations();
	this->cleanLayerGraph();
	this->cleanOneGraph();
	this->cleanZeroGraph();
	this->cleanAdjacencyMatrix();
	this->cleanDegreeSequence();
	this->cleanAdjacencyList();
	this->cleanSetCommonNeighs();
	this->cleanEdgesInSet();
	this->cleanEdgesToSet();
	this->cleanCompletionSets();
	this->cleanCompletionMult();
	this->cleanTwoEdges();
}

/**
 * snapshot()
 */
void SaturationData::snapshot()
{
	this->testSnapshotSize();

	this->snapshots[this->depth] = this->edgestack.size();
	this->snapshotorder[this->depth] = this->n;
	this->augmentorder[this->depth] = this->numAugmentations;

#ifdef _USE_BACKUPS_
	/* make copies */
	this->copyAdjacencyMatrix();
	this->copyDegreeSequences();
	this->copyAdjacencyLists();
	this->copyCommonNeighbors();
	this->copyEdgesInSet();
	this->copyEdgesToSet();
	this->copyCompletionSets();
	this->copyTwoEdges();
#endif

	this->depth = this->depth + 1;
}

/**
 * rollback()
 */
void SaturationData::rollback()
{
	/* unroll the edges */
	this->depth = this->depth - 1;

	if ( this->depth < 0 )
	{
		this->depth = 0;
		return;
	}

	int num_removals = this->edgestack.size() - this->snapshots[this->depth];
	edgetuple* tuples = (edgetuple*) malloc(num_removals * sizeof(edgetuple));
	int t = 0;
	while ( this->edgestack.size() > this->snapshots[this->depth] )
	{
		edgetuple tuple = this->edgestack.top();
		this->edgestack.pop();

		tuples[t] = tuple;
		t++;

		this->removeEdge(tuple.i, tuple.j, tuple.index, tuple.to_type, tuple.orig_type, tuple.inCompletion);
	}

	while ( this->numAugmentations > this->augmentorder[this->depth] )
	{
		this->removeAugmentation();
	}

	free(tuples);

	while ( this->n > this->snapshotorder[this->depth] )
	{
		this->removeVertex();
	}

#ifdef _USE_BACKUPS_
	/* make comparisons */
	bool worked = true;
	worked &= this->compareAdjacencyMatrix();
	worked &= this->compareDegreeSequences();
	worked &= this->compareAdjacencyLists();
	worked &= this->compareCommonNeighbors();
	worked &= this->compareEdgesInSet();
	worked &= this->compareEdgesToSet();
	worked &= this->compareCompletionSets();
	worked &= this->compareTwoEdges();

	if ( !worked )
	{
		printf("--[SaturationData::rollback()] The following tuples were removed:\n");

		for ( int i = 0; i < num_removals; i++ )
		{
			printf("\t %d %d %d %d\n", tuples[i].i, tuples[i].j, tuples[i].to_type, tuples[i].orig_type);
		}
	}
#endif
}

/**
 * augmentEdge()
 *
 * Augment the graph using the following edge change.
 *
 * return FALSE if the constraints are violated.
 */
bool SaturationData::augmentEdge(int i, int j, char type, bool inCompletion)
{
	int index = indexOf(i, j);
	char orig_type = this->adjmat[index];

	/* 1: Verify that the type is compatible */
	if ( type == 2 )
	{
		/* If adding a 2-edge, the edge must have been a 3-edge before */
		if ( orig_type == 3 )
		{
			/* this is a simple change to make */
			this->adjmat[index] = 2;
			edgetuple tuple;
			tuple.i = i;
			tuple.j = j;
			tuple.index = index;
			tuple.orig_type = orig_type;
			tuple.to_type = type;
			tuple.inCompletion = false;

			this->edgestack.push(tuple);
			this->two_edges->add(index);

			return true;
		}
		else
		{
			SD_DEBUG3(2,"--[SaturationData::augmentEdge(%d,%d,%d)] Original type is NOT 3!\n", i,j,(int)type);
			return false;
		}
	}

	if ( type == 1 )
	{
		/* orig_type MUST be 2 */
		if ( orig_type != 2 )
		{
			SD_DEBUG3(2,"--[SaturationData::augmentEdge(%d,%d,%d)] Original type is NOT 2!\n", i,j,(int)type);
			return false;
		}

		/* This edge cannot have a completion */
		if ( this->getCompletionSet(i, j) >= 0 )
		{
			SD_DEBUG3(5,"--[SaturationData::augmentEdge(%d,%d,%d)] Trying to augment a 1 when there is a completion!\n", i,j,(int)type);
			return false;
		}
	}

	if ( type == 0 )
	{
		/* orig_type MUST be 2 */
		if ( orig_type != 2 )
		{
			SD_DEBUG3(2,"--[SaturationData::augmentEdge(%d,%d,%d)] Original type is NOT 2!\n", i,j,(int)type);
			return false;
		}

		if ( inCompletion && this->getCompletionSet(i, j) < 0 )
		{
			SD_DEBUG3(5,"--[SaturationData::augmentEdge(%d,%d,%d)] We need a completion for a 0-edge!\n", i,j,(int)type);
			return false;
		}
	}

	bool success = true;

	/* 2: Generate the new edge types as well as test for other data structures. */
	/* These last things are uniform to all edges, as
	 * 	long as everything worked. */
	if ( type < 2 )
	{
		/* add the edge to adjacency matrix */
		this->adjmat[index] = type;

		/* add to adjacency list */
		success &= this->addNeighbor(i, j, type);
		success &= this->addNeighbor(j, i, type);

		if ( success && type == 1 )
		{
			if ( inCompletion )
			{
				this->addCompletionMult(index);
			}
			/* Modify Edges in Sets  */

			/* Special case Order 2 */
			success &= this->addEdgeInSet(2, index);
			if ( 2 == this->r - 2 )
			{
				/* we now have a completion for ALL common neighbors! */
				Set* common = getSetCommonNeighs(2, index);

				int* array = common->toArray();

				for ( int k = 0; k < common->size(); k++ )
				{
					int vk = array[k];

					for ( int l = k + 1; l < common->size(); l++ )
					{
						int vl = array[l];

						/* need to add completion here!  */
						success &= this->setCompletionSet(vk, vl, index);
					}
				}
			}

			/* Orders 3...r */
			for ( int s = 3; s <= this->getMaxEdgesInSetSize(); s++ )
			{
				/* we need to worry about ALL sets, up to maxN */
				int num_with_ij = numSetsWW(this->maxN, s, i, j);

				int* set_with_ij = (int*) malloc(s * sizeof(int));

				/* over all sets with i and j... */
				for ( int sww_index = 0; sww_index < num_with_ij; sww_index++ )
				{
					getSetWW(s, i, j, sww_index, set_with_ij);
					int full_index = indexOfSet(s, set_with_ij);

					success &= this->addEdgeInSet(s, full_index);
					if ( s == this->r && this->isSetFull(s, full_index) )
					{
						/* we have a K_r! */
						SD_DEBUG3(10, "--[SaturationData::augmentEdge(%d,%d,%d)] Found a K_r!\n", i,j,type);
						success = false;
					}

					if ( s == this->r - 2 && this->isSetFull(s, full_index) )
					{
						/* we now have a completion for ALL common neighbors! */
						Set* common = getSetCommonNeighs(s, full_index);

						int* array = common->toArray();

						for ( int k = 0; k < common->size(); k++ )
						{
							int vk = array[k];

							for ( int l = k + 1; l < common->size(); l++ )
							{
								int vl = array[l];

								/* need to add completion here!  */
								success &= this->setCompletionSet(vk, vl, full_index);
							}
						}
					}
				}

				free(set_with_ij);
			}

			/* Modify Set Common Neighbors and Edges To Sets*/
			/* Are we adding any new common neighbors? */
			/* Do it for i, then j */
			for ( int endpoint = 0; endpoint < 2; endpoint++ )
			{
				int ii = (endpoint == 0) ? i : j;
				int jj = (endpoint == 0) ? j : i;

				/* now, for all sizes */
				for ( int s = 2; s <= this->getMaxCommonNeighsSetSize(); s++ )
				{
					/* Iterate over all sets of size s with ii but not jj */
					/* Only need to worry about sets in the current vertex set */
					int num_with_i_wo_j = numSetsWWO(this->maxN, s, ii, jj);
					int* set_with_i_wo_j = (int*) malloc(s * sizeof(int));

					/* over all sets with ii but not jj */
					for ( int sww_index = 0; sww_index < num_with_i_wo_j; sww_index++ )
					{
						getSetWWO(s, ii, jj, sww_index, set_with_i_wo_j);

						int full_index = indexOfSet(s, set_with_i_wo_j);

						/* Increase the number of edges from set to jj */
						success &= this->addEdgeToSet(s, full_index, jj);

						if ( s <= this->getMaxEdgesToSetSize() && this->doesSetDominate(s, full_index, jj) )
						{
							/* we have a NEW common neighbor! */

							/* Check if this set already is a completion */
							if ( s == this->r - 2 && this->isSetFull(s, full_index) )
							{
								/* we now have a completion for ALL common neighbors! */
								Set* common = this->getSetCommonNeighs(s, full_index);

								for ( common->resetIterator(); common->hasNext(); )
								{
									int ci = common->next();

									/* The edge {ci, jj} has a completion in this set */
									success &= this->setCompletionSet(ci, jj, full_index);
								}
							}

							success &= this->addSetCommonNeigh(s, full_index, jj);
						}
					}

					free(set_with_i_wo_j);
				}
			}
		}

		if ( success && type == 0 )
		{
			/* what to do in the zero case? */
			/* Should there be data structures which remove sets with zeroes? */
		}
	}

	edgetuple tuple;
	tuple.i = i;
	tuple.j = j;
	tuple.index = index;
	tuple.orig_type = orig_type;
	tuple.to_type = type;
	tuple.inCompletion = inCompletion;

	this->edgestack.push(tuple);
	this->two_edges->remove(index);

	return success;
}

/**
 * removeEdge()
 *
 * Remove an edge from the graph using the following edge change.
 *
 * return FALSE if the constraints are violated.
 */
void SaturationData::removeEdge(int i, int j, int index, char to_type, char orig_type, bool inCompletion)
{
	if ( to_type == 2 )
	{
		if ( orig_type == 3 )
		{
			this->adjmat[index] = 3;
			this->two_edges->remove(index);
			return;
		}
		else
		{
			SD_DEBUG3(2,"--[SaturationData::removeEdge(%d,%d,.,.,.)] Reversing %d -> 2, but should be 3->2\n",i,j,orig_type);
			return;
		}
	}

	/* 2: Generate the new edge types as well as test for other data structures. */
	/* These last things are uniform to all edges, as
	 * 	long as everything worked. */
	if ( to_type < 2 )
	{
		if ( to_type == 1 )
		{
			if ( inCompletion )
			{
				this->removeCompletionMult(index);
			}

			/* Modify Set Edges */

			/* Special Case: Order 2 */
			if ( this->r - 2 == 2 )
			{
				/* Get the common neighbors */
				Set* common = getSetCommonNeighs(2, index);

				int* array = common->toArray();

				/* Remove ALL completions within this neighborhood */
				for ( int k = 0; k < common->size(); k++ )
				{
					int vk = array[k];

					for ( int l = k + 1; l < common->size(); l++ )
					{
						int vl = array[l];

						/* need to remove completion here!  */
						this->removeCompletionSet(vk, vl, index);
					}
				}
			}
			this->removeEdgeInSet(2, index);

			/* Modify Set Common Neighbors */
			/* Are we adding any new common neighbors? */
			/* Do it for i, then j */
			for ( int endpoint = 0; endpoint < 2; endpoint++ )
			{
				int ii = (endpoint == 0) ? i : j;
				int jj = (endpoint == 0) ? j : i;

				/* now, for all sizes */
				for ( int s = 2; s <= this->getMaxCommonNeighsSetSize(); s++ )
				{
					/* Iterate over all sets of size s with ii but not jj */
					int num_with_i_wo_j = numSetsWWO(this->maxN, s, ii, jj);
					int* set_with_i_wo_j = (int*) malloc(s * sizeof(int));

					for ( int sww_index = 0; sww_index < num_with_i_wo_j; sww_index++ )
					{
						getSetWWO(s, ii, jj, sww_index, set_with_i_wo_j);

						int full_index = indexOfSet(s, set_with_i_wo_j);

						if ( s <= this->getMaxEdgesToSetSize() && this->doesSetDominate(s, full_index, jj) )
						{
							/* we need to remove a common neighbor! */
							this->removeSetCommonNeigh(s, full_index, jj);

							/* Check if this set already is a completion */
							if ( s == this->r - 2 && this->isSetFull(s, full_index) )
							{
								/* we now have a completion for ALL common neighbors! */
								Set* common = this->getSetCommonNeighs(s, full_index);

								for ( common->resetIterator(); common->hasNext(); )
								{
									int ci = common->next();

									/* The edge {ci, jj} had a completion in this set */
									this->removeCompletionSet(ci, jj, full_index);
								}
							}
						}

						/* Decrease the number of edges from set to jj */
						this->removeEdgeToSet(s, full_index, jj);
					}

					free(set_with_i_wo_j);
				}
			}
			for ( int s = 3; s <= this->getMaxEdgesInSetSize(); s++ )
			{
				int num_with_ij = numSetsWW(this->maxN, s, i, j);

				int* set_with_ij = (int*) malloc(s * sizeof(int));

				for ( int sww_index = 0; sww_index < num_with_ij; sww_index++ )
				{
					getSetWW(s, i, j, sww_index, set_with_ij);
					int full_index = indexOfSet(s, set_with_ij);

					if ( s == this->r - 2 && this->isSetFull(s, full_index) )
					{
						/* we now have a completion for ALL common neighbors! */
						Set* common = getSetCommonNeighs(s, full_index);

						int* array = common->toArray();

						for ( int k = 0; k < common->size(); k++ )
						{
							int vk = array[k];

							for ( int l = k + 1; l < common->size(); l++ )
							{
								int vl = array[l];

								/* need to add completion here!  */
								this->removeCompletionSet(vk, vl, full_index);
							}
						}
					}

					this->removeEdgeInSet(s, full_index);
				}

				free(set_with_ij);
			}
		}

		if ( to_type == 0 )
		{
			/* what to do in the zero case? */
			/* Should there be data structures which remove edges with zeroes? */
		}

		/* Remove the edge from the adjacency matrix */
		this->adjmat[index] = orig_type;

		/* Remove from the adjacency list */
		this->removeNeighbor(i, j, to_type);
		this->removeNeighbor(j, i, to_type);
	}

	if ( orig_type == 2 )
	{
		this->two_edges->add(index);
	}
}

/**
 * augmentVertex()
 *
 * Add a new vertex.
 *
 * return FALSE if there are too many vertices.
 */
bool SaturationData::augmentVertex()
{
	if ( this->n >= this->maxN )
	{
		this->addAugmentation(new Augmentation(ADDVERTEX, -1, 0, 0, this->r, this->n));
		return false;
	}

	this->n = this->n + 1;
	bool success = true;

	for ( int i = 0; success && i < this->n - 1; i++ )
	{
		/* turn those 3-edges into 2-edges */
		success &= this->augmentEdge(i, this->n - 1, 2, false);
	}

	this->addVertexToLayerGraph();
	this->addVertexToOneGraph();
	this->addVertexToZeroGraph();

	this->addAugmentation(new Augmentation(ADDVERTEX, this->n, 0, (int*) 0, this->r, this->n - 1));

	return success;
}

/**
 * removeVertex()
 *
 * Remove the nth vertex.
 */
void SaturationData::removeVertex()
{
	/* edge deletions should have been taken care of during the rollback() */
	this->n = this->n - 1;

	this->removeVertexFromLayerGraph();
	this->removeVertexFromOneGraph();
	this->removeVertexFromZeroGraph();
}

/**
 * augmentFill()
 *
 * Augment by performing a fill operation.
 */
bool SaturationData::augmentFill(int index)
{
	if ( index >= 0 && this->numAugmentations > 0 )
	{
		/* TODO: implement augmentFill for a specific orbit! */
		Augmentation* augment = this->augmentations[this->numAugmentations - 1];

		if ( augment->numOrbits < 0 )
		{
			computeUnassignedOrbits(this, augment);
		}

		if ( index >= augment->numOrbits )
		{
			/* BAD */
			this->addAugmentation(new Augmentation(FILL, -1, 0, 0, this->r, this->n));
			return false;
		}

		int i = 0;

		int* edges = (int*) malloc(this->two_edges->size() * sizeof(int));
		int e_base = 0;

		bool success = true;
		while ( success && augment->orbits[index][i] >= 0 )
		{
			int fill_edge = augment->orbits[index][i];
			int vi, vj;
			indexToPair(fill_edge, vi, vj);

			edges[e_base] = fill_edge;
			e_base++;

			/* not in a completion */
			success = this->augmentEdge(vi, vj, 1, false);
			i++;
		}

		this->addAugmentation(new Augmentation(FILL, e_base, 0, edges, this->r, this->n));
		free(edges); /* was deeply copied in Augmentation() */

		return success;
	}

	/* for all 2-edges, try to make them 1-edges until you fail */
	bool success = true;
	if ( this->two_edges->size() == 0 )
	{
		this->addAugmentation(new Augmentation(FILL, -1, 0, 0, this->r, this->n));
		return false;
	}

	if ( this->has2TypeWithCompletion() )
	{
		this->addAugmentation(new Augmentation(FILL, -1, 0, 0, this->r, this->n));
		return false;
	}

	int* edges = (int*) malloc(this->two_edges->size() * sizeof(int));
	int e_base = 0;
	for ( this->two_edges->resetIterator(); success && this->two_edges->hasNext(); )
	{
		int index = this->two_edges->next();
		int i, j;
		indexToPair(index, i, j);

		success &= this->augmentEdge(i, j, 1, false);
		edges[e_base] = index;
		e_base++;
	}

	this->addAugmentation(new Augmentation(FILL, e_base, 0, edges, this->r, this->n));
	free(edges); /* was deeply copied in Augmentation() */

	return success;
}

/**
 * augmentCompletion()
 *
 * Augment by making a 0-edge with a completion.
 */
bool SaturationData::augmentCompletion(int i, int j, int completionIndex)
{
	int index = indexOf(i, j);
	int oldN = this->n;

	/* simple tests that this will work, such as a 2-edge on {i,j} */
	if ( this->adjmat[index] != 2 )
	{
		this->addAugmentation(new Augmentation(NONEDGE, -1, 0, 0, this->r, this->n));
		return false;
	}

	bool success = true;

	int* completionSet = (int*) malloc((this->r - 2) * sizeof(int));
	indexToSet(this->r - 2, completionIndex, completionSet);

	/* add new vertices, if necessary */
	for ( int ii = 0; success && ii < this->r - 2; ii++ )
	{
		while ( success && completionSet[ii] >= this->n )
		{
			if ( completionSet[ii] >= this->maxN )
			{
				success = false;
			}
			else
			{
				success = this->augmentVertex();
			}
		}
	}

	/* Add all 1-edges in this section, if they are 2-edges */

	/* First, the between */
	for ( int k = 0; success && k < this->r - 2; k++ )
	{
		int vk = completionSet[k];

		int ik_index = indexOf(i, vk);
		if ( this->adjmat[ik_index] == 2 )
		{
			success &= this->augmentEdge(i, vk, 1, true);
		}
		else if ( this->adjmat[ik_index] == 0 )
		{
			/* this is bad! */
			success = false;
		}

		if ( success )
		{
			int jk_index = indexOf(j, vk);
			if ( this->adjmat[jk_index] == 2 )
			{
				success &= this->augmentEdge(j, vk, 1, true);
			}
			else if ( this->adjmat[jk_index] == 0 )
			{
				/* this is bad! */
				success = false;
			}
		}

		/* Second, within the completion set */
		for ( int l = 0; success && l < k; l++ )
		{
			int vl = completionSet[l];

			int kl_index = indexOf(vk, vl);

			if ( this->adjmat[kl_index] == 2 )
			{
				success &= this->augmentEdge(vk, vl, 1, true);
			}
			else if ( this->adjmat[kl_index] == 0 )
			{
				/* this is bad! */
				success = false;
			}
		}
	}

	if ( success )
	{
		/* Add the 0-edge between i,j */
		success &= this->augmentEdge(i, j, 0, true);
	}

	this->addAugmentation(new Augmentation(NONEDGE, i, j, completionSet, this->r, oldN));
	free(completionSet); /* deep copied by Augmentation() */

	return success;
}

/**
 * augmentAutoCompletions
 *
 * Augment by performing all auto-completions.
 */
bool SaturationData::augmentAutoCompletions()
{
	bool success = true;

	int nChoose2 = nChooseK(this->n, 2);

	for ( int index = 0; success && index < nChoose2; index++ )
	{
		if ( this->getAdjacency(index) == 2 )
		{
			int cindex = this->getCompletionSet(index);

			if ( cindex >= 0 )
			{
				int i, j;

				indexToPair(index, i, j);
				success = this->augmentCompletion(i, j, cindex);

				if ( !success )
				{
					/* why could this possibly fail? */
					SD_DEBUG3(1,"--[SaturationData::augmentAutoCompletions()] Failure for pair %d, %d, and completion %d.\n", i,j,cindex);
				}
			}
		}
	}

	return success;
}

/**
 * has0TypeWithCompletion
 *
 * Does the current graph have edges of 2-type
 * 	with automatic completions?
 */
bool SaturationData::has2TypeWithCompletion()
{
	/* test if there are any 2-type edges with automatic completions */
	for ( this->two_edges->resetIterator(); this->two_edges->hasNext(); )
	{
		int index = this->two_edges->next();
		int i, j;

		indexToPair(index, i, j);

		if ( this->getCompletionSet(i, j) >= 0 )
		{
			return true;
		}
	}

	return false;
}

/**
 * getLayeredGraph()
 */
sparsegraph* SaturationData::getLayeredGraph()
{
	return this->layer_graph;
}

/**
 * getOneGraph()
 */
sparsegraph* SaturationData::getOneGraph()
{
	return this->one_graph;
}

/**
 * getZeroGraph()
 */
sparsegraph* SaturationData::getZeroGraph()
{
	return this->zero_graph;
}

/**
 * getOneGraphString()
 *
 * NOTE: DO NOT FREE THIS STRING!
 */
char* SaturationData::getOneGraphString()
{
	return sgtos6(this->one_graph);
}

/**
 * isUniquelyKrSaturated
 */
bool SaturationData::isUniquelyKrSaturated()
{
	if ( this->n < 4 )
	{
		return false;
	}

	/* barring all other failures,
	 * 	we only need to see if there are
	 *  any 2-edges. */
	if ( this->two_edges->size() == 0 )
	{
		/* also check full sets */
		/* TODO: remove when done debugging */
		for ( int i = 0; i < nChooseK(this->n, this->r); i++ )
		{
			if ( this->isSetFull(this->r, i) )
			{
				return false;
			}
		}

		return true;
	}

	return false;
}

/**
 * get2Orbits(numOrbits)
 *
 * Returns an array of orbit representatives and
 * 	places the number of orbits in the numOrbits parameter.
 */
int* SaturationData::get2Orbits(int& numOrbits)
{
#ifdef _NO_SYMMETRY_
	/* if not using symmetry... */
	numOrbits = this->two_edges->size();

	if ( this->has2TypeWithCompletion() )
	{
		/* only return those with automatic completions */
		int* reps = (int*) malloc(numOrbits * sizeof(int));

		int k = 0;
		this->two_edges->resetIterator();
		while ( this->two_edges->hasNext() )
		{
			int index = this->two_edges->next();

			int i, j;
			indexToPair(index, i, j);

			if ( this->getCompletionSet(i, j) >= 0 )
			{
				reps[k] = index;
				k++;
			}
		}

		numOrbits = k;

		return reps;
	}

	if ( numOrbits > 0 )
	{
		int* reps = (int*) malloc(numOrbits * sizeof(int));

		this->two_edges->resetIterator();
		for ( int i = 0; i < numOrbits && this->two_edges->hasNext(); i++ )
		{
			reps[i] = this->two_edges->next();
		}

		return reps;
	}

	return 0;
#else

	if ( this->numAugmentations <= 0 )
	{
		numOrbits = this->two_edges->size();

		if ( numOrbits > 0 )
		{
			int* reps = (int*) malloc(numOrbits * sizeof(int));

			this->two_edges->resetIterator();
			for ( int i = 0; i < numOrbits && this->two_edges->hasNext(); i++ )
			{
				reps[i] = this->two_edges->next();
			}

			return reps;
		}
	}

	if ( this->two_edges->size() == 0 )
	{
		numOrbits = 0;
		return 0;
	}

	/* Compute orbits and put them in the most recent Augmentation */
	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	if ( augment->numOrbits < 0 )
	{
		computeUnassignedOrbits(this, augment);
	}

	/* Now, get stuff out of augment */
	numOrbits = augment->numOrbits;
	int* reps = (int*) malloc(numOrbits * sizeof(int));

	bool restrict_to_completion = false;

	if ( 0 )
	{
		restrict_to_completion = this->has2TypeWithCompletion();
	}

	int ri = 0;
	for ( int i = 0; i < numOrbits; i++ )
	{
		int orb_rep = augment->orbits[i][0];

		if ( this->getAdjacency(orb_rep) == 2 )
		{
			if ( !restrict_to_completion || this->getCompletionSet(orb_rep) >= 0 )
			{
				reps[ri] = orb_rep;
				ri++;
			}
		}
	}
	numOrbits = ri;

	return reps;
#endif

}

/**
 * get2OrbitByIndex(int index)
 *
 * Returns all 2-edges which are in	the
 * 	given orbit (by orbit order)
 */
int* SaturationData::get2OrbitByIndex(int index)
{
	if ( this->numAugmentations <= 0 )
	{
		return 0;
	}

	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	return augment->orbits[index];
}

/**
 * get2OrbitIndexByRepresentative(int rep)
 *
 * Returns the orbit index for the 2-edges which are in	the
 * 	given orbit (by an orbit representative)
 */
int SaturationData::get2OrbitIndexByRepresentative(int rep)
{
	if ( this->numAugmentations <= 0 )
	{
		return 0;
	}

	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	for ( int i = 0; i < augment->numOrbits; i++ )
	{
		if ( augment->orbits[i][0] == augment->orbitLabels[rep] )
		{
			return i;
		}
	}

	return 0;
}

/**
 * get2OrbitByRepresentative(int rep)
 *
 * Returns all 2-edges which are in	the
 * 	given orbit (by an orbit representative)
 */
int* SaturationData::get2OrbitByRepresentative(int rep)
{
	if ( this->numAugmentations <= 0 )
	{
		return 0;
	}

	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	for ( int i = 0; i < augment->numOrbits; i++ )
	{
		if ( augment->orbits[i][0] == augment->orbitLabels[rep] )
		{
			return augment->orbits[i];
		}
	}

	return 0;
}

/**
 * getCompletionOrbits(index, numCompleters)
 *
 * Returns an array of completion indices corresponding
 * 	to possible completions of the given index.
 * The number of completion orbits is placed in numCompleters.
 */
int* SaturationData::getCompletionOrbits(int index, int& numCompleters)
{
	/* index is an ORBIT index? */
	if ( index < 0 || index >= nChooseK(this->n, 2) || this->adjmat[index] != 2 )
	{
		SD_DEBUG2(1,"--[SaturationData::getCompletionOrbits(%d,%d)] Bad input!\n", index, numCompleters);

		if ( index >= 0 && index < nChooseK(this->n, 2) )
		{
			SD_DEBUG1(1,"--[SaturationData::getCompletionOrbits(.,.)] Adjmat value is %d\n", this->adjmat[index]);
		}

		return 0;
	}

#ifdef _NO_SYMMETRY_
	int i, j;
	indexToPair(index, i, j);
	int completer = this->getCompletionSet(i, j);

	if ( completer >= 0 )
	{
		numCompleters = 1;
		int* reps = (int*) malloc(2 * sizeof(int));
		reps[0] = completer;
		reps[1] = -1;
		return reps;
	}

	numCompleters = 0;

	for ( int s = 0; s <= this->r - 2; s++ )
	{
		if ( this->n + this->r - 2 - s <= this->maxN )
		{
			/* choose s IN current graph (but not i,j) and r-2-s OUT */
			numCompleters += numSetsWOWO(this->n, s, i, j);
		}
	}

	if ( numCompleters > 0 )
	{
		int* reps = (int*) malloc(numCompleters * sizeof(int));

		int* set = (int*) malloc((this->r - 2) * sizeof(int));

		int rindex = 0;
		for ( int s = 0; s <= this->r - 2; s++ )
		{
			if ( this->n + this->r - 2 - s <= this->maxN )
			{
				for ( int t = s; t < this->r - 2; t++ )
				{
					set[t] = this->n + (t - s);
				}

				for ( int sindex = 0; sindex < numSetsWOWO(this->n, s, i, j); sindex++ )
				{
					getSetWOWO(s, i, j, sindex, set);

					reps[rindex] = indexOfSet(this->r - 2, set);
					rindex++;
				}
			}
		}

		numCompleters = rindex;
		free(set);
		return reps;
	}
#else

	if ( this->numAugmentations <= 0 )
	{
		/* we must be in the case of a single 2-edge */
		numCompleters = 1;

		int* set = (int*) malloc((this->r - 2) * sizeof(int));

		for ( int i = 0; i < this->r - 2; i++ )
		{
			set[i] = i + 2;
		}

		int* reps = (int*) malloc(2 * sizeof(int));
		reps[0] = indexOfSet(this->r - 2, set);
		free(set);
		reps[1] = -1;

		return reps;
	}

	/* Compute orbits and put them in the most recent Augmentation */
	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	if ( augment->numOrbits < 0 )
	{
		computeUnassignedOrbits(this, augment);
	}

	if ( augment->stabilizer != index )
	{
		if ( augment->completionOrbitReps != 0 )
		{
			for ( int i = 0; i < augment->numStabilizedOrbits; i++ )
			{
				if ( augment->completionOrbitReps[i] != 0 )
				{
					free(augment->completionOrbitReps[i]);
					augment->completionOrbitReps[i] = 0;
				}
			}

			free(augment->completionOrbitReps);
			augment->completionOrbitReps = 0;
		}
	}

	int completionIndex = this->getCompletionSet(index);
	if ( completionIndex >= 0 )
	{
		/* there's only one! */
		int* reps = (int*) malloc(2 * sizeof(int));
		reps[0] = completionIndex;
		reps[1] = -1;

		numCompleters = 1;
		return reps;
	}

	computeStabilizedOrbits(this->r, this, index, augment);

	numCompleters = augment->numStabilizedOrbits;
	int* reps = (int*) malloc(numCompleters * sizeof(int));

	int ri = 0;
	for ( int i = 0; i < numCompleters; i++ )
	{
		int orb_rep = indexOfSet(this->r - 2, augment->completionOrbitReps[i]);
		reps[ri] = orb_rep;
		ri++;
	}
	numCompleters = ri;

	return reps;
#endif
}

/**
 * Is the most recent set of augmentations canonical?
 */
bool SaturationData::isCanonical()
{
	if ( this->numAugmentations <= 0 )
	{
		return true;
	}

	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	if ( this->numAugmentations == 1 )
	{
		if ( augment->type != NONEDGE )
		{
			return false;
		}

		/* only one other option */
		return true;
	}

	Augmentation* paugment = this->augmentations[this->numAugmentations - 2];

	if ( paugment->type == FILL )
	{
		return false;
	}

	if ( augment->type != this->canonicalType() )
	{
		return false;
	}

	if ( augment->type != NONEDGE )
	{
		if ( paugment->type == augment->type )
		{
			return false;
		}

		return true;
	}

	/* check all of the 0-edges! */
	Set* possibly_canonical = new TreeSet();

	this->zero_edges->resetIterator();

	while ( this->zero_edges->hasNext() )
	{
		possibly_canonical->add(this->zero_edges->next());
	}

	int aug_index = indexOf(augment->i, augment->j);
	int rank = 0;
	while ( rank <= _CANONICAL_RESTRICT_LAYERCANON && possibly_canonical->contains(aug_index) == 1
			&& possibly_canonical->size() > 1 )
	{
		this->restrictCanonical(possibly_canonical, augment, rank);
		rank++;
	}

	if ( possibly_canonical->size() == 0 )
	{
		/* probably is a disconnected zero graph */
		delete possibly_canonical;
		return false;
	}

	bool isIn = false;
	if ( possibly_canonical->contains(aug_index) == 1 )
	{
		isIn = true;
	}

	delete possibly_canonical;

	return isIn;
}

/**
 * isZeroCanonical()
 *
 * Is the given 2-edge going to be canonical
 * 	in the zero graph if augmented to be zero?
 */
bool SaturationData::isZeroCanonical(int index)
{
	int vi, vj;
	indexToPair(index, vi, vj);
	Augmentation* augment = new Augmentation(NONEDGE, vi, vj, 0, this->r, this->n);

	bool success = this->augmentEdge(vi, vj, 0, false);

	if ( !success )
	{
		this->removeEdge(vi, vj, index, 0, 2, false);
		return false;
	}

	/* compute symmetry for augment */
	computeZeroGraphSymmetry(this, augment);

	/* check all of the 0-edges! */
	Set* possibly_canonical = new TreeSet();

	this->zero_edges->resetIterator();

	while ( this->zero_edges->hasNext() )
	{
		possibly_canonical->add(this->zero_edges->next());
	}

	int aug_index = indexOf(augment->i, augment->j);
	int rank = 0;
	while ( possibly_canonical->size() > 1 && rank <= _CANONICAL_RESTRICT_ZEROCANON )
	{
		this->restrictCanonical(possibly_canonical, augment, rank);
		rank++;
	}

	bool isIn = false;
	if ( possibly_canonical->contains(aug_index) == true )
	{
		isIn = true;
	}

	int minIndex = possibly_canonical->min();

	delete possibly_canonical;
	this->removeEdge(vi, vj, index, 0, 2, false);

	if ( isIn || augment->zeroGraphOrbitLabels[minIndex] == augment->zeroGraphOrbitLabels[aug_index] )
	{
		return true;
	}

	return false;
}

/**
 * restrictCanonical()
 *
 * Given a set of potential 0-edges for deletion, restrict
 * 	the set using the given level of specificity.
 */
void SaturationData::restrictCanonical(Set* potential, Augmentation* augment, int rank)
{
	if ( potential->size() <= 1 )
	{
		return;
	}

	int minLabel1 = this->maxN;
	int minLabel2 = this->maxN;
	int minIndex = this->maxN * this->maxN;

	/* find minimum non-zero tuple index */
	int minTupleIndex = -1;
	int minTupleMult = this->maxN * this->maxN;

	int zeroDegTuple[2];
	int oneDegTuple[2];
	int* zeroDegMult = 0;
	int* oneDegMult = 0;
	bool hasAutoComplete = false;
	Set* auto_set = 0;

	switch ( rank )
	{
#ifdef _DELETE_ALWAYS_CONNECTED_
	case _CANONICAL_RESTRICT_CONNECTIVITY:
	/* Restrict by connectivity */
	if ( this->zero_edges->size() <= 1 )
	{
		return;
	}

	potential->resetIterator();

	while ( potential->hasNext() )
	{
		int index = potential->next();
		int vi, vj;
		indexToPair(index, vi, vj);

		/* test for connectivity */
		if ( this->getDegree(vi, 0) == 1 || this->getDegree(vj, 0) == 1 )
		{
			if ( this->getDegree(vi, 0) == 1 && this->getDegree(vj, 0) == 1 )
			{
				/* this is a disconnected graph! */
				potential->remove(index);
			}
		}
		else
		{
			/* no leaf, so there would be two non-trivial components if it continued */
			/* just need to test getting from vi to vj without the edge {vi,vj} */
			if ( isCutEdge(this->zero_graph, vi, vj) )
			{
				potential->remove(index);
			}
		}
	}

	break;
#endif

#ifdef _DELETE_AUTOMATIC_
	case _CANONICAL_RESTRICT_AUTOMATIC:
		/* check for an automatic completion */

		auto_set = new TreeSet();
		hasAutoComplete = false;
		while ( potential->hasNext() )
		{
			int index = potential->next();
			int vi, vj;
			indexToPair(index, vi, vj);

			int completion = this->getCompletionSet(index);
			int* cset = (int*) malloc((this->r - 2) * sizeof(int));
			indexToSet(this->r - 2, completion, cset);

			bool isAutoComplete = true;
			for ( int k = 0; isAutoComplete && k < this->r - 2; k++ )
			{
				int vk = cset[k];
				int ik_index = indexOf(vi, vk);
				int jk_index = indexOf(vj, vk);

				if ( this->getCompletionMult(ik_index) == 1 )
				{
					isAutoComplete = false;
				}
				else if ( this->getCompletionMult(jk_index) == 1 )
				{
					isAutoComplete = false;
				}
				else
				{
					for ( int l = 0; isAutoComplete && l < k; l++ )
					{
						int vl = cset[l];
						int lk_index = indexOf(vl, vk);

						if ( this->getCompletionMult(lk_index) == 1 )
						{
							isAutoComplete = false;
						}
					}
				}
			}

			if ( isAutoComplete )
			{
				hasAutoComplete = true;
				auto_set->add(index);
			}
		}

		if ( hasAutoComplete )
		{
			potential->clear();

			auto_set->resetIterator();

			while ( auto_set->hasNext() )
			{
				potential->add(auto_set->next());
			}
		}

		delete auto_set;

		break;
#endif

	case _CANONICAL_RESTRICT_ZERO_MAX_DEGREES:
		if ( _DO_DEGREE_CHECK_ )
		{
			int maxDegSum = 0;
			potential->resetIterator();

			while ( potential->hasNext() )
			{
				int index = potential->next();
				int vi, vj;
				indexToPair(index, vi, vj);

				int degSum = this->getDegree(vi, 0) + this->getDegree(vj, 0);

				if ( degSum > maxDegSum )
				{
					maxDegSum = degSum;
				}
			}

			potential->resetIterator();

			while ( potential->hasNext() )
			{
				int index = potential->next();
				int vi, vj;
				indexToPair(index, vi, vj);

				int degSum = this->getDegree(vi, 0) + this->getDegree(vj, 0);

				if ( degSum < maxDegSum )
				{
					potential->remove(index);
				}
			}

		}

		break;

	case _CANONICAL_RESTRICT_ZERO_MULT_DEGREES:

		if ( 1 )
		{
			return;
		}

		/* Restrict by 0-degrees */
		zeroDegMult = (int*) malloc(this->maxN * this->maxN * sizeof(int));
		bzero(zeroDegMult, this->maxN * this->maxN * sizeof(int));

		potential->resetIterator();

		while ( potential->hasNext() )
		{
			int index = potential->next();
			int vi, vj;
			indexToPair(index, vi, vj);

			zeroDegTuple[0] = this->getDegree(vi, 0);
			zeroDegTuple[1] = this->getDegree(vj, 0);
			int zeroDegIndex1 = indexOfTuple(this->maxN, 2, zeroDegTuple);

			zeroDegTuple[0] = this->getDegree(vj, 0);
			zeroDegTuple[1] = this->getDegree(vi, 0);
			int zeroDegIndex2 = indexOfTuple(this->maxN, 2, zeroDegTuple);

			if ( zeroDegIndex1 < zeroDegIndex2 )
			{
				(zeroDegMult[zeroDegIndex1])++;
			}
			else
			{
				(zeroDegMult[zeroDegIndex2])++;
			}
		}

		for ( int zi = 0; zi < this->maxN * this->maxN; zi++ )
		{
			if ( zeroDegMult[zi] > 0 && zeroDegMult[zi] < minTupleMult )
			{
				minTupleIndex = zi;
				minTupleMult = zeroDegMult[zi];
			}
		}

		free(zeroDegMult);
		zeroDegMult = 0;

		potential->resetIterator();

		while ( potential->hasNext() )
		{
			int index = potential->next();
			int vi, vj;
			indexToPair(index, vi, vj);

			zeroDegTuple[0] = this->getDegree(vi, 0);
			zeroDegTuple[1] = this->getDegree(vj, 0);
			int zeroDegIndex1 = indexOfTuple(this->maxN, 2, zeroDegTuple);

			zeroDegTuple[0] = this->getDegree(vj, 0);
			zeroDegTuple[1] = this->getDegree(vi, 0);
			int zeroDegIndex2 = indexOfTuple(this->maxN, 2, zeroDegTuple);

			if ( zeroDegIndex1 < zeroDegIndex2 && zeroDegIndex1 != minTupleIndex )
			{
				potential->remove(index);
			}
			else if ( zeroDegIndex2 != minTupleIndex )
			{
				potential->remove(index);
			}
		}

		break;

	case _CANONICAL_RESTRICT_ZEROCANON:

		if ( _AUGMENT_0_EDGE_CHECK_ )
		{
			/* Restrict by canonical 0-edge */
			if ( augment->numZeroGraphOrbits < 0 )
			{
				computeZeroGraphSymmetry(this, augment);
			}

			potential->resetIterator();
			while ( potential->hasNext() )
			{
				int index = potential->next();
				int i, j;

				indexToPair(index, i, j);

				int label1 = augment->zeroGraphCanonLabels[i];
				int label2 = augment->zeroGraphCanonLabels[j];

				if ( label2 < label1 )
				{
					int temp = label1;
					label1 = label2;
					label2 = temp;
				}

				if ( label1 < minLabel1 )
				{
					minLabel1 = label1;
					minLabel2 = label2;
					minIndex = index;
				}
				else if ( label1 == minLabel1 && label2 < minLabel2 )
				{
					minLabel1 = label1;
					minLabel2 = label2;
					minIndex = index;
				}
			}

			potential->resetIterator();
			while ( potential->hasNext() )
			{
				int index = potential->next();

				if ( augment->zeroGraphOrbitLabels[index] != augment->zeroGraphOrbitLabels[minIndex] )
				{
					/* not in the canonical orbit */
					potential->remove(index);
				}
			}
		}
		break;

	case _CANONICAL_RESTRICT_ONE_MAX_DEGREES:
		if ( _DO_DEGREE_CHECK_ )
		{
			int maxDegSum = 0;
			potential->resetIterator();

			while ( potential->hasNext() )
			{
				int index = potential->next();
				int vi, vj;
				indexToPair(index, vi, vj);

				int degSum = this->getDegree(vi, 1) + this->getDegree(vj, 1);

				if ( degSum > maxDegSum )
				{
					maxDegSum = degSum;
				}
			}

			potential->resetIterator();

			while ( potential->hasNext() )
			{
				int index = potential->next();
				int vi, vj;
				indexToPair(index, vi, vj);

				int degSum = this->getDegree(vi, 1) + this->getDegree(vj, 1);

				if ( degSum < maxDegSum )
				{
					potential->remove(index);
				}
			}

		}

		break;

	case _CANONICAL_RESTRICT_ONE_MULT_DEGREES:
		/* Restrict by 1-degree Conditions */
		break;

	case _CANONICAL_RESTRICT_LAYERCANON:
		/* Restrict by canonical labels */
		if ( augment->numOrbits < 0 )
		{
			computeUnassignedOrbits(this, augment);
		}

		potential->resetIterator();
		while ( potential->hasNext() )
		{
			int index = potential->next();
			int i, j;

			indexToPair(index, i, j);

			int label1 = augment->canonicalLabels[i];
			int label2 = augment->canonicalLabels[j];

			if ( label2 < label1 )
			{
				int temp = label1;
				label1 = label2;
				label2 = temp;
			}

			if ( label1 < minLabel1 )
			{
				minLabel1 = label1;
				minLabel2 = label2;
				minIndex = index;
			}
			else if ( label1 == minLabel1 && label2 < minLabel2 )
			{
				minLabel1 = label1;
				minLabel2 = label2;
				minIndex = index;
			}
		}

		potential->resetIterator();
		while ( potential->hasNext() )
		{
			int index = potential->next();

			if ( augment->zeroOrbitLabels[index] != augment->zeroOrbitLabels[minIndex] )
			{
				potential->remove(index);
			}
		}

		/* by this point, we should only have a single canonical orbit */
		break;

	default:
		break;
	}
}

/**
 * canonicalType()
 *
 * Given the current graph, what kind of augmentation
 * 	did we expect?
 */
AUGMENT_TYPE SaturationData::canonicalType()
{

#ifndef _AUGMENT_ALL_VERTICES_FIRST_
	/* check for any all 2's vertex */
	for ( int i = 0; i < this->n; i++ )
	{
		if ( this->getDegree(i, 0) + this->getDegree(i, 1) == 0 )
		{
			//			SD_DEBUG1(1,"--Canonical Type is ADDVERTEX for vertex %d.\n",i);
			return ADDVERTEX;
		}
	}
#endif

	/* check for a 1-edge with 0 completion multiplicity */
	for ( int i = 0; i < nChooseK(this->n, 2); i++ )
	{
		if ( this->adjmat[i] == 1 && this->completion_mult[i] == 0 )
		{
			return FILL;
		}
	}

	return NONEDGE;
}

/**
 * getN Return the current number of vertices.
 */
int SaturationData::getN()
{
	return this->n;
}

/**
 * getMaxN Return the maximum number of vertices.
 */
int SaturationData::getMaxN()
{
	return this->maxN;
}

/**
 * getAdjacency return the adjacency for this pair
 * index or vertex pair.
 */
char SaturationData::getAdjacency(int index)
{
	if ( index < 0 || index >= nChooseK(this->n, 2) )
	{
		return -1;
	}

	return this->adjmat[index];
}

char SaturationData::getAdjacency(int i, int j)
{
	return this->getAdjacency(indexOf(i, j));
}

/**
 * getNum2Edges
 */
int SaturationData::getNum2Edges()
{
	return this->two_edges->size();
}

/**
 * printEdgeStack
 */
void SaturationData::printEdgeStack()
{
	Augmentation* augment = this->augmentations[this->numAugmentations - 1];

	if ( augment->numOrbits < 0 )
	{
		computeUnassignedOrbits(this, augment);
	}

	/* check all of the 0-edges! */
	Set* possibly_canonical = new TreeSet();

	this->zero_edges->resetIterator();

	while ( this->zero_edges->hasNext() )
	{
		possibly_canonical->add(this->zero_edges->next());
	}

	int aug_index = indexOf(augment->i, augment->j);
	int minLabel1 = this->maxN;
	int minLabel2 = this->maxN;
	int minIndex = -1;

	possibly_canonical->resetIterator();
	while ( possibly_canonical->hasNext() )
	{
		int index = possibly_canonical->next();
		int i, j;

		indexToPair(index, i, j);

		int label1 = augment->canonicalLabels[i];
		int label2 = augment->canonicalLabels[j];

		if ( label2 < label1 )
		{
			int temp = label1;
			label1 = label2;
			label2 = temp;
		}

		if ( label1 < minLabel1 )
		{
			minLabel1 = label1;
			minLabel2 = label2;
			minIndex = index;
		}
		else if ( label1 == minLabel1 && label2 < minLabel2 )
		{
			minLabel1 = label1;
			minLabel2 = label2;
			minIndex = index;
		}
	}

	int i, j;
	indexToPair(minIndex, i, j);
	printf("-- Canonical edge to delete is %d %d.\n", i, j);

	printf("-- Zero Orbit Labels: ");

	possibly_canonical->resetIterator();
	while ( possibly_canonical->hasNext() )
	{
		int index = possibly_canonical->next();
		int lab = augment->zeroOrbitLabels[index];

		indexToPair(index, i, j);
		printf("(%d,%d)[%d] ", i, j, lab);
	}
	printf("\n");

	delete possibly_canonical;

	std::stack<edgetuple> tempstack;

	int sdepth = this->depth - 1;
	while ( this->edgestack.size() > 0 )
	{
		if ( sdepth >= 0 && this->edgestack.size() <= this->snapshots[sdepth] )
		{
			printf("----------------------------------------------------------\n");
			sdepth--;
		}

		edgetuple tuple = this->edgestack.top();
		this->edgestack.pop();
		tempstack.push(tuple);

		if ( tuple.orig_type != 3 )
		{
			printf("i %2d\tj %2d\tindex %2d\torig_type %2d\tto_type %2d\n", tuple.i, tuple.j, tuple.index,
					tuple.orig_type, tuple.to_type);
		}
	}

	while ( tempstack.size() > 0 )
	{
		this->edgestack.push(tempstack.top());
		tempstack.pop();
	}

}
