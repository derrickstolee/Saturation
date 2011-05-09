/*
 * SaturationGraph.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: stolee
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "nausparse.h"
#include "gtools.h"

#include "Augmentation.hpp"
#include "SaturationGraph.hpp"
#include "Set.hpp"
#include "translation.hpp"
#include "TreeSet.hpp"
#include "SaturationSymmetry.hpp"

/**
 * regenerateOpenEdges
 *
 * THIS IS A TEMPORARY METHOD
 */
void SaturationGraph::regenerateOpenEdges()
{
	clock_t start_c = clock();
	//	/* Regenerate the 2-type Edge Set */

	int Nchoose2 = nChooseK(this->n, 2);
	if ( 0 )
	{
		this->openedges->clear();

		for ( int i = 0; i < Nchoose2; i++ )
		{
			if ( this->adjmat[i] == 2 )
			{
				this->openedges->add(i);
			}
		}

		/* Double-Check every 0-type edge */

		for ( int i = 0; i < Nchoose2; i++ )
		{
			if ( this->adjmat[i] == 0 )
			{
				if ( this->completions[i] == 0 )
				{
					printf("--[SaturationGraph] There is a 0-type edge with no completion: %d\n", i);
				}
				else
				{
					/* be sure all completion edges are 1-type */
					int vi, vj;
					this->indexToPair(i, vi, vj);

					for ( int k = 0; k < this->r - 2; k++ )
					{
						int vk = this->completions[i][k];
						int ik_index = this->indexOf(vi, vk);
						int jk_index = this->indexOf(vj, vk);

						if ( this->adjmat[ik_index] != 1 )
						{
							printf(
							       "\t\t\t--[SaturationGraph] What should be a 1-type edge is of %d-type: index %d vi=%d vk=%d\n",
							       this->adjmat[ik_index], ik_index, vi, vk);
						}
						if ( this->adjmat[jk_index] != 1 )
						{
							printf(
							       "\t\t\t--[SaturationGraph] What should be a 1-type edge is of %d-type: index %d vj=%d vk=%d\n",
							       this->adjmat[jk_index], jk_index, vj, vk);
						}

						for ( int l = 0; l < k; l++ )
						{
							int vl = this->completions[i][l];

							int il_index = this->indexOf(vi, vl);
							int jl_index = this->indexOf(vj, vl);
							int kl_index = this->indexOf(vl, vk);

							if ( this->adjmat[il_index] != 1 )
							{
								printf(
								       "\t\t\t--[SaturationGraph] What should be a 1-type edge is of %d-type: index %d vi=%d vl=%d\n",
								       this->adjmat[il_index], il_index, vi, vl);
							}
							if ( this->adjmat[jl_index] != 1 )
							{
								printf(
								       "\t\t\t--[SaturationGraph] What should be a 1-type edge is of %d-type: index %d vj=%d vl=%d\n",
								       this->adjmat[jl_index], jl_index, vj, vl);
							}
							if ( this->adjmat[kl_index] != 1 )
							{
								printf(
								       "\t\t\t--[SaturationGraph] What should be a 1-type edge is of %d-type: index %d vk=%d vl=%d\n",
								       this->adjmat[kl_index], kl_index, vk, vl);
							}
						}
					}
				}
			}
		}
	}

	if ( 0 )
	{
		/* NOW: Check that all 1-edges have completemult > 0 */
		for ( int i = 0; i < Nchoose2; i++ )
		{
			if ( this->adjmat[i] == 1 && this->completemult[i] <= 0 )
			{
				printf("\t\t\t--[SaturtionGraph] There is a 1-type edge NOT in any completion! %d\n", i);
			}

		}
	}

	/* fix degree sequences */
	//	bzero(this->zeroDegrees, this->N * sizeof(int));


	for ( int i = 0; i < this->n; i++ )
	{
		int ionedeg = 0;
		int izerodeg = 0;
		for ( int j = 0; j < this->n; j++ )
		{
			if ( j != i )
			{
				int index = indexOf(i, j);

				if ( this->adjmat[index] == 0 )
				{
					izerodeg++;
				}
				else if ( this->adjmat[index] == 1 )
				{
					ionedeg++;
				}
			}
		}

		if ( izerodeg != this->zeroDegrees[i] )
		{
			this->zeroDegrees[i] = izerodeg;
		}

		if ( ionedeg != this->oneDegrees[i] )
		{
			this->oneDegrees[i] = ionedeg;
		}
	}

	clock_t end_c = clock();
	(this->time_in_regenerate) = this->time_in_regenerate + (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
}

/**
 * compactTheGraph() Compact the graph g into small_g.
 */
void SaturationGraph::compactTheGraph()
{
	this->regenerateOpenEdges();

	if ( this->small_g != 0 )
	{
		/* we need to free small_g */
		SG_FREE((*(this->small_g)));
		free(this->small_g);
		this->small_g = 0;
	}

	clock_t start_c = clock();
	this->small_g = (sparsegraph*) malloc(sizeof(sparsegraph));

	SG_INIT((*(this->small_g)));

	/* small_g has two levels: one for each type (0/1) of edge */
	int sn = 2 * this->n;
	this->small_g->nv = sn;

	this->small_g->vlen = sn;
	this->small_g->dlen = sn;

	this->small_g->v = (int*) malloc(sn * sizeof(int));
	this->small_g->d = (int*) malloc(sn * sizeof(int));

	int vindex = 0;
	for ( int i = 0; i < this->n; i++ )
	{
		this->small_g->v[i] = vindex;
		this->small_g->d[i] = 1 + this->zeroDegrees[i];
		vindex += 1 + this->zeroDegrees[i];
	}
	for ( int i = 0; i < this->n; i++ )
	{
		this->small_g->v[this->n + i] = vindex;
		this->small_g->d[this->n + i] = 1 + this->oneDegrees[i];
		vindex += 1 + this->oneDegrees[i];
	}

	int sde = vindex;
	this->small_g->nde = sde;
	this->small_g->elen = sde;
	this->small_g->e = (int*) malloc(sde * sizeof(int));

	for ( int i = 0; i < this->n; i++ )
	{
		vindex = this->small_g->v[i];
		int vdeg = this->small_g->d[i];

		/* the cross-bar */
		this->small_g->e[vindex] = i + this->n;

		int ve = 1;
		for ( int j = 0; j < this->n; j++ )
		{
			if ( j != i )
			{
				int index = this->indexOf(i, j);

				/* 0-type edges in this layer */
				if ( this->adjmat[index] == 0 )
				{
					this->small_g->e[vindex + ve] = j;
					ve++;
				}
			}
		}
	}
	for ( int i = 0; i < this->n; i++ )
	{
		vindex = this->small_g->v[this->n + i];
		int vdeg = this->small_g->d[this->n + i];

		/* the cross-bar */
		this->small_g->e[vindex] = i;

		int ve = 1;
		for ( int j = 0; j < this->n; j++ )
		{
			if ( j != i )
			{
				int index = this->indexOf(i, j);

				/* 1-type edges in this layer */
				if ( this->adjmat[index] == 1 )
				{
					this->small_g->e[vindex + ve] = this->n + j;
					ve++;
				}
			}
		}
	}

	this->g_updated = false;

	clock_t end_c = clock();
	(this->time_in_compact) = this->time_in_compact + (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
}

/**
 * Initialize all data elements.
 *
 * @param r the saturation parameter.
 * @param n the order of the initial graph.
 */
void SaturationGraph::initData( int r, int N )
{
	this->r = r;
	this->N = N;

	int Nchoose2 = (this->N * (this->N - 1)) / 2;
	this->adjmat = (char*) malloc(Nchoose2);
	this->completemult = (char*) malloc(Nchoose2);
	this->completions = (int**) malloc(Nchoose2 * sizeof(int*));
	this->zeroDegrees = (int*) malloc(this->N * sizeof(int));
	this->oneDegrees = (int*) malloc(this->N * sizeof(int));

	this->openedges = new TreeSet();
	for ( int i = 0; i < Nchoose2; i++ )
	{
		this->openedges->add(i);

		/* 2-edges fill the graph */
		this->adjmat[i] = 2;

		/* multiplicities are zero */
		this->completemult[i] = 0;

		/* completions are null */
		this->completions[i] = 0;
	}

	bzero(this->zeroDegrees, this->N * sizeof(int));
	bzero(this->oneDegrees, this->N * sizeof(int));

	/* Now, we fill in the base graph, a K_r minus an edge, plus an isolated vertex */
	this->n = r + 1;

	this->convertEdge(0, 0);
	this->completions[0] = (int*) malloc((r - 2) * sizeof(int));

	for ( int i = 0; i < r - 2; i++ )
	{
		this->completions[0][i] = i + 2;
	}

	/* fill in rest of K_r with 1-edges */
	for ( int i = 1; i < nChooseK(r, 2); i++ )
	{
		this->convertEdge(i, 1);
	}

	this->g = 0;

	/* this small_g is null until fixed by compactTheGraph() */
	this->small_g = 0;
}

/**
 * convertEdge
 *
 * Take an edge index, and turn it into a new type.
 */
void SaturationGraph::convertEdge( int index, char type )
{
	int i, j;
	this->indexToPair(index, i, j);
	char cur_type = this->adjmat[index];

	if ( index < 0 )
	{
		printf("--[convertEdge] Bad input %d\n", index);
		return;
	}

	if ( index >= nChooseK(this->n, 2) )
	{
		printf("--[convertEdge] Index too large %d\n", index);
		return;
	}

	if ( cur_type == 0 )
	{
		/* can't switch to 1 */
		if ( type == 2 )
		{
			this->adjmat[index] = 2;
			this->openedges->add(index);

			if ( this->completions[index] != 0 )
			{
				/* Free this, now that Augmentation copies completions */
				free(this->completions[index]);
				this->completions[index] = 0;
			}
			else
			{
				printf("--[convertEdge] Trying to delete a completion which is null...(index %d -- %d, %d)\n", index,
				       i, j);
			}

			this->zeroDegrees[i] = this->zeroDegrees[i] - 1;
			this->zeroDegrees[j] = this->zeroDegrees[j] - 1;
		}
		else if ( type == 1 )
		{
			printf("--[convertEdge] trying to convert %d,%d from 0 to 1.\n", i, j);
		}
	}
	else if ( cur_type == 1 )
	{
		/* can't switch to 0 */
		if ( type == 2 )
		{
			int cur_mult = this->completemult[index];

			cur_mult--;
			this->completemult[index] = cur_mult;

			/* only REALLY switch if the count is right */
			if ( cur_mult == 0 )
			{
				this->completemult[index] = 0;
				this->adjmat[index] = 2;
				this->openedges->add(index);
				this->oneDegrees[i] = this->oneDegrees[i] - 1;
				this->oneDegrees[j] = this->oneDegrees[j] - 1;
			}
			else if ( cur_mult < 0 )
			{
				this->completemult[index] = 0;
				this->adjmat[index] = 2;
				this->openedges->add(index);
				printf("--[convertEdge] cur_mult < 0 for index %d (%d, %d)\n", index, i, j);
			}
		}
		else if ( type == 1 )
		{
			/* if ALREADY a 1, just increase multiplicity */
			int cur_mult = this->completemult[index];
			cur_mult++;
			this->completemult[index] = cur_mult;
		}
		else if ( type == 0 )
		{
			printf("--[convertEdge] trying to convert %d,%d from 1 to 0.\n", i, j);
		}
	}
	else if ( cur_type == 2 )
	{
		if ( type != 2 )
		{
			this->openedges->remove(index);
		}

		this->adjmat[index] = type;

		if ( type == 0 )
		{
			this->zeroDegrees[i] = this->zeroDegrees[i] + 1;
			this->zeroDegrees[j] = this->zeroDegrees[j] + 1;
		}
		else if ( type == 1 )
		{
			this->completemult[index] = this->completemult[index] + 1;
			this->oneDegrees[i] = this->oneDegrees[i] + 1;
			this->oneDegrees[j] = this->oneDegrees[j] + 1;
		}
	}
}

/**
 * Initialize a SaturationGraph with saturation parameter r.
 *
 * As a base element, there exists a non-edge with a completion, using r vertices,
 * 	then adds the (r+1)th vertex with all incident edges of type 2.
 */
SaturationGraph::SaturationGraph( int r, int N )
{
	this->initData(r, N);

	this->time_in_regenerate = 0;
	this->time_in_orbits = 0;
	this->time_in_compact = 0;
	this->time_in_stabilized = 0;
	this->time_in_feasible = 0;

	/* change to FALSE to not do implications */
	this->do_zero_implications = false;
}

/**
 * Destructor
 */
SaturationGraph::~SaturationGraph()
{
	printf("T SUM TIME_IN_ORBITS %lf\n", this->time_in_orbits);
	printf("T SUM TIME_IN_COMPACT %lf\n", this->time_in_compact);
	printf("T SUM TIME_IN_REGEN %lf\n", this->time_in_regenerate);
	printf("T SUM TIME_IN_FEASIBLE %lf\n", this->time_in_feasible);
	printf("T SUM TIME_IN_STABILIZE %lf\n", this->time_in_stabilized);

	int Nchoose2 = (this->N * (this->N - 1)) / 2;

	if ( this->adjmat != 0 )
	{
		free(this->adjmat);
		this->adjmat = 0;
	}

	if ( this->completemult != 0 )
	{
		free(this->completemult);
		this->completemult = 0;
	}

	if ( this->completions != 0 )
	{
		for ( int i = 0; i < Nchoose2; i++ )
		{
			if ( this->completions[i] != 0 )
			{
				/* The Augmentation objects COPY these completions */
				free(this->completions[i]);
				this->completions[i] = 0;
			}
		}

		free(this->completions);
		this->completions = 0;
	}

	if ( this->zeroDegrees != 0 )
	{
		free(this->zeroDegrees);
		this->zeroDegrees = 0;
	}

	if ( this->oneDegrees != 0 )
	{
		free(this->oneDegrees);
		this->oneDegrees = 0;
	}

	while ( this->augmentations.size() > 0 )
	{
		Augmentation* augment = this->augmentations.top();
		this->augmentations.pop();

		delete augment;
		augment = 0;
	}

	if ( this->g != 0 )
	{
		SG_FREE((*(this->g)));
		free(this->g);
		this->g = 0;
	}

	if ( this->small_g != 0 )
	{
		SG_FREE((*(this->small_g)));
		free(this->small_g);
		this->small_g = 0;
	}

	if ( this->openedges != 0 )
	{
		delete this->openedges;
	}
}

/****************** AUGMENTATION METHODS ***************************/

/**
 * augment -- Augment by adding a non-edge with completion.
 *
 * Results in modifying the adjacency matrix and pushes information for this augmentation
 * onto the stack.
 *
 * @param noni One endpoint of the non-edge.
 * @param nonj The other endpoint of the non-edge.
 * @param completion An array of r-2 vertices which form the completion.
 *
 * @return True iff the augmentation succeeds without adding a K_r
 * 	or multiple K_r's when adding some assigned non-edge.
 */
bool SaturationGraph::augment( int noni, int nonj, int* completion )
{
	int oldN = this->n;

	this->regenerateOpenEdges();

	if ( this->augmentations.size() > 0 )
	{
		Augmentation* augment = this->augmentations.top();

		if ( augment->type == FILL )
		{
			this->augmentation_succeeded = false;
			/* FILL is the LAST thing to happen! */
			this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
			return false;
		}
	}

	/* TODO: implement in this augmentation */
	int ij_index = this->indexOf(noni, nonj);

	if ( this->adjmat[ij_index] != 2 )
	{
		printf("--[augment] edge {%d %d} is not a 2! it's a %d\n", noni, nonj, this->adjmat[ij_index]);

		if ( this->openedges->contains(ij_index) == 1 )
		{
			printf("\t\t--[augment] but the index %d IS in openEdges.\n", ij_index);
		}
		else
		{
			printf("\t\t--[augment] and the index %d IS NOT in openEdges.\n", ij_index);
		}

		printf("\t-- completion is ");
		for ( int k = 0; k < this->r - 2; k++ )
		{
			printf("%d ", completion[k]);
		}
		printf("\n");

		this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
		char* str = this->getString();
		printf("\t\t\t\t%s", str);
		free(str);
		this->augmentation_succeeded = false;

		return false;
	}

	if ( noni >= this->n || nonj >= this->n )
	{
		printf("--[augment] Trying to add a vertex which is too large: i=%d j=%d n=%d\n", noni, nonj, this->n);
		this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
		this->augmentation_succeeded = false;
		return false;
	}

	/* STEP 1: Check if the completion requires new vertices. */
	/* 1.a. check that these are within range */
	for ( int k = 0; k < this->r - 2; k++ )
	{
		int vk = completion[k];

		if ( vk >= this->N )
		{
			/* too large! */
			this->augmentation_succeeded = false;
			this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
			printf("\t\t--[augment] Trying to augment with an endpoint which is too large.\n");
			return false;
		}
	}

	/* 1.b. Add those new vertices */
	for ( int k = 0; k < this->r - 2; k++ )
	{
		int vk = completion[k];

		while ( vk >= this->n && this->n < this->N )
		{
			//			printf("\t--[augment] Adding a new vertex %d.\n", this->n);
			/* fill new edges with 2's */
			for ( int i = 0; i < this->n; i++ )
			{
				int iton = this->indexOf(i, this->n);

				this->adjmat[iton] = 2;
				this->completemult[iton] = 0;
				this->completions[iton] = 0;
				this->openedges->add(iton);
			}

			this->zeroDegrees[this->n] = 0;
			this->oneDegrees[this->n] = 0;

			this->n = this->n + 1;
		}
	}

	/* STEP 2: Add edges between {noni,nonj} and the completion vertices. */

	/* 2.a: Verify that we can do this augmentation */
	bool working = true;
	for ( int k = 0; working && k < this->r - 2; k++ )
	{
		int vk = completion[k];

		int ki_index = this->indexOf(vk, noni);
		int kj_index = this->indexOf(vk, nonj);

		if ( this->adjmat[ki_index] == 0 )
		{
			working = false;
		}
		else if ( this->adjmat[kj_index] == 0 )
		{
			working = false;
		}
		else
		{
			/* other edges in the completion */
			for ( int l = k + 1; working && l < this->r - 2; l++ )
			{
				int vl = completion[l];
				int kl_index = this->indexOf(vk, vl);

				if ( this->adjmat[kl_index] == 0 )
				{
					working = false;
				}
			} /* end for all l */
		}
	}

	if ( !working )
	{
		this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
		printf("\t\t--[augment] Trying to augment, but there is a 0-edge somewhere!\n");
		this->n = oldN;
		this->augmentation_succeeded = false;
		return false;
	}

	/* from this point on, we may modify the graph */
	/* instead of just returning false, we'll need to  */
	/* fix what we've done */
	for ( int k = 0; k < this->r - 2; k++ )
	{
		int vk = completion[k];

		int ki_index = this->indexOf(vk, noni);
		int kj_index = this->indexOf(vk, nonj);

		this->convertEdge(ki_index, 1);
		this->convertEdge(kj_index, 1);

		for ( int l = k + 1; working && l < this->r - 2; l++ )
		{
			int vl = completion[l];
			int kl_index = this->indexOf(vk, vl);

			this->convertEdge(kl_index, 1);
		} /* end for all l */
	} /* end for all k */

	this->convertEdge(ij_index, 0);
	this->completions[ij_index] = completion;

	Augmentation* augment = new Augmentation(NONEDGE, noni, nonj, completion, this->r, oldN);

	if ( this->do_zero_implications )
	{
		/* what implied zeroes are there? */
		augment->zeroEdges = (int*) malloc(this->n * this->n * sizeof(int));
		bzero(augment->zeroEdges, this->n * this->n * sizeof(int));
		augment->numZeroEdges = 0;

		bool no_doubles = true;

		for ( int index = 0; no_doubles && index < nChooseK(n, 2); index++ )
		{
			if ( index != ij_index && this->adjmat[index] == 2 )
			{
				/* test for implication */
				int i, j;
				indexToPair(index, i, j);

				int ij_completions = 0;

				if ( this->r == 4 )
				{
					for ( int k = 0; no_doubles && k < this->n - 1; k++ )
					{
						if ( k == i || k == j )
						{
							continue;
						}

						int ik_index = this->indexOf(i, k);
						int jk_index = this->indexOf(j, k);

						if ( this->adjmat[ik_index] != 1 || this->adjmat[jk_index] != 1 )
						{
							continue;
						}

						for ( int l = k + 1; no_doubles && l < this->n; l++ )
						{
							if ( l == i || l == j )
							{
								continue;
							}

							int il_index = this->indexOf(i, l);
							int jl_index = this->indexOf(j, l);
							int kl_index = this->indexOf(k, l);

							if ( this->adjmat[il_index] == 1 && this->adjmat[jl_index] == 1 && this->adjmat[kl_index]
							        == 1 )
							{
								ij_completions++;

								if ( ij_completions >= 2 )
								{
									no_doubles = false;
								}
								else
								{
									/* add this edge to the implied set */
									augment->zeroEdges[augment->numZeroEdges] = index;
									(augment->numZeroEdges)++;
									this->completions[index] = (int*) malloc((this->r - 2) * sizeof(int));
									this->completions[index][0] = k;
									this->completions[index][1] = l;
								}
							}
						}
					}
				}
				else
				{
					printf("--[augment] Trying to imply when r = %d.\n", this->r);
				}
			}
		}

		if ( !no_doubles )
		{
			/* we have a problem, let's UNDO what we did! */
			for ( int k = 0; k < this->r - 2; k++ )
			{
				int vk = completion[k];

				int ki_index = this->indexOf(vk, noni);
				int kj_index = this->indexOf(vk, nonj);

				this->convertEdge(ki_index, 2);
				this->convertEdge(kj_index, 2);

				for ( int l = k + 1; working && l < this->r - 2; l++ )
				{
					int vl = completion[l];
					int kl_index = this->indexOf(vk, vl);

					this->convertEdge(kl_index, 2);
				} /* end for all l */
			} /* end for all k */

			for ( int k = 0; k < augment->numZeroEdges; k++ )
			{
				free(this->completions[augment->zeroEdges[k]]);
				this->completions[augment->zeroEdges[k]] = 0;
			}

			this->convertEdge(ij_index, 2);
			this->completions[ij_index] = 0;

			delete augment;
			this->augmentations.push(new Augmentation(NONEDGE, -1, -1, 0, 0, oldN));
			this->n = oldN;
			this->augmentation_succeeded = false;
			return false;
		}
		else
		{
			for ( int k = 0; k < augment->numZeroEdges; k++ )
			{
				int index = augment->zeroEdges[k];
				this->convertEdge(index, 0);

				int ii, jj;
				indexToPair(index, ii, jj);

				for ( int l = 0; l < this->r - 2; l++ )
				{
					int vl = this->completions[index][l];
					this->convertEdge(indexOf(ii, vl), 1);
					this->convertEdge(indexOf(jj, vl), 1);

					for ( int l2 = l + 1; l2 < this->r - 2; l2++ )
					{
						int vl2 = this->completions[index][l2];
						this->convertEdge(indexOf(vl, vl2), 1);
					}
				}
			}
		}
	}

	this->augmentations.push(augment);
	this->augmentation_succeeded = true;
	this->regenerateOpenEdges();

	return true;
}

/**
 * fill -- Fill all unassigned edges with assigned edges.
 *
 * @return True iff the augmentation succeeds without adding a K_r
 * 	or multiple K_r's when adding some assigned non-edge.
 */
bool SaturationGraph::fill()
{
	this->regenerateOpenEdges();

	/* can't fill twice in a row */
	Augmentation* augment = this->augmentations.top();

	if ( this->openedges->size() == 0 || augment->type == FILL )
	{
		this->augmentation_succeeded = false;
		this->augmentations.push(new Augmentation(FILL, -1, 0, 0, 0, this->n));

		return false;
	}

	int oldN = this->n;

	/* fill all 2-edges to be 1-edges */
	int* completion = (int*) malloc(this->openedges->size() * sizeof(int));
	int i = 0;

	for ( this->openedges->resetIterator(); this->openedges->hasNext(); )
	{
		int index = this->openedges->next();
		int vi, vj;
		this->indexToPair(index, vi, vj);

		/* DON'T use convertEdge, as it increases complete mult */
		if ( this->adjmat[index] == 2 )
		{
			this->adjmat[index] = 1;

			this->oneDegrees[vi] = this->oneDegrees[vi] + 1;
			this->oneDegrees[vj] = this->oneDegrees[vj] + 1;

			completion[i] = index;
			i++;
		}
		else
		{
			printf("\t\t--[SaturationGraph::fill] Trying to convert a %d-edge to a 1-edge (index: %d).\n",
			       this->adjmat[index], index);
		}
	}

	this->openedges->clear();
	this->augmentation_succeeded = true;

	this->augmentations.push(new Augmentation(FILL, i, 0, completion, this->r, oldN));
	/* completion is copied by Augmentation */
	free(completion);
	this->regenerateOpenEdges();

	/* TODO: Halt early (and clean up) if finding a K_r or two K_r-completions */
	return true;
}

/**
 * addVertex -- Add a vertex to the graph with all edges unassigned.
 */
bool SaturationGraph::addVertex()
{
	this->regenerateOpenEdges();

	int oldN = this->n;

	/* can't add a vertex twice in a row */
	if ( this->augmentations.size() > 0 )
	{
		Augmentation* augment = this->augmentations.top();

		if ( augment->type == ADDVERTEX )
		{
			this->augmentation_succeeded = false;
			this->augmentations.push(new Augmentation(ADDVERTEX, -1, 0, 0, 0, oldN));
			return false;
		}
	}

	if ( this->n >= this->N )
	{
		/* we cannot make a larger graph */
		this->augmentations.push(new Augmentation(ADDVERTEX, -1, -1, 0, 0, oldN));

		this->augmentation_succeeded = false;
		return false;
	}

	/* fill new edges with 2's */
	for ( int i = 0; i < oldN; i++ )
	{
		int iton = this->indexOf(i, oldN);

		this->adjmat[iton] = 2;
		this->completemult[iton] = 0;
		this->completions[iton] = 0;
		this->openedges->add(iton);
	}

	this->oneDegrees[oldN] = 0;
	this->zeroDegrees[oldN] = 0;

	Augmentation* augment = new Augmentation(ADDVERTEX, this->n, 0, 0, this->r, oldN);
	this->augmentations.push(augment);

	this->n = this->n + 1;
	this->augmentation_succeeded = true;
	this->regenerateOpenEdges();

	return true;
}

/**
 * pop() -- Remove the top of the stack, reverting to the lower state of the graph.
 *
 * WARNING: This will pop off the top of the augmentations stack, so do not call
 * this if an augmentation failed.
 */
void SaturationGraph::pop()
{
	if ( this->augmentations.size() == 0 )
	{
		return;
	}

	/* remove top augmentation from the stack */
	Augmentation* augment = this->augmentations.top();
	this->augmentations.pop();

	if ( this->augmentation_succeeded == false )
	{
		delete augment;
		this->augmentation_succeeded = true;
		return;
	}

	if ( augment->i < 0 )
	{
		delete augment;
		return;
	}

	/* revert the data structure */
	if ( augment->type == NONEDGE )
	{
		int noni, nonj;
		noni = augment->i;
		nonj = augment->j;

		int nonpair = this->indexOf(noni, nonj);

		/* for all vertices in the completion */
		for ( int k = 0; k < this->r - 2; k++ )
		{
			int vertk = this->completions[nonpair][k];

			/* Decrease the multiplicity of the edges ik and jk */
			int ktoi = this->indexOf(noni, vertk);
			int ktoj = this->indexOf(nonj, vertk);

			this->convertEdge(ktoi, 2);
			this->convertEdge(ktoj, 2);

			for ( int l = k + 1; l < this->r - 2; l++ )
			{
				int vertl = this->completions[nonpair][l];

				/* decrease the multiplicity within the completer */
				int ktol = this->indexOf(vertk, vertl);
				this->convertEdge(ktol, 2);
			}
		}

		/* the array is deleted with augment */
		this->convertEdge(nonpair, 2);

		if ( this->do_zero_implications )
		{
			for ( int a = 0; a < augment->numZeroEdges; a++ )
			{
				int index = augment->zeroEdges[a];

				int ii, jj;
				indexToPair(index, ii, jj);

				for ( int l = 0; l < this->r - 2; l++ )
				{
					int vl = this->completions[index][l];
					this->convertEdge(indexOf(ii, vl), 2);
					this->convertEdge(indexOf(jj, vl), 2);

					for ( int l2 = l + 1; l2 < this->r - 2; l2++ )
					{
						int vl2 = this->completions[index][l2];
						this->convertEdge(indexOf(vl, vl2), 2);
					}
				}

				this->convertEdge(index, 2);
			}
		}

		/* check if this augmentation added vertices */
		for ( int k = this->n - 1; k >= augment->oldN; k-- )
		{
			/* this vertex is now 2-only, remove it */
			for ( int i = 0; i < k; i++ )
			{
				int itok = this->indexOf(i, k);
				this->convertEdge(itok, 2);
				this->openedges->remove(itok);
			}
		}

		this->n = augment->oldN;
	}
	else if ( augment->type == FILL )
	{
		for ( int k = 0; k < augment->i; k++ )
		{
			int pairk = augment->completion[k];
			int pi, pj;
			this->indexToPair(pairk, pi, pj);
			this->adjmat[pairk] = 2;
			this->oneDegrees[pi] = this->oneDegrees[pi] - 1;
			this->oneDegrees[pj] = this->oneDegrees[pj] - 1;
			this->openedges->add(pairk);
		}
	}
	else if ( augment->type == ADDVERTEX )
	{
		/* just remove the last vertex */
		/* the adjmat should already be 2's for incident edges */
		/* which means the completions are all null */
		this->n = this->n - 1;
		for ( int i = 0; i < this->n; i++ )
		{
			/* this edge is no longer in the graph */
			int iton = this->indexOf(i, this->n);

			if ( this->adjmat[iton] == 1 )
			{
				/* ?*/
			}

			if ( this->openedges->contains(iton) == 1 )
			{
				this->openedges->remove(iton);
			}
		}
	}

	this->regenerateOpenEdges();

	/* delete the augmentation */
	delete augment;
	augment = 0;
}

/****************** SYMMETRY METHODS ***************************/

/**
 * computeOrbits
 *
 * Use the current assignment of edges to generate a
 * list of pair orbits.
 */
void SaturationGraph::computeOrbits()
{
	/* get the top augmentation */
	Augmentation* augment = this->augmentations.top();

	if ( augment->numOrbits >= 0 )
	{
		/* orbits are calculated */
		return;
	}

	/* First, build the 2-layer graph small_g */
	this->compactTheGraph();

	if ( 1 )
	{
		/* Now, feed it into the symmetry method */
		clock_t start_c = clock();
		computeUnassignedOrbits(this, augment);

		clock_t end_c = clock();
		(this->time_in_orbits) = (this->time_in_orbits) + (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	}
	else
	{
		/* fill in the orbitLabels */
		augment->orbitLabels = (int*) malloc(n * n * sizeof(int));
		for ( int i = 0; i < n * n; i++ )
		{
			augment->orbitLabels[i] = -1;
		}

		augment->orbits = (int**) malloc(n * n * sizeof(int*));
		bzero(augment->orbits, n * n * sizeof(int*));

		augment->zeroOrbitLabels = (int*) malloc(n * n * sizeof(int));
		for ( int i = 0; i < n * n; i++ )
		{
			augment->zeroOrbitLabels[i] = -1;
		}

		int orb = 0;
		int zeroindex = 0;

		for ( int pair_ij = 0; pair_ij < nChooseK(n, 2); pair_ij++ )
		{
			/* Is it a representative? */
			/* we MOSTLY care about 2-orbits */
			if ( this->getAdjacency(pair_ij) == 2 )
			{
				/* this is the orbit representative */
				augment->orbitLabels[pair_ij] = orb;
				augment->orbits[orb] = (int*) malloc(n * n * sizeof(int));

				/* pre-fill with terminating -1's */
				for ( int k = 0; k < n * n; k++ )
				{
					augment->orbits[orb][k] = -1;
				}

				/* fill orbits with pair from g */
				augment->orbits[orb][0] = pair_ij;

				int k = 1;

				augment->orbits[orb][k] = -1;
				augment->orbits[orb][k + 1] = -1;

				/* increase the pair orbit */
				orb++;
			}
			else if ( this->getAdjacency(pair_ij) == 0 )
			{
				/* but also care about 1-orbits */
				/* a 1-type orbit rep! */
				/* this is the orbit representative */
				augment->zeroOrbitLabels[pair_ij] = zeroindex;

				/* increase the pair orbit */
				zeroindex++;
			}
		}

		/* set numOrbits */
		augment->numOrbits = orb;

		/* compute canonical labels */
		augment->canonicalLabels = (int*) malloc(n * sizeof(int));

		for ( int i = 0; i < n; i++ )
		{
			augment->canonicalLabels[i] = i;
		}
	}
}

/**
 * numOrbits -- Return the number of pair orbits.
 *
 * @return the number of pair orbits.
 */
int SaturationGraph::numOrbits()
{
	if ( this->augmentations.size() == 0 )
	{
		return 3;
	}

	Augmentation* augment = this->augmentations.top();

	this->computeOrbits();

	return augment->numOrbits;
}

/**
 * getOrbit -- Get the (-1)-terminated array containing the
 * 	ith pair orbit.
 *
 * @param i from 0 to numOrbits()-1, the index of the orbit.
 * @return an array of pair indices, terminated by -1.
 */
int* SaturationGraph::getOrbit( int i )
{
	if ( this->augmentations.size() == 0 )
	{
		return 0;
	}

	Augmentation* augment = this->augmentations.top();

	this->computeOrbits();

	if ( i < 0 || i >= augment->numOrbits )
	{
		/* bad index */
		return 0;
	}

	return augment->orbits[i];
}

/**
 * getOrbitIndexForEdge
 *
 * Given a pair {i,j}, return the orbit label for
 * the orbit containing {i,j}.
 *
 * @param i the first endpoint.
 * @param j the second endpoint.
 * @return an index between 0 and numOrbits() which gives {i,j} in that orbit.
 */
int SaturationGraph::getOrbitIndexForEdge( int i, int j )
{
	if ( this->augmentations.size() == 0 )
	{
		return -1;
	}

	Augmentation* augment = this->augmentations.top();

	this->computeOrbits();

	if ( i < 0 || 2 * i >= this->n * (this->n - 1) )
	{
		/* bad index */
		return -1;
	}

	return augment->orbitLabels[i];
}

/**
 * computeStabilizer
 *
 * Compute and store the automorphism group when a single pair is stabilized.
 *
 * @param i one endpoint
 * @param j another endpoint
 */
void SaturationGraph::computeStabilizer( int i, int j )
{
	/* First: check if the current augmentation HAS stabilizer data */
	Augmentation* augment = this->augmentations.top();

	this->computeOrbits();

	int index = this->indexOf(i, j);

	if ( augment->stabilizer == index )
	{
		return;
	}
	else
	{
		if ( 1 )
		{
			clock_t start_c = clock();
			/* we WANT this to be the right way */
			this->compactTheGraph();
			computeStabilizedOrbits(this->r, this, index, augment);

			clock_t end_c = clock();
			(this->time_in_stabilized) = (this->time_in_stabilized) + (double) (end_c - start_c)
			        / (double) CLOCKS_PER_SEC;
		}
		else
		{
			int nChooseSum = 1;
			for ( int s = 1; s <= r - 2; s++ )
			{
				nChooseSum += nChooseK(n, s);
			}

			augment->completionOrbitReps = (int**) malloc(nChooseSum * sizeof(int*));
			bzero(augment->completionOrbitReps, nChooseSum * sizeof(int*));

			int orb = 0;
			int* temp_set = (int*) malloc((r - 2) * sizeof(int));

			int pairi = i;
			int pairj = j;

			for ( int s = 1; s <= r - 2; s++ )
			{
				if ( n + (r - 2 - s) > this->getMaxN() )
				{
					/* there are too many vertices! */
					continue;
				}

				int nChooseS = nChooseK(this->n, s);

				int base_index = index;
				for ( int set_index = 0; set_index < nChooseS; set_index++ )
				{
					/* add the orbit iff it can be a completion for the pair */

					bool can_complete = true;
					indexToSet(s, set_index, temp_set);

					for ( int k = 0; can_complete && k < s; k++ )
					{
						int ki_index = indexOf(pairi, temp_set[k]);
						int kj_index = indexOf(pairj, temp_set[k]);

						if ( this->getAdjacency(ki_index) == 0 || this->getAdjacency(kj_index) == 0 )
						{
							can_complete = false;
						}

						for ( int l = k + 1; can_complete && l < s; l++ )
						{
							int kl_index = indexOf(temp_set[k], temp_set[l]);
							if ( this->getAdjacency(kl_index) == 0 )
							{
								can_complete = false;
							}
						}
					}

					if ( can_complete )
					{
						/* this is a good completion */
						augment->completionOrbitReps[orb] = (int*) malloc((r - 2) * sizeof(int));
						indexToSet(s, set_index, augment->completionOrbitReps[orb]);

						/* fill in the rest */
						for ( int t = s; t < r - 2; t++ )
						{
							augment->completionOrbitReps[orb][t] = n + (t - s);
						}

						/* increase the set orbit */
						orb++;
					}
				}
			}

			free(temp_set);

			/* set numOrbits */
			augment->numStabilizedOrbits = orb;
		}

		augment->stabilizer = index;
	}
}

/**
 * numStablizedOrbits
 *
 * @param i one endpoint
 * @param j another endpoint
 * @return the number of stabilized completion-orbits for the pair {i,j}
 */
int SaturationGraph::numStabilizedOrbits( int i, int j )
{
	int index = this->indexOf(i, j);

	this->computeStabilizer(i, j);

	Augmentation* augment = this->augmentations.top();

	return augment->numStabilizedOrbits;
}

/**
 * getStabilizedCompletion
 *
 * Get a representative completion for the kth completion-orbit
 * 	within the stabilizer of {i,j}
 *
 * @param i one endpoint
 * @param j another endpoint
 * @param k the orbit index
 * @return a representative of the kth stabilized completion-orbit for the pair {i,j}
 * 	which is NOT to be free()'d by the caller.
 */
int* SaturationGraph::getStabilizedCompletion( int i, int j, int k )
{
	if ( k < 0 )
	{
		/* bad input */
		return 0;
	}

	int index = this->indexOf(i, j);

	this->computeStabilizer(i, j);

	Augmentation* augment = this->augmentations.top();

	if ( k >= augment->numStabilizedOrbits )
	{
		/* index too large */
		return 0;
	}

	return augment->completionOrbitReps[k];
}

/**
 * isCanonical
 *
 * Check if the most recent assignment of a non-edge is a canonical deletion, corresponding to
 * a canonical augmentation at that index. It selects a canonical deletion and tests if this
 * edge is in the same orbit.
 *
 * The recent assignment is stored in the augmentation stack.
 *
 * Note: the other augmentations are immediately canonical.
 *
 * @return True if the given augmentation is canonical for this graph.
 */
bool SaturationGraph::isCanonical()
{
	/* check if the CURRENT augmentation is canonical */
	Augmentation* augment = this->augmentations.top();

	if ( this->augmentations.size() >= 2 )
	{
		this->augmentations.pop();

		Augmentation* augment2 = this->augmentations.top();
		this->augmentations.push(augment);

		if ( augment2->type == FILL )
		{
			return false;
		}
	}
	/* Use a domain-reduction strategy, so that we can shortcut the calculation
	 * 	if we ever eliminate the given augmentation
	 */
	this->regenerateOpenEdges();

	bool has_zero_degree = false;

	for ( int i = this->n - 1; !has_zero_degree && i >= 0; i-- )
	{
		if ( this->zeroDegrees[i] == 0 && this->oneDegrees[i] == 0 )
		{
			has_zero_degree = true;
		}
	}

	if ( has_zero_degree )
	{
		/* this MUST be an ADDVERTEX type */
		if ( augment->type == ADDVERTEX )
		{
			return true;
		}

		return false;
	}

	if ( this->openedges->size() == 0 )
	{
		/* check for a completemult == 0  */
		bool has_zero_mult = false;
		int zero_index = 0;
		for ( int i = 0; !has_zero_mult && i < this->n * (this->n - 1) / 2; i++ )
		{
			if ( this->adjmat[i] == 1 && this->completemult[i] == 0 )
			{
				has_zero_mult = true;
				zero_index = i;
			}
		}

		if ( has_zero_mult && augment->type != FILL )
		{
			/* this can ONLY occur when a FILL happens */
			return false;
		}

		return true;
	}

	if ( augment->type != NONEDGE )
	{
		/* If it did not fit the above cases, it should be a NONEDGE augmentation */
		return false;
	}

	int aug_index = indexOf(augment->i, augment->j);

	/* so, we are of nonedge type */
	/* let us reduce the set of things to consider */
	if ( 0 )
	{

		int* aug_degree_tuple = (int*) malloc(4 * sizeof(int));
		aug_degree_tuple[0] = this->zeroDegrees[augment->i];
		aug_degree_tuple[1] = this->oneDegrees[augment->i];
		aug_degree_tuple[2] = this->zeroDegrees[augment->j];
		aug_degree_tuple[3] = this->oneDegrees[augment->j];

		int aug_deg_index = indexOfTuple(this->n, 4, aug_degree_tuple);

		/* try other ORIENTATION of edge */
		aug_degree_tuple[0] = this->zeroDegrees[augment->j];
		aug_degree_tuple[1] = this->oneDegrees[augment->j];
		aug_degree_tuple[2] = this->zeroDegrees[augment->i];
		aug_degree_tuple[3] = this->oneDegrees[augment->i];

		int temp_index = indexOfTuple(this->n, 4, aug_degree_tuple);

		if ( temp_index < aug_deg_index )
		{
			aug_deg_index = temp_index;
		}
		else
		{
			aug_degree_tuple[0] = this->zeroDegrees[augment->i];
			aug_degree_tuple[1] = this->oneDegrees[augment->i];
			aug_degree_tuple[2] = this->zeroDegrees[augment->j];
			aug_degree_tuple[3] = this->oneDegrees[augment->j];
		}

		int* degree_tuple = (int*) malloc(4 * sizeof(int));

		/* compute degrees at each endpoint */
		int npow4 = pow(this->n, 4);
		int* degreeMult = (int*) malloc(npow4 * sizeof(int));
		bzero(degreeMult, npow4 * sizeof(int));

		Set* posIndices = new TreeSet();

		for ( int index = 0; index < nChooseK(this->n, 2); index++ )
		{
			if ( this->adjmat[index] == 0 )
			{
				int ii, ij;
				indexToPair(index, ii, ij);

				degree_tuple[0] = this->zeroDegrees[ii];
				degree_tuple[1] = this->oneDegrees[ii];
				degree_tuple[2] = this->zeroDegrees[ij];
				degree_tuple[3] = this->oneDegrees[ij];

				int degree_index = indexOfTuple(this->n, 4, degree_tuple);

				degree_tuple[0] = this->zeroDegrees[ij];
				degree_tuple[1] = this->oneDegrees[ij];
				degree_tuple[2] = this->zeroDegrees[ii];
				degree_tuple[3] = this->oneDegrees[ii];

				temp_index = indexOfTuple(this->n, 4, degree_tuple);

				if ( temp_index < degree_index )
				{
					degree_index = temp_index;
				}

				posIndices->add(degree_index);

				degreeMult[degree_index] = degreeMult[degree_index] + 1;
			}
		}

		int min_deg_mult = this->n * this->n;
		int min_deg_index = npow4 + 1;
		for ( posIndices->resetIterator(); posIndices->hasNext(); )
		{
			int degree_index = posIndices->next();
			if ( (degreeMult[degree_index] > 0 && degreeMult[degree_index] < min_deg_mult) || (degreeMult[degree_index]
			        == min_deg_mult && degree_index < min_deg_index) )
			{
				min_deg_mult = degreeMult[degree_index];
				min_deg_index = degree_index;
			}
		}
		delete posIndices;

		if ( 0 && min_deg_index != aug_deg_index )
		{
			printf("--NOT CANONICAL since (%d) min_deg_index %d != %d aug_deg_index (%d)\n", min_deg_mult,
			       min_deg_index, aug_deg_index, degreeMult[aug_deg_index]);

			indexToTuple(this->n, 4, min_deg_index, degree_tuple);

			printf("-- aug_degree_tuple %d %d %d %d\n", aug_degree_tuple[0], aug_degree_tuple[1], aug_degree_tuple[2],
			       aug_degree_tuple[3]);
			printf("--     degree_tuple %d %d %d %d\n", degree_tuple[0], degree_tuple[1], degree_tuple[2],
			       degree_tuple[3]);

			free(degreeMult);
			free(degree_tuple);
			free(aug_degree_tuple);
			return false;
		}

		//	if ( 1 )
		//	{
		//		printf("-- IS CANONICAL since (%d) min_deg_index %d == %d aug_deg_index (%d)\n", min_deg_mult, min_deg_index,
		//		       aug_deg_index, degreeMult[aug_deg_index]);
		//
		//		free(degree_tuple);
		//		free(aug_degree_tuple);
		//		return true;
		//	}

		free(degreeMult);
		/* still in the running */

		//	if ( min_deg_mult == 1 )
		//	{
		//		/* decision made AND it's this one! */
		//		//		printf("--IS CANONICAL since min_deg_index %d = %d aug_deg_index AND unique multiplicity!\n", min_deg_index,
		//		//		       aug_deg_index);
		//		free(degree_tuple);
		//		free(aug_degree_tuple);
		//		return true;
		//	}
	}
	else if ( 0 )
	{
		/* compute based on 1-degrees */
		int onedegi = this->oneDegrees[augment->i];
		int onedegj = this->oneDegrees[augment->j];
		int aug_deg_index = 0;
		if ( onedegi < onedegj )
		{
			aug_deg_index = this->n * onedegi + onedegj;
		}
		else
		{
			aug_deg_index = this->n * onedegj + onedegi;
		}

		int* degreeMult = (int*) malloc(this->n * this->n * sizeof(int));
		bzero(degreeMult, this->n * this->n * sizeof(int));

		Set* posIndices = new TreeSet();

		for ( int index = 0; index < nChooseK(this->n, 2); index++ )
		{
			if ( this->adjmat[index] == 0 )
			{
				int ii, jj;
				indexToPair(index, ii, jj);

				int onedegii = this->oneDegrees[ii];
				int onedegjj = this->oneDegrees[jj];
				int deg_index = 0;
				if ( onedegii < onedegjj )
				{
					aug_deg_index = this->n * onedegii + onedegjj;
				}
				else
				{
					aug_deg_index = this->n * onedegjj + onedegii;
				}

				posIndices->add(deg_index);

				degreeMult[deg_index] = degreeMult[deg_index] + 1;
			}
		}

		int min_deg_mult = this->n * this->n;
		int min_deg_index = this->n * this->n + 1;

		for ( posIndices->resetIterator(); posIndices->hasNext(); )
		{
			int deg_index = posIndices->next();

			if ( degreeMult[deg_index] > 0 && degreeMult[deg_index] < min_deg_mult )
			{
				min_deg_index = deg_index;
				min_deg_mult = degreeMult[deg_index];
			}
			else if ( degreeMult[deg_index] == min_deg_mult && deg_index > min_deg_index )
			{
				min_deg_index = deg_index;
			}
		}

		delete posIndices;
		free(degreeMult);

		if ( min_deg_index != aug_deg_index )
		{
			return false;
		}
		else if ( min_deg_mult == 1 )
		{
			return true;
		}
	}

	/* compute canonical labels now */
	this->computeOrbits();

	/* within this set, find a canonical label */
	/* if this label is ever beaten, then not canonical! */
	int aug_orbit = augment->zeroOrbitLabels[aug_index];
	int aug_label = indexOf(augment->canonicalLabels[augment->i], augment->canonicalLabels[augment->j]);
	int min_label = aug_label;
	int min_index = aug_index;

	for ( int index = 0; index < nChooseK(this->n, 2); index++ )
	{
		if ( this->adjmat[index] == 0 )
		{
			int ii, ij;
			indexToPair(index, ii, ij);

			/* try min canon label */
			int canon_index = indexOf(augment->canonicalLabels[ii], augment->canonicalLabels[ij]);

			if ( canon_index < min_label )
			{
				min_label = canon_index;
				min_index = index;
			}
			//			}
		}
	}

	//	free(degree_tuple);
	//	free(aug_degree_tuple);

	if ( augment->zeroOrbitLabels[min_index] != aug_orbit )
	{
		return false;
	}

	return true;
}

/****************** ACCESSOR METHODS ***************************/

/**
 * getN() -- Get the order of the current graph.
 *
 * @return the number of vertices in the current graph.
 */
int SaturationGraph::getN()
{
	return this->n;
}

/**
 * getMaxN() -- Get the largest possible order
 *
 * @return N
 */
int SaturationGraph::getMaxN()
{
	return this->N;
}

/**
 * indexOf -- Get the index of the i-j edge in the adjacency matrix.
 *
 * @param i one of the vertices
 * @param j another vertex
 * @return a position in the adjacency matrix for the i-j edge. If i==j or either is negative, returns -1.
 */
int SaturationGraph::indexOf( int i, int j )
{
	if ( i == j )
	{
		/* bad input */
		return -1;
	}

	if ( j < i )
	{
		/* wrong order */
		return this->indexOf(j, i);
	}

	/* CO-LEX ORDER */

	/* (j choose 2) + i */
	return ((j * (j - 1)) / 2) + i;
}

/**
 * indexToPair -- Get the pair {i,j} corresponding to the given co-lex index.
 *
 * @param index the co-lex index of a pair.
 * @param i the first vertex.
 * @param j the second vertex.
 */
void SaturationGraph::indexToPair( int index, int& i, int& j )
{
	/* Need to find integer j so that (j*(j-1))/2 = (1/2)j^2 - (1/2)j = index */
	/* That is, x^2 - x - 2*index = 0 */
	/* x = (1 + sqrt(1-4(-2*index)))/2 */
	/* Then, round down */
	j = (int) ((1.0 + sqrt(1.0 + 8.0 * (double) index)) * 0.5);

	/* i is the "remainder" */
	i = index - (j * (j - 1)) / 2;
}

/**
 * isUniquelySaturated
 *
 * @return True iff the current graph (minus a 2-dominating vertex)
 * 	is fully assigned and uniquely K_r-saturated.
 */
bool SaturationGraph::isUniquelySaturated()
{
	this->regenerateOpenEdges();

	/* TODO: Test if the current graph IS uniquely K_r-saturated */
	/* It requires all edges to be filled (so no 2-edges) */
	if ( this->openedges->size() > 0 )
	{
		return false;
	}

	/* and requires no K_r's (should be noticed during augmentation) */
	/* and require all 0-edges to have a completion (necessary for no 2-edges) */
	if ( this->isFeasible() == false )
	{
		return false;
	}

	/* every 0-edge has a completion by definition */
	/* the completion is unique by isFeasible() */
	return true;
}

/**
 * isFeasible
 *
 * @return True iff there is no copy of K_r and there is no
 * 		non-edge with more than one completion.
 */
bool SaturationGraph::isFeasible()
{
	/* this is a very stateful way to do it */
	/* but 'isFeasible' should only be called */
	/* after the corresponding augmentation method */
	if ( this->augmentation_succeeded == false )
	{
		return false;
	}

	clock_t start_c = clock();
	if ( r == 4 )
	{
		/* search for a K_4 */
		for ( int i = 0; i < this->n - 3; i++ )
		{
			/* search among the neighbors of i */
			for ( int j = i + 1; j < this->n - 2; j++ )
			{
				int ij_index = this->indexOf(i, j);
				if ( this->adjmat[ij_index] != 1 )
				{
					continue;
				}

				for ( int k = j + 1; k < this->n - 1; k++ )
				{
					int ik_index = this->indexOf(i, k);
					int jk_index = this->indexOf(j, k);

					if ( this->adjmat[ik_index] != 1 || this->adjmat[jk_index] != 1 )
					{
						continue;
					}

					for ( int l = k + 1; l < this->n; l++ )
					{

						int il_index = this->indexOf(i, l);
						int jl_index = this->indexOf(j, l);
						int kl_index = this->indexOf(k, l);

						if ( this->adjmat[il_index] == 1 && this->adjmat[jl_index] == 1 && this->adjmat[kl_index] == 1 )
						{
							clock_t end_c = clock();
							(this->time_in_feasible) = (this->time_in_feasible) + (double) (end_c - start_c)
							        / (double) CLOCKS_PER_SEC;
							return false;
						}
					}
				}
			}
		}

		/* now, search for MULTIPLE completions */
		for ( int i = 0; i < this->n - 1; i++ )
		{
			for ( int j = i + 1; j < this->n; j++ )
			{
				int ij_index = this->indexOf(i, j);
				if ( this->adjmat[ij_index] == 1 )
				{
					/* want ij to be not a 1-type edge */
					continue;
				}

				int ij_completions = 0;

				for ( int k = 0; k < this->n - 1; k++ )
				{
					if ( k == i || k == j )
					{
						continue;
					}

					int ik_index = this->indexOf(i, k);
					int jk_index = this->indexOf(j, k);

					if ( this->adjmat[ik_index] != 1 || this->adjmat[jk_index] != 1 )
					{
						continue;
					}

					for ( int l = k + 1; l < this->n; l++ )
					{
						if ( l == i || l == j )
						{
							continue;
						}

						int il_index = this->indexOf(i, l);
						int jl_index = this->indexOf(j, l);
						int kl_index = this->indexOf(k, l);

						if ( this->adjmat[il_index] == 1 && this->adjmat[jl_index] == 1 && this->adjmat[kl_index] == 1 )
						{
							ij_completions++;
							if ( ij_completions >= 2 )
							{
								clock_t end_c = clock();
								(this->time_in_feasible) = (this->time_in_feasible) + (double) (end_c - start_c)
								        / (double) CLOCKS_PER_SEC;
								return false;
							}
						}
					}
				}
			}
		}
	}

	clock_t end_c = clock();
	(this->time_in_feasible) = (this->time_in_feasible) + (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	return true;
}

/**
 * getString
 *
 * Get a sparse6 representation of the 1-graph.
 *
 * @return a buffer containing the string. It will be free()'d by the caller.
 */
char* SaturationGraph::getString()
{
	sparsegraph sg;
	SG_INIT(sg);

	int sn = this->n;
	sg.nv = sn;

	sg.vlen = sn;
	sg.dlen = sn;

	sg.v = (int*) malloc(sn * sizeof(int));
	sg.d = (int*) malloc(sn * sizeof(int));

	int vindex = 0;

	for ( int i = 0; i < this->n; i++ )
	{
		sg.v[i] = vindex;
		sg.d[i] = this->oneDegrees[i];
		vindex = vindex + this->oneDegrees[i];
	}

	int sde = vindex;
	sg.nde = sde;
	sg.elen = sde;

	sg.e = (int*) malloc(sde * sizeof(int));

	int ve = 0;
	for ( int i = 0; i < this->n; i++ )
	{
		sg.v[i] = ve;

		for ( int j = 0; j < this->n; j++ )
		{
			if ( j != i )
			{
				int index = this->indexOf(i, j);

				/* 1-type edges in this layer */
				if ( this->adjmat[index] == 1 )
				{
					sg.e[ve] = j;
					ve++;
				}
			}
		}
	}

	char* s6 = sgtos6(&sg);

	/* Also, print adjmat and completemult */
	int Nchoose2 = (this->n * (this->n - 1)) / 2;
	char* buff = (char*) malloc(strlen(s6) + 6 * Nchoose2 + 6 * this->n + 5);

	strcpy(buff, s6);

	free(sg.v);
	free(sg.d);
	free(sg.e);

	return buff;
}

/**
 * getDepth
 *
 * Return the number of augmentations.
 */
int SaturationGraph::getDepth()
{
	return this->augmentations.size();
}

/**
 * getSmallG
 */
sparsegraph* SaturationGraph::getSmallG()
{
	return this->small_g;
}

/**
 * getAdjacency
 */
char SaturationGraph::getAdjacency( int index )
{
	return this->adjmat[index];
}
