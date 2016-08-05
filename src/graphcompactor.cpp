/*
 * graphcompactor.c
 *
 *  Created on: Apr 18, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include "nausparse.h"
#include "gtools.h"
#include "graphcompactor.h"

/**
 * compactsg
 *
 * Given a sparsegraph which may have gaps in its representation (the e array),
 * 	compact it to have no gaps.
 */
sparsegraph* compactsg( sparsegraph* g )
{
	sparsegraph* cg = (sparsegraph*) malloc(sizeof(sparsegraph));

	SG_INIT(*cg);

	cg->nv = g->nv;
	cg->vlen = g->nv;
	cg->dlen = g->nv;

	cg->v = (size_t*) malloc(cg->vlen * sizeof(size_t));
	cg->d = (int*) malloc(cg->dlen * sizeof(int));

	int i, j;
	int vindex = 0;
	for ( i = 0; i < g->nv; i++ )
	{
		cg->v[i] = vindex;
		cg->d[i] = g->d[i];
		vindex += g->d[i];
	}

	cg->nde = vindex;
	cg->elen = vindex;
	cg->e = (int*) malloc(cg->elen * sizeof(int));

	for ( i = 0; i < g->nv; i++ )
	{
		vindex = cg->v[i];
		int gvindex = g->v[i];
		int d = g->d[i];

		for ( j = 0; j < d; j++ )
		{
			cg->e[vindex + j] = g->e[gvindex + j];
		}
	}

	return cg;
}
