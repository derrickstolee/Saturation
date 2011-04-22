/*
 * SaturationDataSymmetry.cpp
 *
 *  Created on: Apr 18, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "graphcompactor.h"
#include "SaturationDataSymmetry.hpp"
#include "translation.hpp"
#include "nausparse.h"
#include "nauty.h"
#include "gtools.h"

/**
 * Used for nauty statistics
 */
double nauty_time = 0;
double nauty_layer_time = 0;
double nauty_zero_time = 0;
double nauty_stab_time = 0;
double orbit_layer_time = 0;
double orbit_zero_time = 0;
double orbit_stab_time = 0;
long long int nauty_calls = 0;
long long int nauty_stab_calls = 0;
long long int nauty_layer_calls = 0;
long long int nauty_zero_calls = 0;

/******** PARTITION JOINING METHODS ***********/

int root(int i, int* p1)
{
	int j = p1[i];

	if ( j == i )
	{
		return j;
	}
	else
	{
		p1[i] = root(j, p1);
		return p1[i];
	}
}

/**
 * Take the two partitions and join them, modifying p1.
 */
void joinPartitions2(int* p1, int* p2, int ne)
{
	bool changed = true;
	while ( changed )
	{
		changed = false;

		for ( int i = 0; i < ne; i++ )
		{
			int rooti = root(i, p1);

			if ( p2[i] == i )
			{
				/* look through a 2-part */
				for ( int j = i; j < ne; j++ )
				{
					int rootj = root(j, p1);

					if ( p2[j] == i && rooti != rootj )
					{
						/* merge the i-part and j-parts */
						if ( rooti < rootj )
						{
							p1[rootj] = rooti;
						}
						else if ( rootj < rooti )
						{
							p1[rooti] = rootj;
						}

						changed = true;
					}
				}
			}
		}
	} /* while changed */
}

/******** GENERATOR METHODS ***********/

int gen_n = 0;
int nchoose2 = 0;
int* pair_orbits = 0;
int* new_partition;

void init_generator_data_pairs(int n)
{
	gen_n = n;
	nchoose2 = nChooseK(n, 2);
	pair_orbits = (int*) malloc(nchoose2 * sizeof(int));
	new_partition = (int*) malloc(nchoose2 * sizeof(int));
	for ( int i = 0; i < nchoose2; i++ )
	{
		pair_orbits[i] = i;
	}
}

void use_generator_pairs(int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n)
{

	for ( int i = 0; i < nchoose2; i++ )
	{
		new_partition[i] = nchoose2 + 1;
	}

	for ( int start_pindex = 0; start_pindex < nchoose2; start_pindex++ )
	{
		int i, j;
		indexToPair(start_pindex, i, j);

		if ( new_partition[start_pindex] >= nchoose2 )
		{
			/* loop with generator! */
			new_partition[start_pindex] = start_pindex;

			/* use even indices */
			int ni = perm[i << 1] >> 1;
			int nj = perm[j << 1] >> 1;

			int pindex = indexOf(ni, nj);

			/* wait until ordered pair match */
			while ( ni != i || nj != j )
			{
				new_partition[pindex] = start_pindex;
				ni = perm[ni << 1] >> 1;
				nj = perm[nj << 1] >> 1;
				pindex = indexOf(ni, nj);
			}
		}
	}

	joinPartitions2(pair_orbits, new_partition, nchoose2);
}

void clean_generator_data_pairs()
{
	nchoose2 = 0;
	if ( pair_orbits != 0 )
	{
		free(pair_orbits);
		pair_orbits = 0;
	}

	if ( new_partition != 0 )
	{
		free(new_partition);
	}
}

/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeUnassignedOrbits(SaturationData* satdata, Augmentation* augment)
{
	sparsegraph* g = compactsg(satdata->getLayeredGraph());

	//	printf("%s", sgtos6(g));
	//
	//	printf("d: ");
	//	for ( int i = 0; i < g->nv; i++ )
	//	{
	//		printf("%d ", g->d[i]);
	//	}
	//	printf("\nv: ");
	//	for ( int i = 0; i < g->nv; i++ )
	//	{
	//		printf("%d ", g->v[i]);
	//	}
	//
	//	printf("\ne: ");
	//	for ( int i = 0; i < g->elen; i++ )
	//	{
	//		if ( g->e[i] < 0 )
	//		{
	//			g->e[i] = 0;
	//		}
	//
	//		printf("%d ", g->e[i]);
	//	}
	//	printf("\n");

	/* clean up old augment stuff, if necessary */
	if ( augment->orbits != 0 )
	{
		for ( int i = 0; i < augment->numOrbits; i++ )
		{
			if ( augment->orbits[i] != 0 )
			{
				free(augment->orbits[i]);
				augment->orbits[i] = 0;
			}
		}

		/* TODO: consider realloc? */
		free(augment->orbits);
		augment->orbits = 0;
	}

	if ( augment->orbitLabels != 0 )
	{
		/* TODO: consider realloc? */
		free(augment->orbitLabels);
		augment->orbitLabels = 0;
	}

	if ( augment->zeroOrbitLabels != 0 )
	{
		free(augment->zeroOrbitLabels);
		augment->zeroOrbitLabels = 0;
	}

	if ( augment->canonicalLabels != 0 )
	{
		free(augment->canonicalLabels);
		augment->canonicalLabels = 0;
	}

	int nv = satdata->getN() << 1; /* 2*n */
	int n = satdata->getN(); /* nv / 2 */
	int m = (nv + WORDSIZE - 1) / WORDSIZE;
	nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	/* labels alternate between even and odd */
	for ( int i = 0; i < n; i++ )
	{
		lab[i] = 2 * i;
		lab[n + i] = 2 * i + 1;
		ptn[i] = 1;
		ptn[n + i] = 1;
	}

	/* the partition splits {0,2,..,2n-2},{1,3,...,2n-1} */
	ptn[n - 1] = 0;
	ptn[nv - 1] = 0;

	static DEFAULTOPTIONS_SPARSEGRAPH( options);

	init_generator_data_pairs(n);

	options.defaultptn = FALSE; /* we DO need colors */
	options.getcanon = TRUE; /* gets labels */
	options.digraph = FALSE;
	options.userautomproc = &use_generator_pairs;

	options.invarproc = &adjacencies_sg; /* sparse version */

	statsblk stats; /* we'll use this at the end */
	DYNALLSTAT(setword, workspace, worksize);
	DYNALLOC1(setword, workspace, worksize, 50 * m, "malloc");

	sparsegraph canon_g;
	SG_INIT(canon_g);

	/* call nauty */
	clock_t start_c = clock();
	time_t start_t = time(NULL);
	nauty_calls++;
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();
	time_t end_t = time(NULL);

	nauty_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_layer_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_layer_calls++;
	start_c = clock();

	int index = 0;

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

	int zeroindex = 0;

	for ( int pair_ij = 0; pair_ij < nChooseK(n, 2); pair_ij++ )
	{
		/* Is it a representative? */
		if ( pair_orbits[pair_ij] == pair_ij )
		{
			/* we MOSTLY care about 2-orbits */
			if ( satdata->getAdjacency(pair_ij) == 2 )
			{
				/* this is the orbit representative */
				augment->orbitLabels[pair_ij] = index;
				augment->orbits[index] = (int*) malloc(n * n * sizeof(int));

				/* pre-fill with terminating -1's */
				for ( int k = 0; k < n * n; k++ )
				{
					augment->orbits[index][k] = -1;
				}

				/* fill orbits with pair from g */
				augment->orbits[index][0] = pair_ij;

				int k = 1;
				for ( int pair_kl = pair_ij + 1; pair_kl < nChooseK(n, 2); pair_kl++ )
				{
					if ( pair_orbits[pair_kl] == pair_ij )
					{
						/* add to orbit AND set label */
						augment->orbits[index][k] = pair_kl;
						augment->orbitLabels[pair_kl] = index;
						k++;
					}
				}

				augment->orbits[index][k] = -1;
				augment->orbits[index][k + 1] = -1;

				/* increase the pair orbit */
				index++;
			}
			else if ( satdata->getAdjacency(pair_ij) == 0 )
			{
				/* but also care about 0-orbits */
				/* a 0-type orbit rep! */
				/* this is the orbit representative */
				augment->zeroOrbitLabels[pair_ij] = zeroindex;

				for ( int pair_kl = pair_ij + 1; pair_kl < nChooseK(n, 2); pair_kl++ )
				{
					if ( pair_orbits[pair_kl] == pair_ij )
					{
						/* add to orbit AND set label */
						augment->zeroOrbitLabels[pair_kl] = zeroindex;
					}
				}

				/* increase the pair orbit */
				zeroindex++;
			}
		}
	}

	/* set numOrbits */
	augment->numOrbits = index;
	augment->numZeroEdges = zeroindex;

	/* compute canonical labels */
	augment->canonicalLabels = (int*) malloc(n * sizeof(int));

	/* compress the canonical labels */
	int cur_lab = 0;
	for ( int i = 0; i < n; i++ )
	{
		if ( lab[i] % 2 == 0 )
		{
			augment->canonicalLabels[lab[i] >> 1] = cur_lab;
			cur_lab++;
		}
	}

	/* free workspace */
	DYNFREE(workspace, worksize);
	DYNFREE(lab, lab_n);
	DYNFREE(ptn, ptn_n);
	DYNFREE(orbits, orbits_n);
	SG_FREE(canon_g);
	SG_FREE(*g);
	free(g);
	clean_generator_data_pairs();

	end_c = clock();

	orbit_layer_time += (end_c - start_c) / (double) CLOCKS_PER_SEC;
}

void use_generator_zeros(int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n)
{

	for ( int i = 0; i < nchoose2; i++ )
	{
		new_partition[i] = nchoose2 + 1;
	}

	for ( int start_pindex = 0; start_pindex < nchoose2; start_pindex++ )
	{
		int i, j;
		indexToPair(start_pindex, i, j);

		if ( new_partition[start_pindex] >= nchoose2 )
		{
			/* loop with generator! */
			new_partition[start_pindex] = start_pindex;

			/* use regular indices */
			int ni = perm[i];
			int nj = perm[j];

			int pindex = indexOf(ni, nj);

			/* wait until ordered pair match */
			while ( ni != i || nj != j )
			{
				new_partition[pindex] = start_pindex;
				ni = perm[ni];
				nj = perm[nj];
				pindex = indexOf(ni, nj);
			}
		}
	}

	joinPartitions2(pair_orbits, new_partition, nchoose2);
}

/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeZeroGraphSymmetry(SaturationData* satdata, Augmentation* augment)
{
	sparsegraph* g = compactsg(satdata->getZeroGraph());

	/* clean up old augment stuff, if necessary */
	if ( augment->numZeroGraphOrbits > 0 )
	{
		/* just return */
		return;
	}

	if ( augment->zeroGraphCanonLabels != 0 )
	{
		/* TODO: consider realloc? */
		free(augment->zeroGraphCanonLabels);
		augment->zeroGraphCanonLabels = 0;
	}

	if ( augment->zeroGraphOrbitLabels != 0 )
	{
		free(augment->zeroGraphOrbitLabels);
		augment->zeroGraphOrbitLabels = 0;
	}

	int nv = satdata->getN(); /* 2*n */
	int n = satdata->getN(); /* nv / 2 */
	int m = (nv + WORDSIZE - 1) / WORDSIZE;
	nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	for ( int i = 0; i < n; i++ )
	{
		lab[i] = i;
		ptn[i] = 1;
	}

	static DEFAULTOPTIONS_SPARSEGRAPH( options);

	init_generator_data_pairs(n);

	options.defaultptn = TRUE; /* we don't need colors */
	options.getcanon = TRUE; /* gets labels */
	options.digraph = FALSE;
	options.userautomproc = &use_generator_zeros; /* Special zeros generator method */

	options.invarproc = &adjacencies_sg; /* sparse version */

	statsblk stats; /* we'll use this at the end */
	DYNALLSTAT(setword, workspace, worksize);
	DYNALLOC1(setword, workspace, worksize, 50 * m, "malloc");

	sparsegraph canon_g;
	SG_INIT(canon_g);

	/* call nauty */
	clock_t start_c = clock();
	time_t start_t = time(NULL);
	nauty_calls++;
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();
	time_t end_t = time(NULL);

	nauty_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_zero_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_zero_calls++;

	start_c = clock();

	int index = 0;

	/* fill in the orbitLabels */
	augment->zeroGraphOrbitLabels = (int*) malloc(n * n * sizeof(int));
	for ( int i = 0; i < n * n; i++ )
	{
		augment->zeroGraphOrbitLabels[i] = -1;
	}

	int zeroindex = 0;

	for ( int pair_ij = 0; pair_ij < nChooseK(n, 2); pair_ij++ )
	{
		/* Is it a representative? */
		if ( pair_orbits[pair_ij] == pair_ij )
		{
			/* we MOSTLY care about 2-orbits */
			if ( satdata->getAdjacency(pair_ij) == 0 )
			{
				/* this is the orbit representative */
				augment->zeroGraphOrbitLabels[pair_ij] = index;

				/* increase the pair orbit */
				index++;
			}
		}
	}

	/* set numOrbits */
	augment->numZeroGraphOrbits = index;

	/* compute canonical labels */
	augment->zeroGraphCanonLabels = (int*) malloc(n * sizeof(int));

	/* compress the canonical labels */
	int cur_lab = 0;
	for ( int i = 0; i < n; i++ )
	{
		augment->zeroGraphCanonLabels[lab[i]] = cur_lab;
		cur_lab++;
	}

	/* free workspace */
	DYNFREE(workspace, worksize);
	DYNFREE(lab, lab_n);
	DYNFREE(ptn, ptn_n);
	DYNFREE(orbits, orbits_n);
	SG_FREE(canon_g);
	SG_FREE(*g);
	free(g);
	clean_generator_data_pairs();

	end_c = clock();
	orbit_zero_time += (end_c - start_c) / (double)CLOCKS_PER_SEC;
}

/*******************************************************************/
/*********** STABILIZED COMPLETION ORBIT CODE **********************/

int complete_n;
int complete_size;
int* nchooses = 0;
int** completion_orbits;
int** new_partitions;

void init_generator_data_completion(int n, int r)
{
	complete_n = n;
	complete_size = r - 2;
	nchooses = (int*) malloc((complete_size + 1) * sizeof(int));
	completion_orbits = (int**) malloc((complete_size + 1) * sizeof(int*));

	new_partitions = (int**) malloc((complete_size + 1) * sizeof(int*));

	for ( int k = 0; k <= complete_size; k++ )
	{
		nchooses[k] = nChooseK(n, k);
		completion_orbits[k] = (int*) malloc(nchooses[k] * sizeof(int));
		new_partitions[k] = (int*) malloc(nchooses[k] * sizeof(int));

		for ( int i = 0; i < nchooses[k]; i++ )
		{
			completion_orbits[k][i] = i;
		}
	}
}

void use_generator_completion(int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n)
{
	for ( int s = 1; s <= complete_size; s++ )
	{
		for ( int i = 0; i < nchooses[s]; i++ )
		{
			new_partitions[s][i] = nchooses[s] + 1;
		}

		int* start_set = (int*) malloc(s * sizeof(int));
		int* perm_set = (int*) malloc(s * sizeof(int));

		for ( int start_index = 0; start_index < nchooses[s]; start_index++ )
		{
			if ( new_partitions[s][start_index] >= nchooses[s] )
			{
				indexToSet(s, start_index, start_set);

				/* loop with generator! */
				new_partitions[s][start_index] = start_index;

				/* Use even indices */
				for ( int i = 0; i < s; i++ )
				{
					/* use perm(2*i)/2 */
					perm_set[i] = (perm[(start_set[i]) << 1]) >> 1;
				}

				int pindex = indexOfSet(s, perm_set);

				/* wait until ordered pair match */
				while ( pindex != start_index )
				{
					new_partitions[s][pindex] = start_index;

					for ( int i = 0; i < s; i++ )
					{
						/* use perm(2*i)/2 */
						perm_set[i] = (perm[(perm_set[i]) << 1]) >> 1;
					}

					pindex = indexOfSet(s, perm_set);
				}
			}
		}

		free(start_set);
		free(perm_set);

		joinPartitions2(completion_orbits[s], new_partitions[s], nchooses[s]);
	}
}

void clean_generator_data_completion()
{
	if ( nchooses != 0 )
	{
		free(nchooses);
		nchooses = 0;
	}
	if ( completion_orbits != 0 )
	{
		for ( int i = 0; i <= complete_size; i++ )
		{
			if ( completion_orbits[i] != 0 )
			{
				free(completion_orbits[i]);
				completion_orbits[i] = 0;
			}
		}
		free(completion_orbits);
		completion_orbits = 0;
	}

	if ( new_partitions != 0 )
	{
		for ( int i = 0; i <= complete_size; i++ )
		{
			if ( new_partitions[i] != 0 )
			{
				free(new_partitions[i]);
				new_partitions[i] = 0;
			}
		}
		free(new_partitions);
		new_partitions = 0;
	}
}

/**
 * Given a graph and an unassigned edge, compute all completion orbits
 * 	for that edge.
 *
 * BONUS: the completion orbits may be filtered by removing completion sets
 * 	where there is a 0-edge which would need to be a 1-edge.
 */
void computeStabilizedOrbits(int r, SaturationData* satdata, int pairindex, Augmentation* augment)
{
	sparsegraph* g = compactsg(satdata->getLayeredGraph());

	int pairi, pairj;

	indexToPair(pairindex, pairi, pairj);

	/* clean up old augment stuff, if necessary */
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

		/* TODO: consider realloc? */
		free(augment->completionOrbitReps);
		augment->completionOrbitReps = 0;
	}

	augment->stabilizer = pairindex;

	int nv = g->nv;
	int n = satdata->getN();
	int m = (nv + WORDSIZE - 1) / WORDSIZE;
	nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	/* the partition splits {i,j},{0,...(skip i,j)...,n-1},{n+i,n+j},{n,...(skip n_i, n+j)...,nv-1} */
	for ( int i = 0; i < n; i++ )
	{
		lab[i] = 2 * i;
		lab[n + i] = 2 * i + 1;

		ptn[i] = 1;
		ptn[n + i] = 1;
	}

	/* rearrange {2i,2j} and {2i+1,2j+1} */
	int temp0 = lab[0];
	lab[0] = 2 * pairi;
	lab[pairi] = temp0;
	int temp = lab[1];
	lab[1] = 2 * pairj;
	lab[pairj] = temp;

	temp0 = lab[n];
	lab[n] = 2 * pairi + 1;
	lab[n + pairi] = temp0;
	temp = lab[n + 1];
	lab[n + 1] = 2 * pairj + 1;
	lab[n + pairj] = temp;

	ptn[1] = 0;
	ptn[n - 1] = 0;
	ptn[n + 1] = 0;
	ptn[nv - 1] = 0;

	static DEFAULTOPTIONS_SPARSEGRAPH( options);

	init_generator_data_completion(n, r);

	options.defaultptn = FALSE; /* we DO need colors */
	options.getcanon = FALSE; /* we don't need the stabilized labels */
	options.digraph = FALSE;
	options.userautomproc = &use_generator_completion;

	options.invarproc = &adjacencies_sg; /* sparse version */

	statsblk stats; /* we'll use this at the end */
	DYNALLSTAT(setword, workspace, worksize);
	DYNALLOC1(setword, workspace, worksize, 50 * m, "malloc");

	sparsegraph canon_g;
	SG_INIT(canon_g);

	/* call nauty */
	clock_t start_c = clock();
	time_t start_t = time(NULL);
	nauty_calls++;
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();
	time_t end_t = time(NULL);

	nauty_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_stab_time += (double) (end_c - start_c) / (double) CLOCKS_PER_SEC;
	nauty_stab_calls++;

	start_c = clock();

	/* fill in the orbitLabels */
	int bigOrbSize = 1;

	for ( int s = 1; s <= r - 2; s++ )
	{
		bigOrbSize += nchooses[s];
	}

	augment->completionOrbitReps = (int**) malloc(bigOrbSize * sizeof(int*));
	for ( int i = 0; i < bigOrbSize; i++ )
	{
		augment->completionOrbitReps[i] = 0;
	}

	int index = 0;

	if ( n + r - 2 <= satdata->getN() )
	{
		/* we could add ALL new vertices! */
		augment->completionOrbitReps[index] = (int*) malloc((r - 2) * sizeof(int));

		for ( int i = 0; i < r - 2; i++ )
		{
			augment->completionOrbitReps[index][i] = n + i;
		}

		index++;
	}

	//	printf("\tResulting in orbits ");
	int* temp_set = (int*) malloc((r - 2) * sizeof(int));

	for ( int s = 1; s <= r - 2; s++ )
	{
		if ( n + (r - 2 - s) > satdata->getMaxN() )
		{
			/* there are too many vertices! */
			continue;
		}
		else
		{
			int base_index = index;
			int num_wo_ij = numSetsWOWO(n, s, pairi, pairj);

			for ( int set_index = 0; set_index < num_wo_ij; set_index++ )
			{
				getSetWOWO(s, pairi, pairj, set_index, temp_set);
				int full_index = indexOfSet(s, temp_set);

				/* add the orbit iff it is a representative AND can be a completion for the pair */
				if ( completion_orbits[s][full_index] == full_index )
				{
					bool can_complete = true;

					for ( int k = 0; can_complete && k < s; k++ )
					{
						int ki_index = indexOf(pairi, temp_set[k]);
						int kj_index = indexOf(pairj, temp_set[k]);

						if ( satdata->getAdjacency(ki_index) == 0 || satdata->getAdjacency(kj_index) == 0 )
						{
							can_complete = false;
						}

						for ( int l = 0; can_complete && l < k; l++ )
						{
							int kl_index = indexOf(temp_set[k], temp_set[l]);
							if ( satdata->getAdjacency(kl_index) == 0 )
							{
								can_complete = false;
							}
						}
					}

					if ( can_complete )
					{
						/* this is a good completion */
						augment->completionOrbitReps[index] = (int*) malloc((r - 2) * sizeof(int));
						indexToSet(s, full_index, augment->completionOrbitReps[index]);

						/* fill in the rest */
						for ( int t = s; t < r - 2; t++ )
						{
							augment->completionOrbitReps[index][t] = n + (t - s);
						}

						/* increase the set orbit */
						index++;
					}
				}
			}
		}
	}

	/* set numOrbits */
	augment->numStabilizedOrbits = index;
	free(temp_set);
	temp_set = 0;

	/* free workspace */
	DYNFREE(workspace, worksize);
	DYNFREE(lab, lab_n);
	DYNFREE(ptn, ptn_n);
	DYNFREE(orbits, orbits_n);
	SG_FREE(canon_g);
	SG_FREE(*g);
	free(g);
	clean_generator_data_completion();

	end_c = clock();

	orbit_stab_time += (end_c - start_c) / (double) CLOCKS_PER_SEC;
}
