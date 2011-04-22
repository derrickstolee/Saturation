/*
 * SaturationSymmetry.cpp
 *
 *  Created on: Apr 10, 2011
 *      Author: stolee
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "SaturationSymmetry.hpp"
#include "translation.hpp"
#include "nausparse.h"
#include "nauty.h"
#include "gtools.h"

/******** PARTITION JOINING METHODS ***********/

int root( int i, int* p1 )
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
void joinPartitions2( int* p1, int* p2, int ne )
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

void init_generator_data_pairs( int n )
{
	gen_n = n;
	nchoose2 = nChooseK(n, 2);
	pair_orbits = (int*) malloc(nchoose2 * sizeof(int));

	for ( int i = 0; i < nchoose2; i++ )
	{
		pair_orbits[i] = i;
	}
}

void use_generator_pairs( int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n )
{
	int* new_partition = (int*) malloc(nchoose2 * sizeof(int));

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

	free(new_partition);
}

void clean_generator_data_pairs()
{
	nchoose2 = 0;
	if ( pair_orbits != 0 )
	{
		free(pair_orbits);
		pair_orbits = 0;
	}
}

/**
 * Given a graph, fill in the unassigned orbits for
 */
void computeUnassignedOrbits( SaturationGraph* satgraph, Augmentation* augment )
{
	sparsegraph* g = satgraph->getSmallG();

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

	int nv = g->nv;
	int n = satgraph->getN();
	int m = (nv + WORDSIZE - 1) / WORDSIZE;
	nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	for ( int i = 0; i < nv; i++ )
	{
		lab[i] = i;
		ptn[i] = 1;
	}

	/* the partition splits {0,..,n-1},{n,...,nv-1} */
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
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();
	time_t end_t = time(NULL);

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
			if ( satgraph->getAdjacency(pair_ij) == 2 )
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
			else if ( satgraph->getAdjacency(pair_ij) == 0 )
			{
				/* but also care about 1-orbits */
				/* a 1-type orbit rep! */
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

	/* compute canonical labels */
	augment->canonicalLabels = (int*) malloc(n * sizeof(int));

	for ( int i = 0; i < n; i++ )
	{
		augment->canonicalLabels[lab[i]] = i;
	}

	/* free workspace */
	DYNFREE(workspace, worksize);
	DYNFREE(lab, lab_n);
	DYNFREE(ptn, ptn_n);
	DYNFREE(orbits, orbits_n);
	SG_FREE(canon_g);
	clean_generator_data_pairs();
}

/*******************************************************************/
/*********** STABILIZED COMPLETION ORBIT CODE **********************/

int complete_n;
int complete_size;
int* nchooses = 0;
int** completion_orbits;

void init_generator_data_completion( int n, int r )
{
	complete_n = n;
	complete_size = r - 2;
<<<<<<< HEAD
	nchooses = (int*) malloc((complete_size + 1) * sizeof(int));
	completion_orbits = (int**) malloc((complete_size + 1) * sizeof(int*));
=======
	nchooses = (int*) malloc((complete_size+1) * sizeof(int));
	completion_orbits = (int**) malloc((complete_size+1) * sizeof(int*));
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

	for ( int k = 0; k <= complete_size; k++ )
	{
		nchooses[k] = nChooseK(n, k);
		completion_orbits[k] = (int*) malloc(nchooses[k] * sizeof(int));
		for ( int i = 0; i < nchooses[k]; i++ )
		{
			completion_orbits[k][i] = i;
		}
	}
}

void use_generator_completion( int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n )
{
	for ( int s = 1; s <= complete_size; s++ )
	{
		int* new_partition = (int*) malloc(nchooses[s] * sizeof(int));

		for ( int i = 0; i < nchooses[s]; i++ )
		{
			new_partition[i] = nchooses[s] + 1;
		}

		int* start_set = (int*) malloc(s * sizeof(int));
		int* perm_set = (int*) malloc(s * sizeof(int));

		for ( int start_index = 0; start_index < nchooses[s]; start_index++ )
		{
			if ( new_partition[start_index] >= nchooses[s] )
			{
				indexToSet(s, start_index, start_set);

				/* loop with generator! */
				new_partition[start_index] = start_index;

				for ( int i = 0; i < s; i++ )
				{
					perm_set[i] = perm[start_set[i]];
				}

				int pindex = indexOfSet(s, perm_set);
<<<<<<< HEAD
				//
				//				bool are_equal = true;
				//
				//				for ( int i = 0;are_equal && i < s; i++ )
				//				{
				//					if ( perm_set[i] != start_set[i] )
				//					{
				//						are_equal = false;
				//					}
				//				}
=======
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

				/* wait until ordered pair match */
				while ( pindex != start_index )
				{
					new_partition[pindex] = start_index;

					for ( int i = 0; i < s; i++ )
					{
						perm_set[i] = perm[perm_set[i]];
					}

					pindex = indexOfSet(s, perm_set);
				}
			}
		}

		free(start_set);
		free(perm_set);

		joinPartitions2(completion_orbits[s], new_partition, nchooses[s]);

		free(new_partition);
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
}

/**
 * Given a graph and an unassigned edge, compute all completion orbits
 * 	for that edge.
 *
 * BONUS: the completion orbits may be filtered by removing completion sets
 * 	where there is a 0-edge which would need to be a 1-edge.
 */
void computeStabilizedOrbits( int r, SaturationGraph* satgraph, int pairindex, Augmentation* augment )
{
	sparsegraph* g = satgraph->getSmallG();

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
	int n = satgraph->getN();
	int m = (nv + WORDSIZE - 1) / WORDSIZE;
	nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);
<<<<<<< HEAD

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	/* the partition splits {i,j},{0,...(skip i,j)...,n-1},{n+i,n+j},{n,...(skip n_i, n+j)...,nv-1} */
	for ( int i = 0; i < nv; i++ )
	{
		lab[i] = i;
		ptn[i] = 1;
	}

	/* rearrange {i,j} and {n+i,n+j} */
	lab[0] = pairi;
	lab[pairi] = 0;
	int temp = lab[1];
	lab[1] = pairj;
	lab[pairj] = temp;

	lab[n] = n + pairi;
	lab[n + pairi] = n;
	temp = lab[n + 1];
	lab[n + 1] = n + pairj;
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
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();

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

	if ( n + r - 2 <= satgraph->getN() )
	{
		/* we could add ALL vertices! */
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
		if ( n + (r - 2 - s) > satgraph->getMaxN() )
		{
			//			printf("|{%d}| = X  ", s);
			/* there are too many vertices! */
			continue;
		}

		int base_index = index;
		for ( int set_index = 0; set_index < nchooses[s]; set_index++ )
		{
			/* add the orbit iff it is a representative AND can be a completion for the pair */
			if ( completion_orbits[s][set_index] == set_index )
			{
				bool can_complete = true;
				indexToSet(s, set_index, temp_set);

				/* check if the set contains i or j */
				for ( int k = 0; can_complete && k < s; k++ )
				{
					if ( temp_set[k] == pairi || temp_set[k] == pairj )
					{
						can_complete = false;
					}
				}

				for ( int k = 0; can_complete && k < s; k++ )
				{
					int ki_index = indexOf(pairi, temp_set[k]);
					int kj_index = indexOf(pairj, temp_set[k]);

					if ( satgraph->getAdjacency(ki_index) == 0 || satgraph->getAdjacency(kj_index) == 0 )
					{
						can_complete = false;
					}

=======

	DYNALLSTAT(int, lab, lab_n);
	DYNALLSTAT(int, ptn, ptn_n);
	DYNALLSTAT(int, orbits, orbits_n);

	DYNALLOC1(int, lab, lab_n, nv, "malloc");
	DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
	DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

	/* create labels and partitions */
	/* the partition splits {i,j},{0,...(skip i,j)...,n-1},{n+i,n+j},{n,...(skip n_i, n+j)...,nv-1} */
	for ( int i = 0; i < nv; i++ )
	{
		lab[i] = i;
		ptn[i] = 1;
	}

	/* rearrange {i,j} and {n+i,n+j} */
	lab[0] = pairi;
	lab[pairi] = 0;
	int temp = lab[1];
	lab[1] = pairj;
	lab[pairj] = temp;

	lab[n] = n + pairi;
	lab[n + pairi] = n;
	temp = lab[n + 1];
	lab[n + 1] = n + pairj;
	lab[n + pairj] = temp;

	ptn[1] = 0;
	ptn[n - 1] = 0;
	ptn[n + 1] = 0;
	ptn[nv - 1] = 0;

//	printf("--[computeStabilizedOrbits] Stabilizing %d %d with Labels ", pairi, pairj);

//	for ( int i = 0; i < nv; i++ )
//	{
//		printf("%2d", lab[i]);
//		if ( ptn[i] == 0 )
//		{
//			printf("|");
//		}
//		else
//		{
//			printf(" ");
//		}
//	}

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
	nauty((graph*) g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50 * m, m, nv, (graph*) &canon_g);
	clock_t end_c = clock();
	time_t end_t = time(NULL);

	/* fill in the orbitLabels */
	int bigOrbSize = 1;

	for ( int s = 1; s <= r-2; s++ )
	{
		bigOrbSize += nchooses[s];
	}

	augment->completionOrbitReps = (int**) malloc(bigOrbSize * sizeof(int*));
	for ( int i = 0; i < bigOrbSize; i++ )
	{
		augment->completionOrbitReps[i] = 0;
	}

	int index = 0;

	if ( n + r - 2 <= satgraph->getN() )
	{
		/* we could add ALL vertices! */
		augment->completionOrbitReps[index] = (int*)malloc((r-2)*sizeof(int));

		for ( int i = 0; i < r-2; i++ )
		{
			augment->completionOrbitReps[index][i] = n + i;
		}

		index++;
	}

//	printf("\tResulting in orbits ");
	int* temp_set = (int*) malloc((r - 2) * sizeof(int));

	for ( int s = 1; s <= r - 2; s++ )
	{
		if ( n + (r-2-s) > satgraph->getMaxN() )
		{
//			printf("|{%d}| = X  ", s);
			/* there are too many vertices! */
			continue;
		}

		int base_index = index;
		for ( int set_index = 0; set_index < nchooses[s]; set_index++ )
		{
			/* add the orbit iff it is a representative AND can be a completion for the pair */
			if ( completion_orbits[s][set_index] == set_index )
			{
				bool can_complete = true;
				indexToSet(s, set_index, temp_set);

				/* check if the set contains i or j */
				for ( int k = 0; can_complete && k < s; k++ )
				{
					if ( temp_set[k] == pairi || temp_set[k] == pairj )
					{
						can_complete = false;
					}
				}

				for ( int k = 0; can_complete && k < s; k++ )
				{
					int ki_index = indexOf(pairi, temp_set[k]);
					int kj_index = indexOf(pairj, temp_set[k]);

					if ( satgraph->getAdjacency(ki_index) == 0 || satgraph->getAdjacency(kj_index) == 0 )
					{
						can_complete = false;
					}

>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
					for ( int l = k + 1; can_complete && l < s; l++ )
					{
						int kl_index = indexOf(temp_set[k], temp_set[l]);
						if ( satgraph->getAdjacency(kl_index) == 0 )
						{
							can_complete = false;
						}
					}
				}

				if ( can_complete )
				{
					/* this is a good completion */
					augment->completionOrbitReps[index] = (int*) malloc((r - 2) * sizeof(int));
					indexToSet(s, set_index, augment->completionOrbitReps[index]);

					/* fill in the rest */
<<<<<<< HEAD
					for ( int t = s; t < r - 2; t++ )
					{
						augment->completionOrbitReps[index][t] = n + (t - s);
=======
					for ( int t = s; t < r-2; t++ )
					{
						augment->completionOrbitReps[index][t] = n + (t-s);
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
					}

					/* increase the set orbit */
					index++;
				}
			}
		}

<<<<<<< HEAD
		//		printf("|{%d}| = %d | ", s, index - base_index );
=======
//		printf("|{%d}| = %d | ", s, index - base_index );
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
	}

	/* set numOrbits */
	augment->numStabilizedOrbits = index;
	free(temp_set);
	temp_set = 0;

<<<<<<< HEAD
	//	printf("(%d total)\n", index);
=======
//	printf("(%d total)\n", index);
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

	/* free workspace */
	DYNFREE(workspace, worksize);
	DYNFREE(lab, lab_n);
	DYNFREE(ptn, ptn_n);
	DYNFREE(orbits, orbits_n);
	SG_FREE(canon_g);
	clean_generator_data_completion();
}
