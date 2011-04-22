/*
 * translation.cpp
 *
 *  Created on: Apr 10, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "translation.hpp"

/**
 * indexOf
 *
 * The co-lex index of the pair (i,j).
 */
int indexOf(int i, int j)
{
	if ( i == j )
	{
		/* bad input */
		return -1;
	}

	if ( j < i )
	{
		/* wrong order */
		return indexOf(j, i);
	}

	/* CO-LEX ORDER */

	/* THERE ARE (j choose 2) SETS WITH NUMBERS AT MOST j */
	return nChooseK(j, 2) + i;
}

/**
 * indexToPair
 */
void indexToPair(int index, int& i, int& j)
{
	/* Need to find largest integer j so that (j*(j-1))/2 = (1/2)j^2 - (1/2)j <= index */
	/* That is, x^2 - x - 2*index = 0 */
	/* x = (1 + sqrt(1-4(-2*index)))/2 */
	/* Then, round down */
	for ( j = 1; nChooseK(j, 2) <= index; j++ )
		;
	j--;

	/* i is the "remainder" */
	i = index - nChooseK(j, 2);
}

int recurseIndexOfSet(int size, int* set)
{
	if ( size <= 0 )
	{
		return 0;
	}
	if ( size == 1 )
	{
		return set[0];
	}
	else if ( size == 2 )
	{
		return indexOf(set[0], set[1]);
	}

	/* count all sets of this size with smaller values from the last one */
	int base = nChooseK(set[size - 1], size);
	int offset = recurseIndexOfSet(size - 1, set);

	return base + offset;
}

/**
 * indexOfSet
 *
 * Get the co-lex order of the set of a given size.
 */
int indexOfSet(int size, int* set)
{
	/* ensure sortedness */
	sortSet(size, set);

	return recurseIndexOfSet(size, set);
}

/**
 * indexToSet
 */
void indexToSet(int size, int index, int* set)
{
	for ( int i = size - 1; i >= 0; i-- )
	{
		/* find ith bit by considering largest portion of index */
		/* then lower for the next position */

		/* we need to solve for largest s with (s choose (i+1)) <= index */
		int min_elt = i;
		int max_elt = 2 * (i + 1);

		/* find the right frame */
		while ( nChooseK(max_elt, i + 1) <= index )
		{
			min_elt = max_elt;
			max_elt = max_elt << 1; /* double */
		}

		/* do binary search */
		while ( min_elt <= max_elt )
		{
			int half = (max_elt + min_elt) >> 1;

			if ( nChooseK(half, i + 1) <= index )
			{
				min_elt = half + 1;
			}
			else
			{
				max_elt = half - 1;
			}
		}

		/* place this value */
		set[i] = min_elt - 1;

		/* modify index */
		index -= nChooseK(set[i], i + 1);
	}
}

void testAndSwap(int* set, int i, int j)
{
	if ( i < j && set[i] > set[j] )
	{
		int temp = set[i];
		set[i] = set[j];
		set[j] = temp;
	}
}

/**
 * sortSet
 *
 * This is implemented only to size <= 5, for the case R <= 7.
 */
void sortSet(int size, int* set)
{
	if ( size <= 1 )
	{
		/* do nothing */
	}
	else if ( size == 2 )
	{
		testAndSwap(set, 0, 1);
	}
	else if ( size == 3 )
	{
		testAndSwap(set, 0, 1);
		testAndSwap(set, 0, 2);
		testAndSwap(set, 1, 2);
	}
	else if ( size == 4 )
	{
		testAndSwap(set, 0, 1);
		testAndSwap(set, 0, 2);
		testAndSwap(set, 0, 3);
		testAndSwap(set, 1, 2);
		testAndSwap(set, 1, 3);
		testAndSwap(set, 2, 3);
	}
	else if ( size == 5 )
	{
		testAndSwap(set, 0, 1);
		testAndSwap(set, 0, 2);
		testAndSwap(set, 0, 3);
		testAndSwap(set, 0, 4);
		testAndSwap(set, 1, 2);
		testAndSwap(set, 1, 3);
		testAndSwap(set, 1, 4);
		testAndSwap(set, 2, 3);
		testAndSwap(set, 2, 4);
		testAndSwap(set, 3, 4);
	}
	else if ( size == 6 )
	{

		testAndSwap(set, 0, 1);
		testAndSwap(set, 0, 2);
		testAndSwap(set, 0, 3);
		testAndSwap(set, 0, 4);
		testAndSwap(set, 0, 5);
		testAndSwap(set, 1, 2);
		testAndSwap(set, 1, 3);
		testAndSwap(set, 1, 4);
		testAndSwap(set, 1, 5);
		testAndSwap(set, 2, 3);
		testAndSwap(set, 2, 4);
		testAndSwap(set, 2, 5);
		testAndSwap(set, 3, 4);
		testAndSwap(set, 3, 5);
		testAndSwap(set, 4, 5);
	}
	else
	{
		printf("--[sortSet] NOT IMPLEMENTED FOR SIZE %d.\n", size);
		exit(1);
	}
}

/**
 * nChooseK
 */
int nChooseK(int n, int k)
{
	double N = n;
	double K = k;

	double nCk = 1.0;

	for ( double i = 0; i < k; i++ )
	{
		nCk *= ((N - i) / (K - i));
	}

	return (int) round(nCk);
}

/**
 * indexOfTuple
 *
 * Get the product order of the set of a given size.
 */
<<<<<<< HEAD
int indexOfTuple(int n, int size, int* tuple)
=======
int indexOfTuple( int n, int size, int* tuple )
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
{
	if ( size == 0 )
	{
		return 0;
	}
	else if ( size == 1 )
	{
		return tuple[0];
	}
	else if ( size == 2 )
	{
		return n * tuple[1] + tuple[0];
	}
	else
	{
		int base = 0;
		if ( tuple[size - 1] > 0 )
		{
			base = pow(n, size - 1) * tuple[size - 1];
		}

		int offset = indexOfTuple(n, size - 1, tuple);

		return base + offset;
	}
}

/**
 * indexToTuple
 */
<<<<<<< HEAD
void indexToTuple(int n, int size, int index, int* tuple)
=======
void indexToTuple( int n, int size, int index, int* tuple )
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
{
	if ( size == 0 )
	{
		return;
	}
	else if ( size == 1 )
	{
		tuple[0] = index;
	}
	else
	{
		int npow = pow(n, size - 1);
<<<<<<< HEAD
		tuple[size - 1] = (int) floor(index / npow);
=======
		tuple[size - 1] = (int)floor(index / npow);
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
		int base = npow * tuple[size - 1];
		int offset = index - base;
		indexToTuple(n, size - 1, offset, tuple);
	}
}
<<<<<<< HEAD

/**
 * numSetsW
 *
 * Get the number of sets with i
 */
int numSetsW(int n, int s, int i)
{
	return nChooseK(n - 1, s - 1);
}

/**
 * numSetsWW
 *
 * Get the number of sets of size s
 * 	within [n] with BOTH entries
 */
int numSetsWW(int n, int s, int i, int j)
{
	return nChooseK(n - 2, s - 2);
}

/**
 * numSetsWWO
 *
 * Get the number of sets with first entry, but NOT second.
 */
int numSetsWWO(int n, int s, int i, int j)
{
	return nChooseK(n - 2, s - 1);
}

/**
 * numSetsWOWO
 *
 * Get the number of sets without either entry.
 */
int numSetsWOWO(int n, int s, int i, int j)
{
	return nChooseK(n - 2, s);
}

/**
 * getSetW
 */
void getSetW(int s, int i, int index, int* set)
{
	set[0] = i;

	indexToSet(s - 1, index, &(set[1]));

	for ( int k = 1; k < s; k++ )
	{
		if ( set[k] >= i )
		{
			set[k] = set[k] + 1;
		}
	}
}

/**
 * getSetWW
 */
void getSetWW(int s, int i, int j, int index, int* set)
{
	set[0] = i;
	set[1] = j;

	getSetWOWO(s - 2, i, j, index, &(set[2]));

	//	printf("--getSetWW: ");
	//	for ( int k = 0; k < s; k++ )
	//	{
	//		printf("%d ", set[k]);
	//	}
	//	printf("\n");
}

/**
 * getSetWWO
 */
void getSetWWO(int s, int i, int j, int index, int* set)
{
	set[0] = i;

	getSetWOWO(s - 1, i, j, index, &(set[1]));
}

/**
 * getSetWOWO
 */
void getSetWOWO(int s, int i, int j, int index, int* set)
{
	indexToSet(s, index, set);

	if ( i > j )
	{
		int t = i;
		i = j;
		j = t;
	}

	for ( int k = 0; k < s; k++ )
	{
		if ( set[k] >= i && (k == 0 || set[k - 1] < i) )
		{
			for ( int l = k; l < s; l++ )
			{
				set[l] = set[l] + 1;
			}
		}

		if ( set[k] >= j && (k == 0 || set[k - 1] < j) )
		{
			for ( int l = k; l < s; l++ )
			{
				set[l] = set[l] + 1;
			}
		}
	}
}
=======
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
