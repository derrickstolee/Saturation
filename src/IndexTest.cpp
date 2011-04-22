/*
 * IndexTest.cpp
 *
 *  Created on: Apr 8, 2011
 *      Author: stolee
 */

#include <math.h>
#include <stdio.h>

/**
 * indexOf -- Get the index of the i-j edge in the adjacency matrix.
 *
 * @param i one of the vertices
 * @param j another vertex
 * @return a position in the adjacency matrix for the i-j edge. If i==j or either is negative, returns -1.
 */
int indexOf( int i, int j )
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
void indexToPair( int index, int& i, int& j )
{
	/* Need to find integer j so that (j*(j-1))/2 = (1/2)j^2 - (1/2)j = index */
	/* That is, x^2 - x - 2*index = 0 */
	/* x = (1 + sqrt(1-4(-2*index)))/2 */
	/* Then, round down */
	j = (int) ((1.0 + sqrt(1.0 + 8.0 * (double) index)) * 0.5);

	/* i is the "remainder" */
	i = index - (j * (j - 1)) / 2;
}

int main( void )
{
	int n = 10000;

	for ( int i = 0; i < n; i++ )
	{
		for ( int j = 0; j < n; j++ )
		{
			if ( i != j )
			{
				int index = indexOf(i, j);

				int k, l;

				indexToPair(index, k, l);

				if ( (i != k || j != l) && (i != l || j != k) )
				{
					/* error in translation! */
					printf("ERROR: %2d,%2d\t->\t%2d\t%2d,%2d\n", i, j, index, k, l);
				}
			}
		}
	}
//
//	for ( int index = 0; index < n * (n - 1) / 2; index++ )
//	{
//		int i, j;
//		indexToPair(index,i, j);
//		printf("%2d\t%2d,%2d\n", index, i, j);
//	}

	return 0;
}
