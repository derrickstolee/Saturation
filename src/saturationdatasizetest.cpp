/*
 * saturationdatasizetest.cpp
 *
 *  Created on: Apr 15, 2011
 *      Author: stolee
 */

#include "SaturationData.hpp"

int main( void )
{
	int n = 10;
	int N = 30;

	int r = 4;
	int R = 7;

	for ( int j = r; j <= R; j++ )
	{
		printf("\t\tr = %d\n", j);
		for ( int i = n; i <= N; i++ )
		{
			delete new SaturationData(j, i);
		}
		printf("\n");
	}

	return 0;
}
