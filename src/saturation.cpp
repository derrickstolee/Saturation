/*
 * saturation.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: stolee
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "SaturationManager.hpp"
#include "SaturationAugmenter.hpp"
#include "SaturationData.hpp"

int main( int argc, char** argv )
{
	int r = -1;
	int N = -1;
	int mode = _MODE_CANONICAL_;
	//	SaturationManager* manager = 0;

	for ( int i = 1; i < argc - 1; i++ )
	{
		if ( argv[i][0] == '-' && argv[i][1] != 0 && argv[i][2] == 0 )
		{
			if ( argv[i][1] == 'r' )
			{
				r = atoi(argv[i + 1]);
			}
			else if ( argv[i][1] == 'N' )
			{
				N = atoi(argv[i + 1]);
			}
		}
		else if ( strcmp(argv[i], "--mode") == 0 )
		{
			if ( strcmp(argv[i + 1], "canon") == 0 )
			{
				mode = _MODE_CANONICAL_;
			}
			else if ( strcmp(argv[i + 1], "orbital") == 0 )
			{
				mode = _MODE_ORBITAL_;
			}
		}
	}

	if ( r < 0 || N < 0 )
	{
		printf("Usage: saturation.exe -r # -N # [TreeSearch args]\n");
		exit(1);
	}

	//	SaturationData* satdata = new SaturationData(r, N);

	//	manager = new SaturationManager(r, N);

	SaturationAugmenter* augmenter = new SaturationAugmenter(r, N, mode);
	augmenter->importArguments(argc, argv);

	while ( augmenter->readJob(stdin) >= 0 )
	{
		augmenter->doSearch();
	}

//	delete manager;
//	delete satdata;
	delete augmenter;

	return 0;
}
