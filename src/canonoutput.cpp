/*
 * canonoutput.c
 *
 *  Created on: Jun 29, 2010
 *      Author: derrickstolee
 *
 *  Take a list of sparse6 representations
 *  	and output the list of sparse6 canonical
 *  	representations.
 */

#include <stdio.h>
#include "nauty.h"
#include "nausparse.h"
#include "nautinv.h"
#include "gtools.h"

int main(int argc, char** argv)
{
	char buffer[1500];
	int print_orig = 0;

	if ( argc > 1 )
	{
		print_orig = 1;
	}

	while ( scanf("%s", buffer) != EOF )
	{
		/* buffer holds an s6 string */
		sparsegraph g;
		SG_INIT(g);

		int num_loops;
		stringtosparsegraph(buffer, &g, &num_loops);

		int nv = g.nv;
		int m = (nv + WORDSIZE - 1) / WORDSIZE;
		nauty_check(WORDSIZE, m, nv, NAUTYVERSIONID);

		DYNALLSTAT(int, lab, lab_n);
		DYNALLSTAT(int, ptn, ptn_n);
		DYNALLSTAT(int, orbits, orbits_n);

		DYNALLOC1(int, lab, lab_n, nv, "malloc");
		DYNALLOC1(int, ptn, ptn_n, nv, "malloc");
		DYNALLOC1(int, orbits, orbits_n, nv, "malloc");

		static DEFAULTOPTIONS_SPARSEGRAPH( options);

		options.defaultptn = TRUE; /* Don't need colors */
		options.getcanon = TRUE; /* gets labels */
		options.digraph = TRUE;

		statsblk stats; /* we'll use this at the end */
		DYNALLSTAT(setword, workspace, worksize);
		DYNALLOC1(setword, workspace, worksize, 50 * m, "malloc");

		sparsegraph canon_g;
		SG_INIT(canon_g);

		/* call nauty */
		nauty((graph*) &g, lab, ptn, NULL, orbits, &options, &stats, workspace,
				50 * m, m, g.nv, (graph*) &canon_g);

		sortlists_sg(&canon_g);

		char* canon_str = sgtos6(&canon_g);
		int canon_len = strlen(canon_str);

		if ( print_orig == 0 )
		{
			printf("%s", canon_str);
		}
		else
		{
			canon_str[canon_len-1] = 0;
			printf("%s", canon_str);
			printf("\t%s\n", buffer);
		}

		/* free workspace */
		DYNFREE(workspace, worksize);
		DYNFREE(lab,lab_n);
		DYNFREE(ptn,ptn_n);
		DYNFREE(orbits,orbits_n);

		SG_FREE(canon_g);
		SG_FREE(g);
	}
}
