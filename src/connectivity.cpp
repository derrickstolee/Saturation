/*
 * connectivity.cpp
 *
 *  Created on: Dec 18, 2010
 *      Author: derrickstolee
 */

#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gtools.h"
#include "connectivity.hpp"

bool findPath( sparsegraph* g, int v1, int v2, char* without_ear, int without_vert )
{
	int len = 0;

	for ( len = 0; without_ear[len] >= 0; len++ )
		;

	if ( without_ear[0] == v2 && without_ear[len - 1] == v1 )
	{
		int t = v1;
		v1 = v2;
		v2 = t;
	}

	std::list<int> q;
	char* visited = (char*) malloc(g->nv);
	bzero(visited, g->nv);

	q.push_back(v1);

	for ( int i = 0; without_ear[i] != v2 && without_ear[i] >= 0; i++ )
	{
		visited[without_ear[i]] = 1;
	}
	visited[without_vert] = 1;

	while ( q.size() > 0 )
	{
		int v = q.front();
		q.pop_front();

		for ( int i = 0; i < g->d[v]; i++ )
		{
			int u = g->e[g->v[v] + i];

			if ( u == v2 )
			{
				if ( v != v1 || len > 2 )
				{
					/* important not to use edge from v1 to v2 */
					free(visited);
					return true;
				}
				/* otherwise, ignore that we were here */
			}
			else if ( visited[u] == 0 )
			{
				visited[u] = 1;
				/* add to queue */
				q.push_back(u);
			}
		}
	}

	free(visited);
	return false;
}

/**
 * is2connected
 *
 * Is the graph still 2-connected without the specified ear?
 *
 * @param g an undirected sparsegraph
 * @param v1 an index of a vertex
 * @param v2 an index of another vertex
 * @param without_ear an array of vertices, -1 terminated.
 */
bool is2connected( sparsegraph* g, int v1, int v2, char* without_ear )
{
	/* TODO: the RIGHT way to do this is Max Flow */

	/* check that there are no cut vertices in g - without_ear */
	for ( int i = 0; i < g->nv; i++ )
	{
		bool in_ear = false;

		for ( int j = 0; without_ear[j] >= 0; j++ )
		{
			if ( i == without_ear[j] )
			{
				in_ear = true;
				break;
			}
		}

		if ( !in_ear && findPath(g, v1, v2, without_ear, i) == false )
		{

			return false;
		}
	}

	return true;
}

/**
 * isCutEdge
 *
 * Is the given edge a cut edge?
 */
bool isCutEdge( sparsegraph* g, int v1, int v2 )
{
	std::list<int> q;
	char* visited = (char*) malloc(g->nv);
	bzero(visited, g->nv);

	q.push_back(v1);
	visited[v1] = 1;

	while ( q.size() > 0 )
	{
		int v = q.front();
		q.pop_front();

		for ( int i = 0; i < g->d[v]; i++ )
		{
			int u = g->e[g->v[v] + i];

			if ( u == v2 )
			{
				if ( v != v1  )
				{
					/* important not to use edge from v1 to v2 */
					free(visited);
					return false;
				}
				/* otherwise, ignore that we were here */
			}
			else if ( visited[u] == 0 )
			{
				visited[u] = 1;
				/* add to queue */
				q.push_back(u);
			}
		}
	}

	free(visited);
	return true;
}
