/*
 * SaturationGraph.hpp
 *
 *  Created on: Apr 6, 2011
 *      Author: stolee
 */

#ifndef SATURATIONGRAPH_HPP_
#define SATURATIONGRAPH_HPP_

#include <stack>

#include "nauty.h"
#include "nausparse.h"
#include "gtools.h"

#include "Augmentation.hpp"
#include "TreeSet.hpp"

/**
 * SaturationGraph class.
 *
 * Store a "stack" of graphs, along with the augmentations to build them
 * step-by-step, for an in-place data structure in the saturation search.
 *
 * The adjacency matrix is stored with positions in co-lex order.
 * The method 'indexOf(int i, int j)' gives the position of the i-j edge, which
 * 	is used for computing pair-orbits as well.
 * It is inverted by the method 'indexToPair(int index, int& i, int& j)'.
 * The adjacency matrix stores characters, which in turn store values of 0, 1, or 2.
 *
 * Value 2: Unassigned edge.
 * Value 1: Assigned edge.
 * Value 0: Assigned non-edge, which has a unique K_r-completion.
 */
class SaturationGraph
{
	protected:
		/**
		 * r -- We are looking for a uniquely K_r-saturated graph.
		 */
		int r;

		/**
		 * n -- The current number of vertices.
		 */
		int n;

		/**
		 * N -- the MAXIMUM number of vertices.
		 */
		int N;

		/**
		 * adjmat -- The 0-1-2 adjacency matrix in co-lex order.
		 */
		char* adjmat;

		/**
		 * completemult -- The multiplicity of each 1-edge's use
		 * 	in a completion.
		 */
		char* completemult;

		/**
		 * completions -- The completion set for each 0-edge.
		 * It is ordered in co-lex order, to be accessible by the pair index.
		 * If the edge is not a 0-edge, the completion is null.
		 */
		int** completions;

		/**
		 * zeroDegress -- The degree sequence of 0-edges at each vertex.
		 */
		int* zeroDegrees;

		/**
		 * oneDegress -- The degree sequence of 1-edges at each vertex.
		 */
		int* oneDegrees;


		/**
		 * augmentations -- The stack of augmentations.
		 */
		std::stack<Augmentation*> augmentations;

		/**
		 * augmentation_succeeded Did the last augmentation succeed?
		 *
		 * False iff something went wrong, like too many vertices
		 * 	or a copy of K_r or a second completion.
		 */
		bool augmentation_succeeded;

		/**
		 * g -- A sparsegraph representation of the currently-assigned edges.
		 *
		 * There are 2N vertices, although 2n are used at a time.
		 *
		 * Vertices are partitioned by 0...(N-1) | N...(2N-1).
		 *
		 * The vertex i is connected to N+i.
		 *
		 * The vertices 0...(N-1) have edges between them corresponding to 0-type
		 * 	edges.
		 *
		 * The vertices N...(2N-1) have edges between them corresponding to 1-type
		 * 	edges.
		 *
		 * The positions v[i] are given by multiples of N, since there
		 * 	are at most N neighbors to each vertex.
		 *
		 * The degree d[i] is 1 (to connect to N+i%2N)
		 * 	plus the current degree of that vertex within that type.
		 *
		 * The most recent edges are in the last parts of each list e[i].
		 */
		sparsegraph* g;

		/**
		 * small_g -- A compacted version of g, without wasted space.
		 *
		 * This is used to compute orbits using nauty.
		 */
		sparsegraph* small_g;

		/**
		 * g_updated -- True iff g and small_g do not agree.
		 */
		bool g_updated;

		/**
		 * openedges -- A set of pair indices for open edges.
		 */
		TreeSet* openedges;


		bool do_zero_implications;

		/********************* statistics */
		double time_in_regenerate;

		double time_in_compact;

		double time_in_orbits;

		double time_in_stabilized;

		double time_in_feasible;

		/********************* protected methods */

		/**
		 * compactTheGraph() Compact the graph g into small_g.
		 */
		void compactTheGraph();

		/**
		 * Initialize all data elements.
		 *
		 * @param r the saturation parameter.
		 */
		void initData(int r, int N);

		/**
		 * convertEdge
		 *
		 * Take an edge index, and turn it into a new type.
		 */
		void convertEdge(int index, char type);

		/**
		 * regenerateOpenEdges
		 *
		 * THIS IS A TEMPORARY METHOD
		 */
		void regenerateOpenEdges();

	public:
		/**
		 * Initialize a SaturationGraph with saturation parameter r.
		 *
		 * As a base element, there exists a non-edge with a completion, using r vertices,
		 * 	then adds the (r+1)th vertex with all incident edges of type 2.
		 */
		SaturationGraph(int r, int N);

		/**
		 * Destructor
		 */
		virtual ~SaturationGraph();

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
		bool augment(int noni, int nonj, int* completion);

		/**
		 * fill -- Fill all unassigned edges with assigned edges.
		 *
		 * @return True iff the augmentation succeeds without adding a K_r
		 * 	or multiple K_r's when adding some assigned non-edge.
		 */
		bool fill();

		/**
		 * addVertex -- Add a vertex to the graph with all edges unassigned.
		 */
		bool addVertex();

		/**
		 * pop() -- Remove the top of the stack, reverting to the lower state of the graph.
		 */
		void pop();


		/****************** SYMMETRY METHODS ***************************/

		/**
		 * computeOrbits
		 *
		 * Use the current assignment of edges to generate a
		 * list of pair orbits.
		 */
		void computeOrbits();

		/**
		 * numOrbits -- Return the number of pair orbits.
		 *
		 * @return the number of pair orbits.
		 */
		int numOrbits();

		/**
		 * getOrbit -- Get the (-1)-terminated array containing the
		 * 	ith pair orbit.
		 *
		 * @param i from 0 to numOrbits()-1, the index of the orbit.
		 * @return an array of pair indices, terminated by -1.
		 */
		int* getOrbit(int i);

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
		int getOrbitIndexForEdge(int i, int j);


		/**
		 * computeStabilizer
		 *
		 * Compute and store the automorphism group when a single pair is stabilized.
		 *
		 * @param i one endpoint
		 * @param j another endpoint
		 */
		void computeStabilizer(int i, int j);

		/**
		 * numStablizedOrbits
		 *
		 * @param i one endpoint
		 * @param j another endpoint
		 * @return the number of stabilized completion-orbits for the pair {i,j}
		 */
		int numStabilizedOrbits(int i, int j);

		/**
		 * getStabilizedCompletion
		 *
		 * Get a representative completion for the kth completion-orbit
		 * 	within the stabilizer of {i,j}
		 *
		 * @param i one endpoint
		 * @param j another endpoint
		 * @param k the orbit index
		 * @return a representative of
	if ( this->augmentation_succeeded == false )
	{
		return;
	}the kth stabilized completion-orbit for the pair {i,j}
		 * 	which is NOT to be free()'d by the caller.
		 */
		int* getStabilizedCompletion(int i, int j, int k);

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
		bool isCanonical();


		/****************** ACCESSOR METHODS ***************************/

		/**
		 * getN() -- Get the order of the current graph.
		 *
		 * @return the number of vertices in the current graph.
		 */
		int getN();

		/**
		 * getMaxN() -- Get the largest possible order
		 *
		 * @return N
		 */
		int getMaxN();

		/**
		 * indexOf -- Get the index of the i-j edge in the adjacency matrix.
		 *
		 * @param i one of the vertices
		 * @param j another vertex
		 * @return a position in the adjacency matrix for the i-j edge. If i==j or either is negative, returns -1.
		 */
		int indexOf(int i, int j);

		/**
		 * indexToPair -- Get the pair {i,j} corresponding to the given co-lex index.
		 *
		 * @param index the co-lex index of a pair.
		 * @param i the first vertex.
		 * @param j the second vertex.
		 */
		void indexToPair(int index, int& i, int& j);

		/**
		 * isUniquelySaturated
		 *
		 * @return True iff the current graph (minus a 2-dominating vertex)
		 * 	is fully assigned and uniquely K_r-saturated.
		 */
		bool isUniquelySaturated();

		/**
		 * isFeasible
		 *
		 * @return True iff there is no copy of K_r and there is no
		 * 		non-edge with more than one completion.
		 */
		bool isFeasible();

		/**
		 * getString
		 *
		 * Get a sparse6 representation of the graph.
		 *
		 * @return a buffer containing the string. It will be free()'d by the caller.
		 */
		char* getString();


		/**
		 * getDepth
		 *
		 * Return the number of augmentations.
		 */
		int getDepth();


		/**
		 * getSmallG
		 */
		sparsegraph* getSmallG();


		/**
		 * getAdjacency
		 */
		char getAdjacency(int index);
<<<<<<< HEAD

=======
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be
};

#endif /* SATURATIONGRAPH_HPP_ */
