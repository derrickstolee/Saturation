/*
 * SaturationData.hpp
 *
 *  Created on: Apr 14, 2011
 *      Author: stolee
 */

#ifndef SATURATIONDATA_HPP_
#define SATURATIONDATA_HPP_

#include <stack>
#include "nausparse.h"
#include "Set.hpp"
#include "Augmentation.hpp"

//#define _USE_BACKUPS_


/* options for type of augmentations/deletions */

/**
 * This option specifies that all vertices should be added to the graph
 * 	before augmenting.
 */
#define _AUGMENT_ALL_VERTICES_FIRST_ true

/**
 * This option specifies that the zero edges are augmented
 * as a connected graph the entire way.
 *
 * This conflicts with _DELETE_AUTOMATIC_
 */
//#define _DELETE_ALWAYS_CONNECTED_

/**
 * This specifies that we should delete an edge which has
 * 	an automatic completion.
 */
#define _DELETE_AUTOMATIC_ true

/**
 * This specifies that we check 0-edges up to orbit before trying augmentations.
 */
#define _AUGMENT_0_EDGE_CHECK_ false

/**
 * This specifies that we do degree checks before canonical checks.
 */
#define _DO_DEGREE_CHECK_ true


#define _CANONICAL_RESTRICT_CONNECTIVITY 0
#define _CANONICAL_RESTRICT_AUTOMATIC 0
#define _CANONICAL_RESTRICT_ZERO_MAX_DEGREES 1
#define _CANONICAL_RESTRICT_ZERO_MULT_DEGREES 2
#define _CANONICAL_RESTRICT_ZEROCANON 3
#define _CANONICAL_RESTRICT_ONE_MAX_DEGREES 4
#define _CANONICAL_RESTRICT_ONE_MULT_DEGREES 5
#define _CANONICAL_RESTRICT_LAYERCANON 6

/**
 * an edgetuple stores the info necessary for
 * 	changing an edge from one type to another.
 */
typedef struct tuple_struct
{
		int i;
		int j;
		int index;
		char orig_type;
		char to_type;
		bool inCompletion;
} edgetuple;

/**
 * SaturationData is a second attempt at storing data
 * for the Saturation project. It stores all graph information
 * as well as data structures which find complete subgraphs and
 * multiple K_r completions.
 *
 * It behaves in a stack-like fashion, storing a stack of
 * 	operations which can be reversed.
 * However, these operations are stored edge-by-edge.
 *
 * The 'snapshot()' method stores the current data point as
 * 	if it was a stack entry. Any edge operations after this
 * point are stored in the stack.
 *
 * The 'rollback()' method reverses the edge operations
 * 	until reaching the position of the top 'snapshot()' call.
 */
class SaturationData
{
	protected:
		/**
		 * maxN the largest number of vertices possible.
		 */
		int maxN;

		/**
		 * n the current number of vertices.
		 */
		int n;

		/**
		 * r -- Specify which K_r we are using for unique saturation
		 */
		int r;

		/**
		 * edgestack -- The stack of edge augmentations.
		 */
		std::stack<edgetuple> edgestack;

		/**
		 * depth -- The current snapshot depth.
		 */
		int depth;

		/**
		 * maxDepth -- The maximum snapshot size.
		 *
		 * Will increase with enough snapshots.
		 */
		int maxDepth;

		/**
		 * snapshots -- The number of edge augmentations at the current depth.
		 */
		size_t* snapshots;

		/**
		 * snapshotorder -- The number of vertices at each snapshot depth.
		 */
		int* snapshotorder;

		/**
		 * augmentorder -- The number of augmentations at
		 * each snapshot.
		 */
		int* augmentorder;

		int initSnapshots();
		void testSnapshotSize();
		void cleanSnapshots();

		/**
		 * augmentations -- An array of augmentations for
		 * 	storing orbis.
		 */
		Augmentation** augmentations;
		int numAugmentations;
		int augmentationsSize;
		int initAugmentations();
		void testAugmentationSize();
		void cleanAugmentations();
		void addAugmentation( Augmentation* augment );
		void removeAugmentation();

#ifdef _USE_BACKUPS_
		/* These are depth-backups for each data type */
		/* Used for debugging */
		char** bu_adjmat;
		void copyAdjacencyMatrix();
		bool compareAdjacencyMatrix();

		char** bu_degseq;
		void copyDegreeSequences();
		bool compareDegreeSequences();

		char** bu_adjlist;
		void copyAdjacencyLists();
		bool compareAdjacencyLists();

		Set*** bu_set_common_neighs;
		void copyCommonNeighbors();
		bool compareCommonNeighbors();

		char** bu_edges_in_set;
		void copyEdgesInSet();
		bool compareEdgesInSet();

		char** bu_edges_to_set;
		void copyEdgesToSet();
		bool compareEdgesToSet();

		int** bu_completion_sets;
		void copyCompletionSets();
		bool compareCompletionSets();

		Set** bu_two_edges;
		void copyTwoEdges();
		bool compareTwoEdges();
#endif

		/**
		 * layer_graph is the two-layered graph which places edges in each layer
		 * 	for type 0 and type 1 edges.
		 * It is built to grow as edges and vertices are added.
		 * The layers alternate: normal vertex 0 maps to 0 and 1, vert 1 is 2 and 3, etc.
		 *
		 * The graph is modified by the following operations:
		 *  addEdgeToLayerGraph( i,  j,  type)
		 *  removeEdgeFromLayerGraph(i, j, type)
		 *  addVertexToLayerGraph()
		 *  removeVertexFromLayerGraph()
		 */
		sparsegraph* layer_graph;
		int initLayerGraph();
		void cleanLayerGraph();
		void addEdgeToLayerGraph( int i, int j, char type );
		void removeEdgeFromLayerGraph( int i, int j, char type );
		void addVertexToLayerGraph();
		void removeVertexFromLayerGraph();

		/**
		 * one_graph is the standard graph of just the 1-type edges.
		 *
		 * The graph is built to grow with vertex and edge additions.
		 *
		 * The graph is modified by the methods
		 * 	addEdgeToOneGraph(i,j)
		 *  removeEdgeFromOneGraph(i,j)
		 */
		sparsegraph* one_graph;
		int initOneGraph();
		void cleanOneGraph();
		void addEdgeToOneGraph( int i, int j );
		void removeEdgeFromOneGraph( int i, int j );
		void addVertexToOneGraph();
		void removeVertexFromOneGraph();

		/**
		 * zero_graph is the graph containing only the 0-type edges.
		 *
		 * The graph is built to grow with vertex and edge additions.
		 *
		 * The graph is modified by the methods
		 * 	addEdgeToZeroGraph(i,j)
		 *  removeEdgeFromZeroGraph(i,j)
		 */
		sparsegraph* zero_graph;
		int initZeroGraph();
		void cleanZeroGraph();
		void addEdgeToZeroGraph( int i, int j );
		void removeEdgeFromZeroGraph( int i, int j );
		void addVertexToZeroGraph();
		void removeVertexFromZeroGraph();

		/**
		 * The adjacency matrix stores edge types of 0, 1, 2, and 3.
		 *
		 * 0: A known non-edge. It MUST have a completion.
		 * 1: A known edge.
		 * 2: An unknown edge.
		 * 3: An edge pair which uses a vertex that is not currently in the graph.
		 */
		char* adjmat;
		int initAdjacencyMatrix();
		void cleanAdjacencyMatrix();

		/**
		 * The degree sequence stores two types of degrees: 0 and 1.
		 *
		 * It is a 2-dimensional array, maxN x 2, packed into a single character list.
		 *
		 * The data type is char since we will never reach close to even 50 vertices.
		 *
		 * The getDegree(i,type) and setDegree(i,type) methods modify this sequence.
		 */
		char* degseq;
		int initDegreeSequence();
		void cleanDegreeSequence();
		void setDegree( int i, int type, char value );

		/**
		 * the adjacency lists store the (0/1)-neighborhoods of each vertex,
		 * 	in order of augmentations.
		 *
		 * It is a 3-dimensional array, maxN x 2 x maxN, packed into a single list
		 *
		 * The methods addNeighbor(i,j,type) and removeNeighbor(i,j,type) modify
		 * 	the lists, the degrees, AND the adjacency matrix.
		 */
		char* adjlist;
		int getBasePosition( int i, char type );
		int initAdjacencyList();
		void cleanAdjacencyList();
		bool addNeighbor( int i, int j, char type );
		void removeNeighbor( int i, int j, char type );

		/**
		 * set_common_neighs stores sets of common neighbors for each set
		 * of order s (s from 2 to r-2)
		 *
		 * set_common_base_index stores the base index for each size.
		 *
		 * So, the array is (N choose 2) + (N choose 3) + ... + (N choose r-2)
		 */
		Set** set_common_neighs;
		int* set_common_base_index;
		int initSetCommonNeighs();
		void cleanSetCommonNeighs();
		int getMaxCommonNeighsSetSize();
		Set* getSetCommonNeighs( int size, int index );
		bool addSetCommonNeigh( int size, int set_index, int i );
		void removeSetCommonNeigh( int size, int set_index, int i );

		/**
		 * edges_in_set stores the number of edges within a set of order 2...(r-1)
		 *
		 * This is modified by the methods addEdgeInSet(size, index) and
		 * 	removeEdgeInSet(size, index).
		 * Returns FALSE if a K_r is created.
		 */
		char* edges_in_set;
		int* edges_in_set_base_index;
		int initEdgesInSet();
		void cleanEdgesInSet();
		int getMaxEdgesInSetSize();
		bool addEdgeInSet( int size, int index );
		void removeEdgeInSet( int size, int index );
		bool isSetFull( int size, int index );

		/**
		 * edges_to_set stores the number of edges between a set of order 2...(r-1)
		 * 	and another vertex i
		 */
		char* edges_to_set;
		int* edges_to_set_base_index;
		int initEdgesToSet();
		void cleanEdgesToSet();
		int getEdgesToSetPosition( int size, int index, int i );
		int getMaxEdgesToSetSize();
		bool addEdgeToSet( int size, int index, int i );
		void removeEdgeToSet( int size, int index, int i );
		bool doesSetDominate( int size, int index, int i );

		/**
		 * completion_sets stores the index of the (r-2)-set which
		 * 	completes the given edge. Value is -1 if there is no completion.
		 */
		int* completion_sets;
		int initCompletionSets();
		void cleanCompletionSets();
		int getCompletionSet( int i, int j );
		int getCompletionSet( int index );
		bool setCompletionSet( int i, int j, int completionIndex );
		void removeCompletionSet( int i, int j, int completionIndex );

		/**
		 * completion_mult stores the number of times a 1-edge is in
		 * 	a completion.
		 */
		char* completion_mult;
		int initCompletionMult();
		void cleanCompletionMult();
		char getCompletionMult( int index );
		bool addCompletionMult( int index );
		void removeCompletionMult( int index );

		/**
		 * two_edges
		 */
		Set* two_edges;
		Set* zero_edges;
		int initTwoEdges();
		void cleanTwoEdges();
		bool addTwoEdge( int index );
		void removeTwoEdge( int index );

		/**
		 * removeEdge()
		 *
		 * Remove an edge from the graph using the following edge change.
		 */
		void removeEdge( int i, int j, int index, char to_type, char orig_type, bool inCompletion );

		/**
		 * removeVertex()
		 *
		 * Remove the nth vertex.
		 */
		void removeVertex();

	public:
		/**
		 * Constructor
		 */
		SaturationData( int r, int maxN );

		/**
		 * Destructor
		 */
		virtual ~SaturationData();

		/**
		 * snapshot()
		 */
		void snapshot();

		/**
		 * rollback()
		 */
		void rollback();

		/**
		 * augmentEdge()
		 *
		 * Augment the graph using the following edge change.
		 *
		 * return FALSE if the constraints are violated.
		 */
		bool augmentEdge( int i, int j, char type, bool inCompletion );

		/**
		 * augmentVertex()
		 *
		 * Add a new vertex.
		 *
		 * return FALSE if there are too many vertices.
		 */
		bool augmentVertex();

		/**
		 * augmentFill(index)
		 *
		 * Augment by performing a fill operation on the given
		 * 	orbit (by index).
		 * If index < 0, then fill ALL edges!
		 */
		bool augmentFill( int index );

		/**
		 * augmentCompletion()
		 *
		 * Augment by making a 0-edge with a completion.
		 */
		bool augmentCompletion( int i, int j, int completionIndex );

		/**
		 * augmentAutoCompletions
		 *
		 * Augment by performing all auto-completions.
		 */
		bool augmentAutoCompletions();

		/**
		 * has2TypeWithCompletion
		 *
		 * Does the current graph have edges of 2-type
		 * 	with automatic completions?
		 */
		bool has2TypeWithCompletion();

		/**
		 * getLayeredGraph()
		 */
		sparsegraph* getLayeredGraph();

		/**
		 * getOneGraph()
		 */
		sparsegraph* getOneGraph();

		/**
		 * getZeroGraph()
		 */
		sparsegraph* getZeroGraph();

		/**
		 * getOneGraphString()
		 *
		 * NOTE: DO NOT FREE() THIS STRING!
		 */
		char* getOneGraphString();

		/**
		 * isUniquelyKrSaturated
		 */
		bool isUniquelyKrSaturated();

		/**
		 * get2Orbits(numOrbits)
		 *
		 * Returns an array of orbit representatives and
		 * 	places the number of orbits in the numOrbits parameter.
		 */
		int* get2Orbits( int& numOrbits );

		/**
		 * get2OrbitByIndex(int index)
		 *
		 * Returns all 2-edges which are in	the
		 * 	given orbit (by orbit order)
		 */
		int* get2OrbitByIndex( int index );

		/**
		 * get2OrbitIndexByRepresentative(int rep)
		 *
		 * Returns the orbit index for the 2-edges which are in	the
		 * 	given orbit (by an orbit representative)
		 */
		int get2OrbitIndexByRepresentative( int rep );
		/**
		 * get2OrbitByRepresentative(int rep)
		 *
		 * Returns all 2-edges which are in	the
		 * 	given orbit (by an orbit representative)
		 */
		int* get2OrbitByRepresentative( int rep );

		/**
		 * getCompletionOrbits(index, numCompleters)
		 *
		 * Returns an array of completion indices corresponding
		 * 	to possible completions of the given index.
		 * The number of completion orbits is placed in numCompleters.
		 */
		int* getCompletionOrbits( int index, int& numCompleters );

		/**
		 * Is the most recent set of augmentations canonical?
		 */
		bool isCanonical();

		/**
		 * isZeroCanonical()
		 *
		 * Is the given 2-edge going to be canonical
		 * 	in the zero graph if augmented to be zero?
		 */
		bool isZeroCanonical( int index );

		/**
		 * restrictCanonical()
		 *
		 * Given a set of potential 0-edges for deletion, restrict
		 * 	the set using the given level of specificity.
		 */
		void restrictCanonical( Set* potential, Augmentation* augment, int rank );

		/**
		 * canonicalType()
		 *
		 * Given the current graph, what kind of augmentation
		 * 	did we expect?
		 */
		AUGMENT_TYPE canonicalType();

		/**
		 * getN Return the current number of vertices.
		 */
		int getN();

		/**
		 * getMaxN Return the maximum number of vertices.
		 */
		int getMaxN();

		/**
		 * getAdjacency return the adjacency for this pair
		 * index or vertex pair.
		 */
		char getAdjacency( int index );
		char getAdjacency( int i, int j );

		/**
		 * getNum2Edges
		 */
		int getNum2Edges();


		char getDegree( int i, int type );

		/**
		 * printEdgeStack
		 */
		void printEdgeStack();
};

#endif /* SATURATIONDATA_HPP_ */
