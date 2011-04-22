/*
 * SaturationManager.hpp
 *
 *  Created on: Apr 6, 2011
 *      Author: stolee
 */

#ifndef SATURATIONMANAGER_HPP_
#define SATURATIONMANAGER_HPP_

#include "SearchManager.hpp"

#include "SaturationGraph.hpp"

/**
 * SaturationManager
 *
 * This class performs the search for uniquely K_r-saturated graphs.
 *
 * It extends the SearchManager class from the TreeSearch library.
 */
class SaturationManager : public SearchManager
{
	protected:
		/**
		 * r the saturation parameter.
		 */
		int r;

		/**
		 * N the maximum order of a graph in the search.
		 */
		int N;

		/**
		 * satgraph The saturation graph containing most of the data
		 * necessary for the search.
		 */
		SaturationGraph* satgraph;

		/**
		 * augmentation_suceeded True iff the most recent augmentation suceeded.
		 */
		bool augmentation_suceeded;


	public:
		/**
		 * Constructor
		 *
		 * @param r the 'saturation parameter' which determines
		 * 	which K_r is used for unique saturation.
		 * @param N the maximum number of vertices to use.
		 */
		SaturationManager(int r, int N);

		/**
		 * Destructor
		 */
		virtual ~SaturationManager();

		/**
		 * pushNext -- deepen the search to the next child
		 * 	of the current node.
		 *
		 * @return the label for the new node. -1 if none.
		 */
		virtual LONG_T pushNext();

		/**
		 * pushTo -- deepen the search to the specified child
		 * 	of the current node.
		 *
		 * @param child the specified label for the new node
		 * @return the label for the new node. -1 if none, or failed.
		 */
		virtual LONG_T pushTo(LONG_T child);

		/**
		 * pop -- remove the current node and move up the tree.
		 *
		 * @return the label of the node after the pop.
		 * 		This return value is used for validation purposes
		 * 		to check proper implementation of push*() and pop().
		 */
		virtual LONG_T pop();

		/**
		 * prune -- Perform a check to see if this node should be pruned.
		 *
		 * @return 0 if no prune should occur, 1 if prune should occur.
		 */
		virtual int prune();

		/**
		 * isSolution -- Perform a check to see if a solution exists
		 * 		at this point.
		 *
		 * @return 0 if no solution is found, 1 if a solution is found.
		 */
		virtual int isSolution();

		/**
		 * writeSolution -- create a buffer that contains a
		 * 	description of the solution.
		 *
		 * Solution strings start with S.
		 *
		 * @return a string of data allocated with malloc().
		 * 	It will be deleted using free().
		 */
		virtual char* writeSolution();

		/**
		 * writeStatistics -- create a buffer that contains a
		 * 	description of the solution.
		 *
		 * Statistics take the following format in each line:
		 *
		 * T [TYPE] [ID] [VALUE]
		 *
		 * @return a string of data allocated with malloc().
		 * 	It will be deleted using free().
		 */
//		virtual char* writeStatistics();
};

#endif /* SATURATIONMANAGER_HPP_ */
