/*
 * translation.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: stolee
 */

#ifndef TRANSLATION_HPP_
#define TRANSLATION_HPP_

/**
 * indexOf
 *
 * The co-lex index of the pair (i,j).
 */
int indexOf( int i, int j );

/**
 * indexToPair
 */
void indexToPair( int index, int& i, int& j );

/**
 * indexOfSet
 *
 * Get the co-lex order of the set of a given size.
 */
int indexOfSet( int size, int* set );

/**
 * indexToSet
 */
void indexToSet( int size, int index, int* set );

/**
 * sortSet
 */
void sortSet( int size, int* set );

/**
 * nChooseK
 */
int nChooseK( int n, int k );

/**
 * indexOfTuple
 *
 * Get the product order of the set of a given size.
 */
<<<<<<< HEAD
int indexOfTuple( int n, int size, int* tuple );
=======
int indexOfTuple(int n, int size, int* tuple);
>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

/**
 * indexToTuple
 */
<<<<<<< HEAD
void indexToTuple( int n, int size, int index, int* tuple );

/**
 * numSetsW
 *
 * Get the number of sets with i
 */
int numSetsW( int n, int s, int i );

/**
 * numSetsWW
 *
 * Get the number of sets of size s
 * 	within [n] with BOTH entries
 */
int numSetsWW( int n, int s, int i, int j );

/**
 * numSetsWWO
 *
 * Get the number of sets with first entry, but NOT second.
 */
int numSetsWWO( int n, int s, int i, int j );

/**
 * numSetsWOWO
 *
 * Get the number of sets without either entry.
 */
int numSetsWOWO( int n, int s, int i, int j );

/**
 * getSetW
 */
void getSetW( int s, int i, int index, int* set );

/**
 * getSetWW
 */
void getSetWW( int s, int i, int j, int index, int* set );

/**
 * getSetWWO
 */
void getSetWWO( int s, int i, int j, int index, int* set );

/**
 * getSetWOWO
 */
void getSetWOWO( int s, int i, int j, int index, int* set );
=======
void indexToTuple(int n, int size, int index, int* tuple);

>>>>>>> 54e39aa63d051762fc09f22b848b5047ef0e72be

#endif /* TRANSLATION_HPP_ */
