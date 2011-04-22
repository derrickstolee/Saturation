/*
 * TreeSet.cpp
 *
 *  Created on: Feb 3, 2010
 *      Author: derrickstolee
 *
 *  Oct 13, 2010: Refactored from Set.cpp to new object TreeSet.cpp.
 *
 *     Implements the bnode structure to be a B-tree as described in
 *     Chapter 18: B-Trees of Introduction to Algorithms by CLRS
 *
 *     Further, implements the class Set, which contains integer elements
 *     within a B-tree.
 */

#include "Set.hpp"
#include "TreeSet.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * Constructor.
 */
TreeSet::TreeSet()
{
	this->count = 0;

	/* initialize the B tree */
	this->root = binit();

	this->array = 0;
	this->array_up_to_date = false;
}

/**
 * Constructor.
 *
 * @param T - the minimum degree
 */
TreeSet::TreeSet(int T)
{
	this->t = T;

	this->count = 0;

	this->root = binit();

	this->array = 0;
	this->array_up_to_date = false;
}

/**
 * Copy Constructor
 */
TreeSet::TreeSet(TreeSet& s)
{
	this->t = s.t;
	this->count = 0;

	/* initialize the B tree */
	this->root = binit();

	/* naively add the elements */
	int copies = 0;
	for ( s.resetIterator(); s.hasNext(); )
	{
		this->add(s.next());
		copies++;
	}

	this->array = 0;
	this->array_up_to_date = false;
}

/**
 * Destructor
 */
TreeSet::~TreeSet()
{
	if ( this->root != 0 )
	{
		/*  iterate through and delete all elements */
		bfree(this->root);
		this->root = 0;
	}

	if ( this->array != 0 )
	{
		free(this->array);
	}
}

/**
 * size - return the number of elements in the set.
 *
 * @return the number of elements
 */
int TreeSet::size()
{
	return this->count;
}

/**
 * binit - initialize the root of a tree
 */
bnode_ptr binit()
{
	bnode_ptr node = (bnode_ptr) malloc(sizeof(bnode));

	node->leaf = 1;
	node->key_count = 0;

	node->keys = (int*) malloc(NUM_B_KEYS * sizeof(int));
	node->children = (bnode_ptr*) malloc(NUM_B_CHILDREN * sizeof(bnode_ptr));

	bzero(node->keys, NUM_B_KEYS * sizeof(int));
	bzero(node->children, NUM_B_CHILDREN * sizeof(bnode_ptr));

	return node;
}

/**
 * bverify - Check that a node is valid.
 */
int bverify(bnode_ptr z)
{
	int i;
	int res = 0;

	if ( z->key_count < 0 || z->key_count >= NUM_B_KEYS )
	{
		printf("[bverify] we have a bad number of keys! %d.\n", z->key_count);
		res = 1;
	}

	if ( res == 0 )
	{
		for ( i = 0; i < z->key_count; i++ )
		{
			if ( z->keys[i] < 0 )
			{
				printf("[bverify] we have a bad key!\n");
				res = 1;
			}

			if ( (z->leaf == 0) && (z->children[i] == 0 || z->children[i + 1]
					== 0) )
			{
				printf("[bverify] we have bad children! %d %lX %lX\n", i,
						(long int) (z->children[i]), (long int) (z->children[i
								+ 1]));
				res = 1;
			}
			else if ( z->leaf == 0 )
			{
				int max = bmax(z->children[i]);
				int min = bmin(z->children[i + 1]);

				if ( max >= z->keys[i] )
				{
					printf("[bverify] we have a bad order: %d >= %d\n", max,
							z->keys[i]);
					res = 1;
				}
				if ( min <= z->keys[i] )
				{
					printf("[bverify] we have a bad order: %d <= %d\n", min,
							z->keys[i]);
					res = 1;
				}
			}
		}
	}

	if ( res == 1 )
	{
		printf("[bverify] the bad node: %d keys - ", z->key_count);

		for ( i = 0; i < z->key_count; i++ )
		{
			printf("%d ", z->keys[i]);
		}
		printf("\n");

		return 1;
	}

	return 0;
}

int bdeepverify(bnode_ptr node)
{
	int res = bverify(node);

	if ( res == 1 )
	{
		return 1;
	}

	if ( node->leaf == 0 )
	{
		for ( int i = 0; i <= node->key_count; i++ )
		{
			res = bdeepverify(node->children[i]);

			if ( res == 1 )
			{
				return 1;
			}
		}
	}

	return 0;
}

/**
 * bfree - free up the memory from a node, including recursing below.
 */
void bfree(bnode_ptr node)
{
	if ( node->leaf == 0 )
	{
		int i;
		for ( i = 0; i <= node->key_count; i++ )
		{
			/* recurse on the child */
			bfree(node->children[i]);

			/* each child frees themselves, so clear out data */
			node->children[i] = 0;
		}
	}

	free(node->keys);
	free(node->children);

	node->key_count = -1;
	node->leaf = -1;
	free(node);
}

/**
 * bsearch - recursively search the tree for an element.
 *
 * @param elt the element
 * @return elt if it is there, -1 otherwise
 */
int bsearch(bnode_ptr node, int elt)
{
	int i = 0;

	for ( i = 0; i < node->key_count && elt > node->keys[i]; i++ )
		;

	if ( i < node->key_count && elt == node->keys[i] )
	{
		return elt;
	}

	if ( node->leaf == 0 )
	{
		return bsearch(node->children[i], elt);
	}
	else
	{
		return -1;
	}
}

/**
 * binsert - insert an element into the btree.
 *
 * @return the new root of the tree
 */
bnode_ptr binsert(bnode_ptr node, int elt)
{
	if ( node->key_count == 2 * BTREE_T - 1 )
	{
		/* create a new root! */
		bnode_ptr s = binit();
		s->leaf = 0;
		s->key_count = 0; /* no between-keys */
		s->children[0] = node; /* one child */
		bsplit(s, 0); /* now, do the splits! */
		binsert(s, elt);

		// if (  DO_VERIFY bverify(s); #endif

		/* return the new root */
		return s;
	}
	else
	{
		int i = node->key_count - 1;

		if ( node->leaf == 1 )
		{
			/* if it is a leaf, stick it here! */
			/* look for key position */
			while ( i >= 0 && elt < node->keys[i] )
			{
				node->keys[i + 1] = node->keys[i];
				i--;
			}

			/* overshot! */
			i++;

			node->keys[i] = elt;
			node->key_count = node->key_count + 1;

			// if (  DO_VERIFY bverify(node); #endif

			/* return this leaf! */
			return node;
		}
		else
		{
			/* it is not a leaf */
			/* look for step down */
			while ( i >= 0 && elt < node->keys[i] )
			{
				i--;
			}

			/* overshot */
			i++;

			/* split, or recurse? */
			if ( node->children[i]->key_count == (2 * BTREE_T) - 1 )
			{
				bsplit(node, i);
				// if (  DO_VERIFY bverify(node); #endif

				if ( elt > node->keys[i] )
				{
					i++;
				}

				binsert(node->children[i], elt);

				// if (  DO_VERIFY bverify(node->children[i]); #endif
				// if (  DO_VERIFY bverify(node); #endif
				/* return the regular node */
				return node;
			}

			binsert(node->children[i], elt);

			// if (  DO_VERIFY bverify(node->children[i]); #endif
			// if (  DO_VERIFY bverify(node); #endif

			return node;
		}
	}
}

/**
 * bsplit - split a node into parts at position 'pos'
 */
bnode_ptr bsplit(bnode_ptr node, int pos)
{
	bnode_ptr z = binit();
	bnode_ptr y = node->children[pos];

	/* initialize z data */
	z->leaf = y->leaf;
	z->key_count = BTREE_T - 1;

	int j;

	/* put the keys from y into z */
	for ( j = 0; j < BTREE_T - 1; j++ )
	{
		z->keys[j] = y->keys[j + BTREE_T];
		//		y->keys[j + BTREE_T] = -1;
	}

	if ( y->leaf == 0 )
	{
		/* copy children, too! */
		for ( j = 0; j < BTREE_T; j++ )
		{
			z->children[j] = y->children[j + BTREE_T];
		}
	}

	/* insert at 'pos' the z value */
	for ( j = node->key_count; j > pos; j-- )
	{
		node->children[j + 1] = node->children[j];
	}
	node->children[pos + 1] = z;

	/* shuffle the keys */
	for ( j = node->key_count - 1; j >= pos; j-- )
	{
		node->keys[j + 1] = node->keys[j];
	}
	/* finally, place the splitting key into new position */
	node->keys[pos] = y->keys[BTREE_T - 1];

	node->key_count = node->key_count + 1;

	y->keys[BTREE_T - 1] = -1;
	y->key_count = BTREE_T - 1;

	// if (  DO_VERIFY bverify(node); #endif
	// if (  DO_VERIFY bverify(y); #endif
	// if (  DO_VERIFY bverify(z); #endif
	bverify(node);

	return node;
}

/**
 * bremove - remove an element 'elt' from the node
 *
 * @return the new root, if needed
 */
bnode_ptr bremove(bnode_ptr node, int elt)
{
	// if (  DO_VERIFY bverify(node); #endif
	int i = 0, j = 0;

	/* find the key position */
	while ( i < node->key_count && node->keys[i] < elt )
	{
		i++;
	}
	/* now, either node->keys[i] == elt, OR elt is in node->children[i] */

	if ( node->leaf == 1 )
	{
		/* CASE 1 */
		/* delete the key from the leaf! */
		if ( node->keys[i] == elt )
		{
			/* delete it! */
			node->key_count = (node->key_count) - 1;
			for ( j = i; j < node->key_count; j++ )
			{
				node->keys[j] = node->keys[j + 1];
			}
		}

		// if (  DO_VERIFY bverify(node); #endif

		return node;
	}
	else if ( i < node->key_count && node->keys[i] == elt )
	{
		/* CASE 2 */
		/* delete the key from this internal node! */
		if ( node->children[i]->key_count >= BTREE_T )
		{
			/* CASE 2.a. */
			/* if the predecessor has room, take that key! */
			bnode_ptr prechild = node->children[i];

			int key = bmax(prechild);

			/* replace elt by k */
			node->keys[i] = key;

			bremove(prechild, key);

			// if (  DO_VERIFY bverify(prechild); #endif
			// if (  DO_VERIFY bverify(node); #endif

			return node;
		}
		else if ( node->children[i + 1]->key_count >= BTREE_T )
		{
			/* CASE 2.b. */
			/* if the successor has room, take that key! */
			bnode_ptr postchild = node->children[i + 1];

			int key = bmin(postchild);

			bremove(postchild, key);

			/* replace elt by k */
			node->keys[i] = key;

			// if (  DO_VERIFY bverify(postchild); #endif
			// if (  DO_VERIFY bverify(node); #endif

			return node;
		}
		else
		{
			/* CASE 2.c. */
			/* both are small enough for merging! */
			bnode_ptr prechild = node->children[i];
			bnode_ptr postchild = node->children[i + 1];

			/* prechild keys + k + postchild keys */
			/* then delete k */
			prechild->keys[prechild->key_count] = elt;

			int j;
			for ( j = 0; j < postchild->key_count; j++ )
			{
				prechild->keys[prechild->key_count + j + 1]
						= postchild->keys[j];
				prechild->children[prechild->key_count + j + 1]
						= postchild->children[j];
				postchild->children[j] = 0;
			}
			prechild->children[prechild->key_count + j + 1]
					= postchild->children[j];
			postchild->children[j] = 0;

			prechild->key_count = prechild->key_count + postchild->key_count
					+ 1;

			/* merging is complete, delete postchild */
			postchild->leaf = 1;
			bfree(postchild);

			/* shrink the node */
			for ( j = i; j < node->key_count - 1; j++ )
			{
				node->keys[j] = node->keys[j + 1];
				node->children[j + 1] = node->children[j + 2];
			}
			node->children[j + 1] = 0;
			node->key_count = node->key_count - 1;

			/* also, delete elt from prechild */
			bremove(prechild, elt);

			// if (  DO_VERIFY bverify(prechild); #endif
			// if (  DO_VERIFY bverify(node); #endif

			return node;
		}
	}
	else
	{
		/* CASE 3. */
		/* delete the key from a child! */
		/* check for siblings with extra! */
		bnode_ptr child = node->children[i];
		bnode_ptr presibling = 0;
		bnode_ptr postsibling = 0;

		if ( i > 0 )
		{
			presibling = node->children[i - 1];
		}
		if ( i < node->key_count )
		{
			postsibling = node->children[i + 1];
		}

		/* should be in child i */
		while ( child->key_count < BTREE_T )
		{
			if ( presibling != 0 && presibling->key_count >= BTREE_T )
			{
				/* CASE 3.a. */
				/* shuffle! pre -> key -> child */

				/* make room in child */
				child->children[child->key_count + 2]
						= child->children[child->key_count + 1];
				for ( j = child->key_count; j >= 0; j-- )
				{
					child->keys[j + 1] = child->keys[j];
					child->children[j + 1] = child->children[j];
				}
				child->key_count = child->key_count + 1;

				/* rotate keys and children */
				int key = node->keys[i - 1];
				node->keys[i - 1] = presibling->keys[presibling->key_count - 1];
				child->keys[0] = key;

				child->children[0]
						= presibling->children[presibling->key_count];

				/* clear spot in presibling */
				presibling->keys[presibling->key_count - 1] = -1;
				presibling->children[presibling->key_count] = 0;

				/* lower count in presibling */
				presibling->key_count = presibling->key_count - 1;

				// if (  DO_VERIFY bverify(presibling); #endif
				// if (  DO_VERIFY bverify(child); #endif
				// if (  DO_VERIFY bverify(node); #endif
			}
			else if ( postsibling != 0 && postsibling->key_count >= BTREE_T )
			{
				/* CASE 3.a. */
				/* shuffle! post -> key -> child */

				/* rotate keys and child */
				int key = node->keys[i];
				child->keys[child->key_count] = key;
				child->children[child->key_count + 1]
						= postsibling->children[0];
				node->keys[i] = postsibling->keys[0];
				child->key_count = child->key_count + 1;

				postsibling->key_count = postsibling->key_count - 1;
				for ( j = 0; j < postsibling->key_count; j++ )
				{
					postsibling->keys[j] = postsibling->keys[j + 1];
					postsibling->children[j] = postsibling->children[j + 1];
				}
				postsibling->children[j] = postsibling->children[j + 1];

				/* clear out postsibling */
				postsibling->keys[j] = -1;
				postsibling->children[j + 1] = 0;

				// if (  DO_VERIFY bverify(postsibling); #endif
				// if (  DO_VERIFY bverify(child); #endif
				// if (  DO_VERIFY bverify(node); #endif
			}
			else
			{
				/* CASE 3.b. */
				/* MERGE child with a sibling... GROSS! */
				if ( presibling != 0 )
				{
					/* presibling + key + child */
					int key = node->keys[i - 1];
					presibling->keys[presibling->key_count] = key;

					presibling->children[presibling->key_count + 1]
							= child->children[0];
					for ( j = 0; j < child->key_count; j++ )
					{
						presibling->keys[presibling->key_count + j + 1]
								= child->keys[j];
						presibling->children[presibling->key_count + j + 2]
								= child->children[j + 1];
					}

					/* new key count */
					presibling->key_count = presibling->key_count
							+ child->key_count + 1;

					/* delete child */
					child->leaf = 1;
					bfree(child);
					child = 0;

					/* child = presibling */
					child = presibling;
					// if (  DO_VERIFY bverify(child); #endif

					/* time to reduce the node's data */
					node->key_count = node->key_count - 1;
					i--;
					for ( j = i; j < node->key_count; j++ )
					{
						node->keys[j] = node->keys[j + 1];
						node->children[j + 1] = node->children[j + 2];
					}

					node->children[j + 1] = 0;
					node->keys[j] = -1;

					if ( i > 0 )
					{
						presibling = node->children[i - 1];
					}
					else
					{
						presibling = 0;
					}

					// if (  DO_VERIFY bverify(node); #endif
				}
				else if ( postsibling != 0 )
				{
					/* child + key + postsibling */
					int key = node->keys[i];
					child->keys[child->key_count] = key;

					child->children[child->key_count + 1]
							= postsibling->children[0];
					for ( j = 0; j < postsibling->key_count; j++ )
					{
						child->keys[child->key_count + j + 1]
								= postsibling->keys[j];
						child->children[child->key_count + j + 2]
								= postsibling->children[j + 1];
					}

					/* new key count */
					child->key_count = child->key_count
							+ postsibling->key_count + 1;
					// if (  DO_VERIFY bverify(child); #endif

					/* delete child */
					postsibling->leaf = 1;
					bfree(postsibling);
					postsibling = 0;

					/* time to reduce the node's data */
					node->key_count = node->key_count - 1;
					for ( j = i; j < node->key_count; j++ )
					{
						node->keys[j] = node->keys[j + 1];
						node->children[j + 1] = node->children[j + 2];
					}

					node->keys[j] = -1;
					node->children[j + 1] = 0;

					if ( i < node->key_count - 1 )
					{
						postsibling = node->children[i + 1];
					}

					// if (  DO_VERIFY bverify(node); #endif
				} /* End postsibling != 0 */
			} /* END CASE 3.B */
		}

		/* recurse! */
		bnode_ptr root = bremove(child, elt);

		// if (  DO_VERIFY bverify(child); #endif
		// if (  DO_VERIFY bverify(node); #endif

		if ( node->key_count == 0 )
		{
			/* popping off the top! */
			node->leaf = 1;
			bfree(node);
			return root;
		}

		return node;
	}
}

/**
 * bnext - the element in the set strictly larger than this element.
 *
 * @param node - the node to recurse over
 * @param elt - the element to be strictly larger
 * @return the next element in the list.  -1 if no more.
 */
int bnext(bnode_ptr node, int elt)
{
	/* first, find the key that is just larger than this element, if one exists */
	int i;
	for ( i = 0; i < node->key_count; i++ )
	{
		//	printf("%d ", node->keys[i]);
		if ( elt == node->keys[i] )
		{
			if ( node->leaf == 1 )
			{
				if ( i < node->key_count - 1 )
				{
					return node->keys[i + 1];
				}

				return -1;
			}
			else
			{
				return bmin(node->children[i + 1]);
			}
		}
		else if ( elt < node->keys[i] )
		{
			/* we stop here, and recurse! */
			int min = node->keys[i];

			if ( node->leaf == 0 )
			{
				int other = bnext(node->children[i], elt);
				if ( other >= 0 && other < min )
				{
					return other;
				}
			}

			return min;
		}
	}

	/* if it is not a leaf, try the right-most child! */
	if ( node->leaf == 0 )
	{
		return bnext(node->children[node->key_count], elt);
	}

	/* else, */
	return -1;
}

/**
 * bmax - get the maximum element
 */
int bmax(bnode_ptr node)
{
	if ( node->leaf == 0 )
	{
		int cur_max = bmax(node->children[node->key_count]);

		int i = 1;
		while ( cur_max < 0 && i <= node->key_count )
		{
			cur_max = bmax(node->children[node->key_count - i]);
			i++;
		}

		return cur_max;
	}

	if ( node->key_count == 0 )
	{
		return -1;
	}
	return node->keys[node->key_count - 1];
}

/**
 * bmin - get the minimum element
 */
int bmin(bnode_ptr node)
{
	if ( node->leaf == 0 )
	{
		return bmin(node->children[0]);
	}

	return node->keys[0];
}

/**
 * bprint - output everything about this node and its children
 */
void bprint(bnode_ptr node)
{
	printf("\t-> %2d keys:\t", node->key_count);
	int i;
	for ( i = 0; i < node->key_count && i < BTREE_T * 2; i++ )
	{
		printf(" [%d] ", node->keys[i]);
	}

	if ( node->leaf == 0 )
	{
		printf("\n");
		for ( i = 0; i <= node->key_count && i < BTREE_T * 2 + 1; i++ )
		{
			bprint(node->children[i]);
		}

		printf("<- %2d keys:\t", node->key_count);
		for ( i = 0; i < node->key_count && i < BTREE_T * 2; i++ )
		{
			printf(" [%d] ", node->keys[i]);
		}

		printf("\n");
	}
	else
	{
		printf(" <- \n");
	}
}

/**
 * contains - checks if an element is in the set
 *
 * @param elt the element tho find.
 * @return 1 if elt is in the Set. 0 otherwise.
 */
int TreeSet::contains(int elt)
{
	if ( bsearch(this->root, elt) == elt )
	{
		return 1;
	}

	return 0;
}

/**
 * remove - removes an element from the set
 *
 * @param elt the element to find
 * @return elt if the element was found and removed. -1 otherwise.
 */
int TreeSet::remove(int elt)
{
	// if (  DO_VERIFY bdeepverify(this->root); #endif

	if ( bsearch(this->root, elt) == elt )
	{
		this->root = bremove(this->root, elt);

		/* this may be wrong! */
		this->count = this->count - 1;
		this->array_up_to_date = false;
		// if (  DO_VERIFY bverify(this->root);  #endif
	}

	return elt;
}

/**
 * add - adds an element to the set
 *
 * @param elt the element to add
 * @return elt
 */
int TreeSet::add(int elt)
{
	if ( elt < 0 )
	{
		printf("--[TreeSet] Trying to insert %d < 0!\n", elt);
		return elt;
	}

	if ( bsearch(this->root, elt) != elt )
	{
		this->root = binsert(this->root, elt);
		// if (  DO_VERIFY bverify(this->root);  #endif
		this->count = this->count + 1;
		this->array_up_to_date = false;
	}

	return elt;
}

/**
 * join - form a new Set as the union of this and another Set.
 *
 * @param set the other set
 * @return a new Set filled with elements from both sets
 */
Set* TreeSet::join(Set* set)
{
	Set* u = new TreeSet();

	for ( this->resetIterator(); this->hasNext(); )
	{
		u->add(this->next());
	}

	for ( set->resetIterator(); set->hasNext(); )
	{
		u->add(set->next());
	}

	return u;
}

/**
 * intersect - form a new Set as the intersection of this and another Set.
 *
 * @param set the other set
 * @return a new Set filled with elements in both sets
 */
Set* TreeSet::intersect(Set* set)
{
	Set* u = new TreeSet();

	if ( this->size() < set->size() )
	{
		for ( this->resetIterator(); this->hasNext(); )
		{
			int elt = this->next();

			if ( set->contains(elt) )
			{
				u->add(elt);
			}
		}
	}
	else
	{
		for ( set->resetIterator(); set->hasNext(); )
		{
			int elt = set->next();

			if ( this->contains(elt) )
			{
				u->add(elt);
			}
		}
	}

	return u;
}

/**
 * resetIterator - reset the iteration through the set.
 * That is, we want to enumerate the elements in the set, but
 * need to start from the beginning after this call.
 */
void TreeSet::resetIterator()
{
	this->cur_elt = -1;
}

/**
 * next - get the next element in the Set.
 *
 * @return the next element in the traversal
 */
int TreeSet::next()
{
	this->cur_elt = bnext(this->root, this->cur_elt);

	return this->cur_elt;
}

/**
 * hasNext - there is another element in the Set that we have not
 * enumerated yet.
 *
 * @return 1 if there is, 0 if there is not.
 */
int TreeSet::hasNext()
{
	if ( bnext(this->root, this->cur_elt) >= 0 )
	{
		return 1;
	}

	return 0;
}

/**
 * min -- get the minimum key.
 */
int TreeSet::min()
{
	return bmin(this->root);
}

/**
 * max -- get the maximum key.
 */
int TreeSet::max()
{
	return bmax(this->root);
}

/**
 * print
 */
void TreeSet::print()
{
	bprint(this->root);
}

/**
 * toArray -- get the set as an array that the set will clean up.
 */
int* TreeSet::toArray()
{
	if ( this->array_up_to_date == true )
	{
		return this->array;
	}

	this->array = (int*) realloc(this->array, this->size() * sizeof(int));

	int i = 0;
	for ( this->resetIterator(); this->hasNext(); i++ )
	{
		this->array[i] = this->next();
	}

	this->array_up_to_date = true;

	return this->array;
}

/**
 * clear -- remove all elements from the set.
 */
void TreeSet::clear()
{
	for ( this->resetIterator(); this->hasNext(); )
	{
		int i = this->next();
		this->remove(i);
	}
}

