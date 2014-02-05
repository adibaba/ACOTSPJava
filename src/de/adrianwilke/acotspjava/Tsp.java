package de.adrianwilke.acotspjava;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
import de.adrianwilke.acotspjava.InOut.Distance_type;

public class Tsp {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: TSP.c
     * Author: Thomas Stuetzle
     * Purpose: TSP related procedures, distance computation, neighbour lists
     * Check: README and gpl.txt
     * Copyright (C) 2002 Thomas Stuetzle
     */

    /***************************************************************************
     * Program's name: acotsp
     * 
     * Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
     * symmetric TSP
     * 
     * Copyright (C) 2004 Thomas Stuetzle
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or
     * (at your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     * GNU General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
     * 
     * email: stuetzle no@spam informatik.tu-darmstadt.de
     * mail address: Universitaet Darmstadt
     * Fachbereich Informatik
     * Hochschulstr. 10
     * D-64283 Darmstadt
     * Germany
     ***************************************************************************/

    static class point {
	double x;
	double y;
    }

    static class problem {
	String name; /* instance name */
	String edge_weight_type; /* selfexplanatory */
	int optimum; /* optimal tour length if known, otherwise a bound */
	int n; /* number of cities */
	int n_near; /* number of nearest neighbors */
	point[] nodeptr; /* array of structs containing coordinates of nodes */
	int[][] distance; /* distance matrix: distance[i][j] gives distance */
	int[][] nn_list; /* nearest neighbor list; contains for each node i a sorted list of n_near nearest neighbors */
    }

    static int n; /* number of cities in the instance to be solved */

    static problem instance;

    static double dtrunc(double x) {
	int k;

	k = (int) x;
	x = (double) k;
	return x;
    }

    // static int (*distance)(int , int ); /* function pointer */

    /*
     * FUNCTION: the following four functions implement different ways of
     * computing distances for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     */

    static int round_distance(int i, int j)
    /*
     * FUNCTION: compute Euclidean distances between two nodes rounded to next
     * integer for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     * COMMENTS: for the definition of how to compute this distance see TSPLIB
     */
    {
	double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
	double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
	double r = Math.sqrt(xd * xd + yd * yd) + 0.5;

	return (int) r;
    }

    static int ceil_distance(int i, int j)
    /*
     * FUNCTION: compute ceiling distance between two nodes rounded to next
     * integer for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     * COMMENTS: for the definition of how to compute this distance see TSPLIB
     */
    {
	double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
	double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
	double r = Math.sqrt(xd * xd + yd * yd);

	return (int) Math.ceil(r);
    }

    static int geo_distance(int i, int j)
    /*
     * FUNCTION: compute geometric distance between two nodes rounded to next
     * integer for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     * COMMENTS: adapted from concorde code
     * for the definition of how to compute this distance see TSPLIB
     */
    {
	double deg, min;
	double lati, latj, longi, longj;
	double q1, q2, q3;
	int dd;
	double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x, y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

	deg = dtrunc(x1);
	min = x1 - deg;
	lati = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
	deg = dtrunc(x2);
	min = x2 - deg;
	latj = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;

	deg = dtrunc(y1);
	min = y1 - deg;
	longi = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
	deg = dtrunc(y2);
	min = y2 - deg;
	longj = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;

	q1 = Math.cos(longi - longj);
	q2 = Math.cos(lati - latj);
	q3 = Math.cos(lati + latj);
	dd = (int) (6378.388 * Math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
	return dd;

    }

    static int att_distance(int i, int j)
    /*
     * FUNCTION: compute ATT distance between two nodes rounded to next
     * integer for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     * COMMENTS: for the definition of how to compute this distance see TSPLIB
     */
    {
	double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
	double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
	double rij = Math.sqrt((xd * xd + yd * yd) / 10.0);
	double tij = dtrunc(rij);
	int dij;

	if (tij < rij)
	    dij = (int) tij + 1;
	else
	    dij = (int) tij;
	return dij;
    }

    static int[][] compute_distances()
    /*
     * FUNCTION: computes the matrix of all intercity distances
     * INPUT: none
     * OUTPUT: pointer to distance matrix, has to be freed when program stops
     */
    {
	int i, j;
	int matrix[][] = new int[n][n];
	for (i = 0; i < n; i++) {
	    for (j = 0; j < n; j++) {
		if (InOut.distance_type == Distance_type.ATT) {
		    matrix[i][j] = att_distance(i, j);
		} else if (InOut.distance_type == Distance_type.CEIL_2D) {
		    matrix[i][j] = ceil_distance(i, j);
		} else if (InOut.distance_type == Distance_type.EUC_2D) {
		    matrix[i][j] = round_distance(i, j);
		} else if (InOut.distance_type == Distance_type.GEO) {
		    matrix[i][j] = geo_distance(i, j);
		}
	    }
	}
	return matrix;
    }

    static int[][] compute_nn_lists()
    /*
     * FUNCTION: computes nearest neighbor lists of depth nn for each city
     * INPUT: none
     * OUTPUT: pointer to the nearest neighbor lists
     */
    {
	int i, node, nn;
	int[] distance_vector = new int[n];
	int[] help_vector = new int[n];

	// TRACE ( System.out.println("\n computing nearest neighbor lists, "); )

	nn = Math.max(LocalSearch.nn_ls, Ants.nn_ants);
	if (nn >= n)
	    nn = n - 1;

	int[][] m_nnear = new int[n][nn];

	// DEBUG ( assert( n > nn ); )

	// TRACE ( System.out.println("nn = %ld ... \n",nn); )

	for (node = 0; node < n; node++) { /* compute cnd-sets for all node */

	    for (i = 0; i < n; i++) { /* Copy distances from nodes to the others */
		distance_vector[i] = instance.distance[node][i];
		help_vector[i] = i;
	    }
	    distance_vector[node] = Integer.MAX_VALUE; /* city is not nearest neighbour */
	    Utilities.sort2(distance_vector, help_vector, 0, n - 1);
	    for (i = 0; i < nn; i++) {
		m_nnear[node][i] = help_vector[i];
	    }
	}

	// TRACE ( System.out.println("\n    .. done\n"); )
	return m_nnear;
    }

    static int compute_tour_length(int[] t)
    /*
     * FUNCTION: compute the tour length of tour t
     * INPUT: pointer to tour t
     * OUTPUT: tour length of tour t
     */
    {
	int i;
	int tour_length = 0;

	for (i = 0; i < n; i++) {
	    tour_length += instance.distance[t[i]][t[i + 1]];
	}
	return tour_length;
    }

    static boolean tsp_check_tour(int[] t) {
	boolean error = false;

	int i;
	int[] used = new int[n];
	int size = n;

	if (t == null) {
	    System.err.println("error: permutation is not initialized!");
	    System.exit(1);
	}

	for (i = 0; i < size; i++) {
	    if (used[t[i]] != 0) {
		System.err.println("error: solution vector has two times the value " + t[i] + "(last position: " + i
			+ ")");
		error = true;
	    } else
		used[t[i]] = 1; // was true
	}

	if (!error)
	    for (i = 0; i < size; i++) {
		if (used[i] == 0) {
		    System.out.println("error: vector position " + i + " not occupied");
		    error = true;
		}
	    }

	if (!error)
	    if (t[0] != t[size]) {
		System.err.println("error: permutation is not a closed tour.");
		error = true;
	    }

	if (!error)
	    return true;

	// error: (ansi c mark)

	System.err.println("error: solution_vector:");
	for (i = 0; i < size; i++)
	    System.err.println(t[i]);
	System.out.println();

	return false;
    }
}
