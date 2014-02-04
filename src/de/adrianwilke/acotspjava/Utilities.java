package de.adrianwilke.acotspjava;

import java.util.Random;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class Utilities {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: utilities.c
     * Author: Thomas Stuetzle
     * Purpose: some additional useful procedures
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

    public static final int MAXIMUM_NO_TRIES = 100;

    private static Random random;

    static int seed;

    static double mean(int[] values, int max)
    /*
     * FUNCTION: compute the average value of an integer array of length max
     * INPUT: pointer to array, length of array
     * OUTPUT: average
     * (SIDE)EFFECTS: none
     */
    {
	int j;
	double m;

	m = 0.;
	for (j = 0; j < max; j++) {
	    m += (double) values[j];
	}
	m = m / (double) max;
	return m;
    }

    static double meanr(double[] values, int max)
    /*
     * FUNCTION: compute the average value of a floating number array of length max
     * INPUT: pointer to array, length of array
     * OUTPUT: average
     * (SIDE)EFFECTS: none
     */
    {
	int j;
	double m;

	m = 0.;
	for (j = 0; j < max; j++) {
	    m += values[j];
	}
	m = m / (double) max;
	return m;
    }

    static double std_deviation(int[] values, int max, double mean)
    /*
     * FUNCTION: compute the standard deviation of an integer array
     * INPUT: pointer to array, length of array, mean
     * OUTPUT: standard deviation
     * (SIDE)EFFECTS: none
     */
    {
	int j;
	double dev = 0.;

	if (max <= 1)
	    return 0.;
	for (j = 0; j < max; j++) {
	    dev += ((double) values[j] - mean) * ((double) values[j] - mean);
	}
	return Math.sqrt(dev / (double) (max - 1));
    }

    static double std_deviationr(double[] values, int max, double mean)
    /*
     * FUNCTION: compute the standard deviation of a floating number array
     * INPUT: pointer to array, length of array, mean
     * OUTPUT: standard deviation
     * (SIDE)EFFECTS: none
     */
    {
	int j;
	double dev;

	if (max <= 1)
	    return 0.;
	dev = 0.;
	for (j = 0; j < max; j++) {
	    dev += ((double) values[j] - mean) * ((double) values[j] - mean);
	}
	return Math.sqrt(dev / (double) (max - 1));
    }

    static int best_of_vector(int[] values, int l)
    /*
     * FUNCTION: return the minimum value in an integer value
     * INPUT: pointer to array, length of array
     * OUTPUT: smallest number in the array
     * (SIDE)EFFECTS: none
     */
    {
	int min, k;

	k = 0;
	min = values[k];
	for (k = 1; k < l; k++) {
	    if (values[k] < min) {
		min = values[k];
	    }
	}
	return min;
    }

    static int worst_of_vector(int[] values, int l)
    /*
     * FUNCTION: return the maximum value in an integer value
     * INPUT: pointer to array, length of array
     * OUTPUT: largest number in the array
     * (SIDE)EFFECTS: none
     */
    {
	int max, k;

	k = 0;
	max = values[k];
	for (k = 1; k < l; k++) {
	    if (values[k] > max) {
		max = values[k];
	    }
	}
	return max;
    }

    static double quantil(int v[], double q, int l)
    /*
     * FUNCTION: return the q-quantil of an ordered integer array
     * INPUT: one array, desired quantil q, length of array
     * OUTPUT: q-quantil of array
     * (SIDE)EFFECTS: none
     */
    {
	int i, j;
	double tmp;

	tmp = q * (double) l;
	if ((double) ((int) tmp) == tmp) {
	    i = (int) tmp;
	    j = (int) (tmp + 1.);
	    return ((double) v[i - 1] + (double) v[j - 1]) / 2.;
	} else {
	    i = (int) (tmp + 1.);
	    return v[i - 1];
	}
    }

    static void swap(int v[], int i, int j)
    /*
     * FUNCTION: auxiliary routine for sorting an integer array
     * INPUT: array, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of array are swapped
     */
    {
	int tmp;

	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;
    }

    static void sort(int v[], int left, int right)
    /*
     * FUNCTION: recursive routine (quicksort) for sorting an array
     * INPUT: one array, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
	int k, last;

	if (left >= right)
	    return;
	swap(v, left, (left + right) / 2);
	last = left;
	for (k = left + 1; k <= right; k++)
	    if (v[k] < v[left])
		swap(v, ++last, k);
	swap(v, left, last);
	sort(v, left, last);
	sort(v, last + 1, right);
    }

    static void swap2(int v[], int v2[], int i, int j)
    /*
     * FUNCTION: auxiliary routine for sorting an integer array
     * INPUT: two arraya, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
	int tmp;

	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;
	tmp = v2[i];
	v2[i] = v2[j];
	v2[j] = tmp;
    }

    static void sort2(int v[], int v2[], int left, int right)
    /*
     * FUNCTION: recursive routine (quicksort) for sorting one array; second
     * arrays does the same sequence of swaps
     * INPUT: two arrays, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
	int k, last;

	if (left >= right)
	    return;
	swap2(v, v2, left, (left + right) / 2);
	last = left;
	for (k = left + 1; k <= right; k++)
	    if (v[k] < v[left])
		swap2(v, v2, ++last, k);
	swap2(v, v2, left, last);
	sort2(v, v2, left, last);
	sort2(v, v2, last + 1, right);
    }

    static double ran01(long idum)
    /*
     * FUNCTION: generate a random number that is uniformly distributed in [0,1]
     * INPUT: pointer to variable with the current seed
     * OUTPUT: random number uniformly distributed in [0,1]
     * (SIDE)EFFECTS: random number seed is modified (important, this has to be done!)
     * ORIGIN: numerical recipes in C
     */
    {
	if (random == null) {
	    random = new Random(seed);
	}

	return random.nextDouble();
    }

    static int random_number(long idum)
    /*
     * FUNCTION: generate an integer random number
     * INPUT: pointer to variable containing random number seed
     * OUTPUT: integer random number uniformly distributed in {0,2147483647}
     * (SIDE)EFFECTS: random number seed is modified (important, has to be done!)
     * ORIGIN: numerical recipes in C
     */
    {
	if (random == null) {
	    random = new Random(seed);
	}

	return random.nextInt(2147483647);
    }

    static int[][] generate_int_matrix(int n, int m)
    /*
     * FUNCTION: malloc a matrix and return pointer to it
     * INPUT: size of matrix as n x m
     * OUTPUT: pointer to matrix
     * (SIDE)EFFECTS:
     */
    {
	return new int[n][m];
    }

    static double[][] generate_double_matrix(int n, int m)
    /*
     * FUNCTION: malloc a matrix and return pointer to it
     * INPUT: size of matrix as n x m
     * OUTPUT: pointer to matrix
     * (SIDE)EFFECTS:
     */
    {
	return new double[n][m];
    }

    static int aw_best_tour_index() {
	int min, k;

	int[] values = InOut.best_in_try;
	int l = InOut.max_tries;

	k = 0;
	min = values[k];
	int minIndex = 0;
	for (k = 1; k < l; k++) {
	    if (values[k] < min) {
		min = values[k];
		minIndex = k;
	    }
	}
	return minIndex;
    }
}
