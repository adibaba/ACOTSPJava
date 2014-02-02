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
public class Ants {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: ants.c
     * Author: Thomas Stuetzle
     * Purpose: implementation of procedures for ants' behaviour
     * Check: README.txt and legal.txt
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

    static class ant_struct {
	int[] tour;
	boolean[] visited;
	int tour_length;
    }

    public static final int MAX_ANTS = 1024;
    public static final int MAX_NEIGHBOURS = 512;

    static ant_struct ant[];
    static ant_struct best_so_far_ant;
    static ant_struct restart_best_ant;

    static double pheromone[][];
    static double total[][];

    static double prob_of_selection[];

    static int n_ants; /* number of ants */
    static int nn_ants; /*
			 * length of nearest neighbor lists for the ants'
			 * solution construction
			 */

    static double rho; /* parameter for evaporation */
    static double alpha; /* importance of trail */
    static double beta; /* importance of heuristic evaluate */
    static double q_0; /* probability of best choice in tour construction */

    static boolean as_flag; /* ant system */
    static boolean eas_flag; /* elitist ant system */
    static boolean ras_flag; /* rank-based version of ant system */
    static boolean mmas_flag; /* MAX-MIN ant system */
    static boolean bwas_flag; /* best-worst ant system */
    static boolean acs_flag; /* ant colony system */

    static int elitist_ants; /*
			      * additional parameter for elitist
			      * ant system, no. elitist ants
			      */

    static int ras_ranks; /*
			   * additional parameter for rank-based version
			   * of ant system
			   */

    static double trail_max; /* maximum pheromone trail in MMAS */
    static double trail_min; /* minimum pheromone trail in MMAS */
    static int u_gb; /* every u_gb iterations update with best-so-far ant */

    static double trail_0; /* initial pheromone level in ACS and BWAS */

    static double HEURISTIC(int m, int n) {
	return (1.0 / ((double) Tsp.instance.distance[m][n] + 0.1));
    }

    public static final double EPSILON = 0.00000000000000000000000000000001;

    static void allocate_ants()
    /*
     * FUNCTION: allocate the memory for the ant colony, the best-so-far and
     * the iteration best ant
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: allocation of memory for the ant colony and two ants that
     * store intermediate tours
     */
    {
	int i;

	ant = new ant_struct[n_ants];

	for (i = 0; i < n_ants; i++) {
	    ant[i] = new ant_struct();
	    ant[i].tour = new int[Tsp.n + 1];
	    ant[i].visited = new boolean[Tsp.n];
	}
	best_so_far_ant = new ant_struct();

	best_so_far_ant.tour = new int[Tsp.n + 1];
	best_so_far_ant.visited = new boolean[Tsp.n];

	restart_best_ant = new ant_struct();

	restart_best_ant.tour = new int[Tsp.n + 1];
	restart_best_ant.visited = new boolean[Tsp.n];

	prob_of_selection = new double[nn_ants + 1];
	for (i = 0; i < nn_ants + 1; i++) {
	    prob_of_selection[i] = Double.POSITIVE_INFINITY;
	}
    }

    static int find_best()
    /*
     * FUNCTION: find the best ant of the current iteration
     * INPUT: none
     * OUTPUT: index of struct containing the iteration best ant
     * (SIDE)EFFECTS: none
     */
    {
	int min;
	int k, k_min;

	min = ant[0].tour_length;
	k_min = 0;
	for (k = 1; k < n_ants; k++) {
	    if (ant[k].tour_length < min) {
		min = ant[k].tour_length;
		k_min = k;
	    }
	}
	return k_min;
    }

    static int find_worst()
    /*
     * FUNCTION: find the worst ant of the current iteration
     * INPUT: none
     * OUTPUT: pointer to struct containing iteration best ant
     * (SIDE)EFFECTS: none
     */
    {
	int max;
	int k, k_max;

	max = ant[0].tour_length;
	k_max = 0;
	for (k = 1; k < n_ants; k++) {
	    if (ant[k].tour_length > max) {
		max = ant[k].tour_length;
		k_max = k;
	    }
	}
	return k_max;
    }

    /************************************************************
     ************************************************************ 
     Procedures for pheromone manipulation
     ************************************************************ 
     ************************************************************/

    static void init_pheromone_trails(double initial_trail)
    /*
     * FUNCTION: initialize pheromone trails
     * INPUT: initial value of pheromone trails "initial_trail"
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromone matrix is reinitialized
     */
    {
	int i, j;

	// TRACE ( System.out.println(" init trails with %.15f\n",initial_trail); );

	/* Initialize pheromone trails */
	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j <= i; j++) {
		pheromone[i][j] = initial_trail;
		pheromone[j][i] = initial_trail;
		total[i][j] = initial_trail;
		total[j][i] = initial_trail;
	    }
	}
    }

    static void evaporation()
    /*
     * FUNCTION: implements the pheromone trail evaporation
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are reduced by factor rho
     */
    {
	int i, j;

	// TRACE ( System.out.println("pheromone evaporation\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j <= i; j++) {
		pheromone[i][j] = (1 - rho) * pheromone[i][j];
		pheromone[j][i] = pheromone[i][j];
	    }
	}
    }

    static void evaporation_nn_list()
    /*
     * FUNCTION: simulation of the pheromone trail evaporation
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are reduced by factor rho
     * REMARKS: if local search is used, this evaporation procedure
     * only considers links between a city and those cities
     * of its candidate list
     */
    {
	int i, j, help_city;

	// TRACE ( System.out.println("pheromone evaporation nn_list\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < nn_ants; j++) {
		help_city = Tsp.instance.nn_list[i][j];
		pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
	    }
	}
    }

    static void global_update_pheromone(ant_struct a)
    /*
     * FUNCTION: reinforces edges used in ant k's solution
     * INPUT: pointer to ant that updates the pheromone trail
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
     */
    {
	int i, j, h;
	double d_tau;

	// TRACE ( System.out.println("global pheromone update\n"); );

	d_tau = 1.0 / (double) a.tour_length;
	for (i = 0; i < Tsp.n; i++) {
	    j = a.tour[i];
	    h = a.tour[i + 1];
	    pheromone[j][h] += d_tau;
	    pheromone[h][j] = pheromone[j][h];
	}
    }

    static void global_update_pheromone_weighted(ant_struct a, int weight)
    /*
     * FUNCTION: reinforces edges of the ant's tour with weight "weight"
     * INPUT: pointer to ant that updates pheromones and its weight
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones of arcs in the ant's tour are increased
     */
    {
	int i, j, h;
	double d_tau;

	// TRACE ( System.out.println("global pheromone update weighted\n"); );

	d_tau = (double) weight / (double) a.tour_length;
	for (i = 0; i < Tsp.n; i++) {
	    j = a.tour[i];
	    h = a.tour[i + 1];
	    pheromone[j][h] += d_tau;
	    pheromone[h][j] = pheromone[j][h];
	}
    }

    static void compute_total_information()
    /*
     * FUNCTION: calculates heuristic info times pheromone for each arc
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j;

	// TRACE ( System.out.println("compute total information\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < i; j++) {
		total[i][j] = Math.pow(pheromone[i][j], alpha) * Math.pow(HEURISTIC(i, j), beta);
		total[j][i] = total[i][j];
	    }
	}
    }

    static void compute_nn_list_total_information()
    /*
     * FUNCTION: calculates heuristic info times pheromone for arcs in nn_list
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j, h;

	// TRACE ( System.out.println("compute total information nn_list\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < nn_ants; j++) {
		h = Tsp.instance.nn_list[i][j];
		if (pheromone[i][h] < pheromone[h][i])
		    /* force pheromone trails to be symmetric as much as possible */
		    pheromone[h][i] = pheromone[i][h];
		total[i][h] = Math.pow(pheromone[i][h], alpha) * Math.pow(HEURISTIC(i, h), beta);
		total[h][i] = total[i][h];
	    }
	}
    }

    /****************************************************************
     **************************************************************** 
     Procedures implementing solution construction and related things
     **************************************************************** 
     ****************************************************************/

    static void ant_empty_memory(ant_struct a)
    /*
     * FUNCTION: empty the ants's memory regarding visited cities
     * INPUT: ant identifier
     * OUTPUT: none
     * (SIDE)EFFECTS: vector of visited cities is reinitialized to Utilities.FALSE
     */
    {
	int i;

	for (i = 0; i < Tsp.n; i++) {
	    a.visited[i] = false;
	}
    }

    static void place_ant(ant_struct a, int step)
    /*
     * FUNCTION: place an ant on a randomly chosen initial city
     * INPUT: pointer to ant and the number of construction steps
     * OUTPUT: none
     * (SIDE)EFFECT: ant is put on the chosen city
     */
    {
	int rnd;

	rnd = (int) (Utilities.ran01(Utilities.seed) * (double) Tsp.n); /* random number between 0 .. n-1 */

	a.tour[step] = rnd;
	a.visited[rnd] = true;
    }

    static void choose_best_next(ant_struct a, int phase)
    /*
     * FUNCTION: chooses for an ant as the next city the one with
     * maximal value of heuristic information times pheromone
     * INPUT: pointer to ant and the construction step
     * OUTPUT: none
     * (SIDE)EFFECT: ant moves to the chosen city
     */
    {
	int city, current_city, next_city;
	double value_best;

	next_city = Tsp.n;
	assert (phase > 0 && phase < Tsp.n);
	current_city = a.tour[phase - 1];
	value_best = -1.; /* values in total matrix are always >= 0.0 */
	for (city = 0; city < Tsp.n; city++) {
	    if (a.visited[city])
		; /* city already visited, do nothing */
	    else {
		if (total[current_city][city] > value_best) {
		    next_city = city;
		    value_best = total[current_city][city];
		}
	    }
	}
	assert (0 <= next_city && next_city < Tsp.n);
	assert (value_best > 0.0);
	assert (a.visited[next_city] == false);
	a.tour[phase] = next_city;
	a.visited[next_city] = true;
    }

    static void neighbour_choose_best_next(ant_struct a, int phase)
    /*
     * FUNCTION: chooses for an ant as the next city the one with
     * maximal value of heuristic information times pheromone
     * INPUT: pointer to ant and the construction step "phase"
     * OUTPUT: none
     * (SIDE)EFFECT: ant moves to the chosen city
     */
    {
	int i, current_city, next_city, help_city;
	double value_best, help;

	next_city = Tsp.n;
	assert (phase > 0 && phase < Tsp.n);
	current_city = a.tour[phase - 1];
	assert (0 <= current_city && current_city < Tsp.n);
	value_best = -1.; /* values in total matix are always >= 0.0 */
	for (i = 0; i < nn_ants; i++) {
	    help_city = Tsp.instance.nn_list[current_city][i];
	    if (a.visited[help_city])
		; /* city already visited, do nothing */
	    else {
		help = total[current_city][help_city];
		if (help > value_best) {
		    value_best = help;
		    next_city = help_city;
		}
	    }
	}
	if (next_city == Tsp.n)
	    /* all cities in nearest neighbor list were already visited */
	    choose_best_next(a, phase);
	else {
	    assert (0 <= next_city && next_city < Tsp.n);
	    assert (value_best > 0.0);
	    assert (a.visited[next_city] == false);
	    a.tour[phase] = next_city;
	    a.visited[next_city] = true;
	}
    }

    static void choose_closest_next(ant_struct a, int phase)
    /*
     * FUNCTION: Chooses for an ant the closest city as the next one
     * INPUT: pointer to ant and the construction step "phase"
     * OUTPUT: none
     * (SIDE)EFFECT: ant moves to the chosen city
     */
    {
	int city, current_city, next_city, min_distance;

	next_city = Tsp.n;
	assert (phase > 0 && phase < Tsp.n);
	current_city = a.tour[phase - 1];
	min_distance = Integer.MAX_VALUE; /* Search shortest edge */
	for (city = 0; city < Tsp.n; city++) {
	    if (a.visited[city])
		; /* city already visited */
	    else {
		if (Tsp.instance.distance[current_city][city] < min_distance) {
		    next_city = city;
		    min_distance = Tsp.instance.distance[current_city][city];
		}
	    }
	}
	assert (0 <= next_city && next_city < Tsp.n);
	a.tour[phase] = next_city;
	a.visited[next_city] = true;
    }

    static void neighbour_choose_and_move_to_next(ant_struct a, int phase)
    /*
     * FUNCTION: Choose for an ant probabilistically a next city among all
     * unvisited cities in the current city's candidate list.
     * If this is not possible, choose the closest next
     * INPUT: pointer to ant the construction step "phase"
     * OUTPUT: none
     * (SIDE)EFFECT: ant moves to the chosen city
     */
    {
	int i, help;
	int current_city;
	double rnd, partial_sum = 0., sum_prob = 0.0;
	/* double *prob_of_selection; *//*
					 * stores the selection probabilities
					 * of the nearest neighbor cities
					 */
	double prob_ptr[];

	if ((q_0 > 0.0) && (Utilities.ran01(Utilities.seed) < q_0)) {
	    /*
	     * with a probability q_0 make the best possible choice
	     * according to pheromone trails and heuristic information
	     */
	    /*
	     * we first check whether q_0 > 0.0, to avoid the very common case
	     * of q_0 = 0.0 to have to compute a random number, which is
	     * expensive computationally
	     */
	    neighbour_choose_best_next(a, phase);
	    return;
	}

	prob_ptr = prob_of_selection;

	current_city = a.tour[phase - 1]; /* current_city city of ant k */
	assert (current_city >= 0 && current_city < Tsp.n);
	for (i = 0; i < nn_ants; i++) {
	    if (a.visited[Tsp.instance.nn_list[current_city][i]])
		prob_ptr[i] = 0.0; /* city already visited */
	    else {
		assert (Tsp.instance.nn_list[current_city][i] >= 0 && Tsp.instance.nn_list[current_city][i] < Tsp.n);
		prob_ptr[i] = total[current_city][Tsp.instance.nn_list[current_city][i]];
		sum_prob += prob_ptr[i];
	    }
	}

	if (sum_prob <= 0.0) {
	    /* All cities from the candidate set are tabu */
	    choose_best_next(a, phase);
	} else {
	    /*
	     * at least one neighbor is eligible, chose one according to the
	     * selection probabilities
	     */
	    rnd = Utilities.ran01(Utilities.seed);
	    rnd *= sum_prob;
	    i = 0;
	    partial_sum = prob_ptr[i];
	    /* This loop always stops because prob_ptr[nn_ants] == HUGE_VAL */
	    while (partial_sum <= rnd) {
		i++;
		partial_sum += prob_ptr[i];
	    }
	    /*
	     * This may very rarely happen because of rounding if rnd is
	     * close to 1.
	     */
	    if (i == nn_ants) {
		neighbour_choose_best_next(a, phase);
		return;
	    }
	    assert (0 <= i && i < nn_ants);
	    assert (prob_ptr[i] >= 0.0);
	    help = Tsp.instance.nn_list[current_city][i];
	    assert (help >= 0 && help < Tsp.n);
	    assert (a.visited[help] == false);
	    a.tour[phase] = help; /* Tsp.instance.nn_list[current_city][i]; */
	    a.visited[help] = true;
	}
    }

    /****************************************************************
     **************************************************************** 
     Procedures specific to MAX-MIN Ant System
     **************************************************************** 
     ****************************************************************/

    static void mmas_evaporation_nn_list()
    /*
     * FUNCTION: simulation of the pheromone trail evaporation for MMAS
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are reduced by factor rho
     * REMARKS: if local search is used, this evaporation procedure
     * only considers links between a city and those cities
     * of its candidate list
     */
    {
	int i, j, help_city;

	// TRACE ( System.out.println("mmas specific evaporation on nn_lists\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < nn_ants; j++) {
		help_city = Tsp.instance.nn_list[i][j];
		pheromone[i][help_city] = (1 - rho) * pheromone[i][help_city];
		if (pheromone[i][help_city] < trail_min)
		    pheromone[i][help_city] = trail_min;
	    }
	}
    }

    static void check_pheromone_trail_limits()
    /*
     * FUNCTION: only for MMAS without local search:
     * keeps pheromone trails inside trail limits
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
     */
    {
	int i, j;

	// TRACE ( System.out.println("mmas specific: check pheromone trail limits\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < i; j++) {
		if (pheromone[i][j] < trail_min) {
		    pheromone[i][j] = trail_min;
		    pheromone[j][i] = trail_min;
		} else if (pheromone[i][j] > trail_max) {
		    pheromone[i][j] = trail_max;
		    pheromone[j][i] = trail_max;
		}
	    }
	}
    }

    static void check_nn_list_pheromone_trail_limits()
    /*
     * FUNCTION: only for MMAS with local search: keeps pheromone trails
     * inside trail limits
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
     * COMMENTS: currently not used since check for trail_min is integrated
     * mmas_evaporation_nn_list and typically check for trail_max
     * is not done (see FGCS paper or ACO book for explanation
     */
    {
	int i, j, help_city;

	// TRACE ( System.out.println("mmas specific: check pheromone trail limits nn_list\n"); );

	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < nn_ants; j++) {
		help_city = Tsp.instance.nn_list[i][j];
		if (pheromone[i][help_city] < trail_min)
		    pheromone[i][help_city] = trail_min;
		if (pheromone[i][help_city] > trail_max)
		    pheromone[i][help_city] = trail_max;
	    }
	}
    }

    /****************************************************************
     **************************************************************** 
     Procedures specific to Ant Colony System
     **************************************************************** 
     ****************************************************************/

    static void global_acs_pheromone_update(ant_struct a)
    /*
     * FUNCTION: reinforces the edges used in ant's solution as in ACS
     * INPUT: pointer to ant that updates the pheromone trail
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
     */
    {
	int i, j, h;
	double d_tau;

	// TRACE ( System.out.println("acs specific: global pheromone update\n"); );

	d_tau = 1.0 / (double) a.tour_length;

	for (i = 0; i < Tsp.n; i++) {
	    j = a.tour[i];
	    h = a.tour[i + 1];

	    pheromone[j][h] = (1. - rho) * pheromone[j][h] + rho * d_tau;
	    pheromone[h][j] = pheromone[j][h];

	    total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
	    total[j][h] = total[h][j];
	}
    }

    static void local_acs_pheromone_update(ant_struct a, int phase)
    /*
     * FUNCTION: removes some pheromone on edge just passed by the ant
     * INPUT: pointer to ant and number of constr. phase
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
     * COMMENTS: I did not do experiments with with different values of the parameter
     * xi for the local pheromone update; therefore, here xi is fixed to 0.1
     * as suggested by Gambardella and Dorigo for the TSP. If you wish to run
     * experiments with that parameter it may be reasonable to use it as a
     * commandline parameter
     */
    {
	int h, j;

	assert (phase > 0 && phase <= Tsp.n);
	j = a.tour[phase];

	h = a.tour[phase - 1];
	assert (0 <= j && j < Tsp.n);
	assert (0 <= h && h < Tsp.n);
	/* still additional parameter has to be introduced */
	pheromone[h][j] = (1. - 0.1) * pheromone[h][j] + 0.1 * trail_0;
	pheromone[j][h] = pheromone[h][j];
	total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
	total[j][h] = total[h][j];
    }

    /****************************************************************
     **************************************************************** 
     Procedures specific to Best-Worst Ant System
     **************************************************************** 
     ****************************************************************/

    static void bwas_worst_ant_update(ant_struct a1, ant_struct a2)
    /*
     * FUNCTION: uses additional evaporation on the arcs of iteration worst
     * ant that are not shared with the global best ant
     * INPUT: pointer to the worst (a1) and the best (a2) ant
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones on some arcs undergo additional evaporation
     */
    {
	int i, j, h, pos, pred;
	int pos2[]; /* positions of cities in tour of ant a2 */

	// TRACE ( System.out.println("bwas specific: best-worst pheromone update\n"); );

	pos2 = new int[Tsp.n];
	for (i = 0; i < Tsp.n; i++) {
	    pos2[a2.tour[i]] = i;
	}

	for (i = 0; i < Tsp.n; i++) {
	    j = a1.tour[i];
	    h = a1.tour[i + 1];
	    pos = pos2[j];
	    if (pos - 1 < 0)
		pred = Tsp.n - 1;
	    else
		pred = pos - 1;
	    if (a2.tour[pos + 1] == h)
		; /* do nothing, edge is common with a2 (best solution found so far) */
	    else if (a2.tour[pred] == h)
		; /* do nothing, edge is common with a2 (best solution found so far) */
	    else { /* edge (j,h) does not occur in ant a2 */
		pheromone[j][h] = (1 - rho) * pheromone[j][h];
		pheromone[h][j] = (1 - rho) * pheromone[h][j];
	    }
	}
    }

    static void bwas_pheromone_mutation()
    /*
     * FUNCTION: implements the pheromone mutation in Best-Worst Ant System
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j, k;
	int num_mutations;
	double avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

	// TRACE ( System.out.println("bwas specific: pheromone mutation\n"); );

	/* compute average pheromone trail on edges of global best solution */
	for (i = 0; i < Tsp.n; i++) {
	    avg_trail += pheromone[best_so_far_ant.tour[i]][best_so_far_ant.tour[i + 1]];
	}
	avg_trail /= (double) Tsp.n;

	/* determine mutation strength of pheromone matrix */
	/*
	 * FIXME: we add a small value to the denominator to avoid any
	 * potential division by zero. This may not be fully correct
	 * according to the original BWAS.
	 */
	if (InOut.max_time > 0.1)
	    mutation_strength = 4. * avg_trail * (Timer.elapsed_time() - InOut.restart_time)
		    / (InOut.max_time - InOut.restart_time + 0.0001);
	else if (InOut.max_tours > 100)
	    mutation_strength = 4. * avg_trail * (InOut.iteration - InOut.restart_iteration)
		    / (InOut.max_tours - InOut.restart_iteration + 1);
	else
	    System.out.println("apparently no termination condition applied!!\n");

	/* finally use fast version of matrix mutation */
	mutation_rate = mutation_rate / Tsp.n * nn_ants;
	num_mutations = (int) (Tsp.n * mutation_rate / 2);
	/* / 2 because of adjustment for symmetry of pheromone trails */

	if (InOut.restart_iteration < 2)
	    num_mutations = 0;

	for (i = 0; i < num_mutations; i++) {
	    j = (int) (Utilities.ran01(Utilities.seed) * (double) Tsp.n);
	    k = (int) (Utilities.ran01(Utilities.seed) * (double) Tsp.n);
	    if (Utilities.ran01(Utilities.seed) < 0.5) {
		pheromone[j][k] += mutation_strength;
		pheromone[k][j] = pheromone[j][k];
	    } else {
		pheromone[j][k] -= mutation_strength;
		if (pheromone[j][k] <= 0.0) {
		    pheromone[j][k] = EPSILON;
		}
		pheromone[k][j] = pheromone[j][k];
	    }
	}
    }

    /**************************************************************************
     ************************************************************************** 
     Procedures specific to the ant's tour manipulation other than construction
     *************************************************************************** 
     **************************************************************************/

    static void copy_from_to(ant_struct a1, ant_struct a2) {
	/*
	 * FUNCTION: copy solution from ant a1 into ant a2
	 * INPUT: pointers to the two ants a1 and a2
	 * OUTPUT: none
	 * (SIDE)EFFECTS: a2 is copy of a1
	 */
	int i;

	a2.tour_length = a1.tour_length;
	for (i = 0; i < Tsp.n; i++) {
	    a2.tour[i] = a1.tour[i];
	}
	a2.tour[Tsp.n] = a2.tour[0];
    }

    static int nn_tour()
    /*
     * FUNCTION: generate some nearest neighbor tour and compute tour length
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: needs ant colony and one statistic ants
     */
    {
	int phase, help;

	ant_empty_memory(ant[0]);

	phase = 0; /* counter of the construction steps */
	place_ant(ant[0], phase);

	while (phase < Tsp.n - 1) {
	    phase++;
	    choose_closest_next(ant[0], phase);
	}
	phase = Tsp.n;
	ant[0].tour[Tsp.n] = ant[0].tour[0];
	if (LocalSearch.ls_flag != 0) {
	    LocalSearch.two_opt_first(ant[0].tour);
	}
	InOut.n_tours += 1;
	/* copy_from_to( &ant[0], best_so_far_ant ); */
	ant[0].tour_length = Tsp.compute_tour_length(ant[0].tour);

	help = ant[0].tour_length;
	ant_empty_memory(ant[0]);
	return help;
    }

    static int distance_between_ants(ant_struct a1, ant_struct a2)
    /*
     * FUNCTION: compute the distance between the tours of ant a1 and a2
     * INPUT: pointers to the two ants a1 and a2
     * OUTPUT: distance between ant a1 and a2
     */
    {
	int i, j, h, pos, pred;
	int distance;
	int[] pos2; /* positions of cities in tour of ant a2 */

	pos2 = new int[Tsp.n];
	for (i = 0; i < Tsp.n; i++) {
	    pos2[a2.tour[i]] = i;
	}

	distance = 0;
	for (i = 0; i < Tsp.n; i++) {
	    j = a1.tour[i];
	    h = a1.tour[i + 1];
	    pos = pos2[j];
	    if (pos - 1 < 0)
		pred = Tsp.n - 1;
	    else
		pred = pos - 1;
	    if (a2.tour[pos + 1] == h)
		; /* do nothing, edge is common with best solution found so far */
	    else if (a2.tour[pred] == h)
		; /* do nothing, edge is common with best solution found so far */
	    else { /* edge (j,h) does not occur in ant a2 */
		distance++;
	    }
	}
	return distance;
    }

}
