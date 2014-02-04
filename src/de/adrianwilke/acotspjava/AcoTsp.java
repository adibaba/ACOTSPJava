package de.adrianwilke.acotspjava;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class AcoTsp {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: main.c
     * Author: Thomas Stuetzle
     * Purpose: main routines and control for the ACO algorithms
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

    static boolean termination_condition()
    /*
     * FUNCTION: checks whether termination condition is met
     * INPUT: none
     * OUTPUT: 0 if condition is not met, number neq 0 otherwise
     * (SIDE)EFFECTS: none
     */
    {
	return (((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time)) || (Ants.best_so_far_ant.tour_length <= InOut.optimal));
    }

    static void construct_solutions()
    /*
     * FUNCTION: manage the solution construction phase
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: when finished, all ants of the colony have constructed a solution
     */
    {
	int k; /* counter variable */
	int step; /* counter of the number of construction steps */

	// TRACE ( System.out.println("construct solutions for all ants\n"); );

	/* Mark all cities as unvisited */
	for (k = 0; k < Ants.n_ants; k++) {
	    Ants.ant_empty_memory(Ants.ant[k]);
	}

	step = 0;
	/* Place the ants on same initial city */
	for (k = 0; k < Ants.n_ants; k++)
	    Ants.place_ant(Ants.ant[k], step);

	while (step < Tsp.n - 1) {
	    step++;
	    for (k = 0; k < Ants.n_ants; k++) {
		Ants.neighbour_choose_and_move_to_next(Ants.ant[k], step);
		if (Ants.acs_flag)
		    Ants.local_acs_pheromone_update(Ants.ant[k], step);
	    }
	}

	step = Tsp.n;
	for (k = 0; k < Ants.n_ants; k++) {
	    Ants.ant[k].tour[Tsp.n] = Ants.ant[k].tour[0];
	    Ants.ant[k].tour_length = Tsp.compute_tour_length(Ants.ant[k].tour);
	    if (Ants.acs_flag)
		Ants.local_acs_pheromone_update(Ants.ant[k], step);
	}
	InOut.n_tours += Ants.n_ants;
    }

    static void init_try(int ntry)
    /*
     * FUNCTION: initilialize variables appropriately when starting a trial
     * INPUT: trial number
     * OUTPUT: none
     * COMMENTS: none
     */
    {

	// TRACE ( System.out.println("INITUtilities.IALIZE TRUtilities.IAL\n"); );

	Timer.start_timers();
	InOut.time_used = Timer.elapsed_time();
	InOut.time_passed = InOut.time_used;

	if (InOut.comp_report != null) {
	    InOut.printToFile(InOut.comp_report, "Utilities.seed " + Utilities.seed);
	}
	/* Initialize variables concerning statistics etc. */

	InOut.n_tours = 1;
	InOut.iteration = 1;
	InOut.restart_iteration = 1;
	InOut.lambda = 0.05;
	Ants.best_so_far_ant.tour_length = Integer.MAX_VALUE;
	InOut.found_best = 0;

	/*
	 * Initialize the Pheromone trails, only if ACS is used, Ants.pheromones
	 * have to be initialized differently
	 */
	if (!(Ants.acs_flag || Ants.mmas_flag || Ants.bwas_flag)) {
	    Ants.trail_0 = 1. / ((Ants.rho) * Ants.nn_tour());
	    /*
	     * in the original papers on Ant System, Elitist Ant System, and
	     * Rank-based Ant System it is not exactly defined what the
	     * initial value of the Ants.pheromones is. Here we set it to some
	     * small constant, analogously as done in MAX-MIN Ant System.
	     */
	    Ants.init_pheromone_trails(Ants.trail_0);
	}
	if (Ants.bwas_flag) {
	    Ants.trail_0 = 1. / ((double) Tsp.n * (double) Ants.nn_tour());
	    Ants.init_pheromone_trails(Ants.trail_0);
	}
	if (Ants.mmas_flag) {
	    Ants.trail_max = 1. / ((Ants.rho) * Ants.nn_tour());
	    Ants.trail_min = Ants.trail_max / (2. * Tsp.n);
	    Ants.init_pheromone_trails(Ants.trail_max);
	}
	if (Ants.acs_flag) {
	    Ants.trail_0 = 1. / ((double) Tsp.n * (double) Ants.nn_tour());
	    Ants.init_pheromone_trails(Ants.trail_0);
	}

	/* Calculate combined information Ants.pheromone times heuristic information */
	Ants.compute_total_information();

	if (InOut.comp_report != null)
	    InOut.printToFile(InOut.comp_report, "begin try " + ntry);
	if (InOut.stat_report != null)
	    InOut.printToFile(InOut.stat_report, "begin try " + ntry);
    }

    static void local_search()
    /*
     * FUNCTION: manage the local search phase; apply local search to ALL ants; in
     * dependence of LocalSearch.ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
     * is chosen.
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: all ants of the colony have locally optimal tours
     * COMMENTS: typically, best performance is obtained by applying local search
     * to all ants. It is known that some improvements (e.g. convergence
     * speed towards high quality solutions) may be obtained for some
     * ACO algorithms by applying local search to only some of the ants.
     * Overall best performance is typcially obtained by using 3-opt.
     */
    {
	int k;

	// TRACE ( System.out.println("apply local search to all ants\n"); );

	for (k = 0; k < Ants.n_ants; k++) {
	    switch (LocalSearch.ls_flag) {
	    case 1:
		LocalSearch.two_opt_first(Ants.ant[k].tour); /* 2-opt local search */
		break;
	    case 2:
		LocalSearch.two_h_opt_first(Ants.ant[k].tour); /* 2.5-opt local search */
		break;
	    case 3:
		LocalSearch.three_opt_first(Ants.ant[k].tour); /* 3-opt local search */
		break;
	    default:
		System.err.println("type of local search procedure not correctly specified");
		System.exit(1);
	    }
	    Ants.ant[k].tour_length = Tsp.compute_tour_length(Ants.ant[k].tour);
	    if (termination_condition())
		return;
	}
    }

    static void update_statistics()
    /*
     * FUNCTION: manage some statistical information about the trial, especially
     * if a new best solution (best-so-far or restart-best) is found and
     * adjust some parameters if a new best solution is found
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: restart-best and best-so-far ant may be updated; Ants.trail_min
     * and Ants.trail_max used by MMAS may be updated
     */
    {

	int iteration_best_ant;
	double p_x; /* only used by MMAS */

	iteration_best_ant = Ants.find_best(); /* iteration_best_ant is a global variable */

	if (Ants.ant[iteration_best_ant].tour_length < Ants.best_so_far_ant.tour_length) {

	    InOut.time_used = Timer.elapsed_time(); /* best sol found after time_used */
	    Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.best_so_far_ant);
	    Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.restart_best_ant);

	    InOut.found_best = InOut.iteration;
	    InOut.restart_found_best = InOut.iteration;
	    InOut.found_branching = InOut.node_branching(InOut.lambda);
	    InOut.branching_factor = InOut.found_branching;
	    if (Ants.mmas_flag) {
		if (LocalSearch.ls_flag == 0) {
		    p_x = Math.exp(Math.log(0.05) / Tsp.n);
		    Ants.trail_min = 1. * (1. - p_x) / (p_x * (double) ((Ants.nn_ants + 1) / 2));
		    Ants.trail_max = 1. / ((Ants.rho) * Ants.best_so_far_ant.tour_length);
		    Ants.trail_0 = Ants.trail_max;
		    Ants.trail_min = Ants.trail_max * Ants.trail_min;
		} else {
		    Ants.trail_max = 1. / ((Ants.rho) * Ants.best_so_far_ant.tour_length);
		    Ants.trail_min = Ants.trail_max / (2. * Tsp.n);
		    Ants.trail_0 = Ants.trail_max;
		}
	    }
	    InOut.write_report();
	}
	if (Ants.ant[iteration_best_ant].tour_length < Ants.restart_best_ant.tour_length) {
	    Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.restart_best_ant);
	    InOut.restart_found_best = InOut.iteration;
	    System.out.println("restart best: " + Ants.restart_best_ant.tour_length + " restart_found_best "
		    + InOut.restart_found_best + ", time " + Timer.elapsed_time());
	}
    }

    static void search_control_and_statistics()
    /*
     * FUNCTION: occasionally compute some statistics and check whether algorithm
     * is converged
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: restart-best and best-so-far ant may be updated; Ants.trail_min
     * and Ants.trail_max used by MMAS may be updated
     */
    {
	// TRACE ( System.out.println("SEARCH CONTROL AND STATISTICS\n"); );

	if ((InOut.iteration % 100) == 0) {
	    InOut.population_statistics();
	    InOut.branching_factor = InOut.node_branching(InOut.lambda);
	    System.out.println("best so far " + Ants.best_so_far_ant.tour_length + ", iteration: " + InOut.iteration
		    + ", time " + Timer.elapsed_time() + ", b_fac " + InOut.branching_factor);

	    if (Ants.mmas_flag && (InOut.branching_factor < InOut.branch_fac)
		    && (InOut.iteration - InOut.restart_found_best > 250)) {
		/*
		 * MAX-MIN Ant System was the first ACO algorithm to use
		 * Ants.pheromone trail re-initialisation as implemented
		 * here. Other ACO algorithms may also profit from this mechanism.
		 */
		System.out.println("INIT TRAILS!!!\n");
		Ants.restart_best_ant.tour_length = Integer.MAX_VALUE;
		Ants.init_pheromone_trails(Ants.trail_max);
		Ants.compute_total_information();
		InOut.restart_iteration = InOut.iteration;
		InOut.restart_time = Timer.elapsed_time();
	    }
	    System.out.println("try " + InOut.n_try + " iteration " + InOut.iteration + ", b-fac "
		    + InOut.branching_factor);
	}
    }

    static void as_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for Ant System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: all ants deposit Ants.pheromones on matrix "Ants.pheromone"
     */
    {
	int k;

	// TRACE ( System.out.println("Ant System Ants.pheromone deposit\n"); );

	for (k = 0; k < Ants.n_ants; k++)
	    Ants.global_update_pheromone(Ants.ant[k]);
    }

    static void eas_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for Elitist Ant System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: all ants plus elitist ant deposit Ants.pheromones on matrix "Ants.pheromone"
     */
    {
	int k;

	// TRACE ( System.out.println("Elitist Ant System Ants.pheromone deposit\n"); );

	for (k = 0; k < Ants.n_ants; k++)
	    Ants.global_update_pheromone(Ants.ant[k]);
	Ants.global_update_pheromone_weighted(Ants.best_so_far_ant, Ants.elitist_ants);
    }

    static void ras_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for Rank-based Ant System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: the Ants.ras_ranks-1 best ants plus the best-so-far ant deposit Ants.pheromone
     * on matrix "Ants.pheromone"
     * COMMENTS: this procedure could be implemented slightly faster, but it is
     * anyway not critical w.r.t. CPU time given that Ants.ras_ranks is
     * typically very small.
     */
    {
	int i, k, b, target;
	int[] help_b;

	// TRACE ( System.out.println("Rank-based Ant System Ants.pheromone deposit\n"); );

	help_b = new int[Ants.n_ants];
	for (k = 0; k < Ants.n_ants; k++)
	    help_b[k] = Ants.ant[k].tour_length;

	for (i = 0; i < Ants.ras_ranks - 1; i++) {
	    b = help_b[0];
	    target = 0;
	    for (k = 0; k < Ants.n_ants; k++) {
		if (help_b[k] < b) {
		    b = help_b[k];
		    target = k;
		}
	    }
	    help_b[target] = Integer.MAX_VALUE;
	    Ants.global_update_pheromone_weighted(Ants.ant[target], Ants.ras_ranks - i - 1);
	}
	Ants.global_update_pheromone_weighted(Ants.best_so_far_ant, Ants.ras_ranks);

    }

    static void mmas_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for MAX-MIN Ant System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: either the iteration-best or the best-so-far ant deposit Ants.pheromone
     * on matrix "Ants.pheromone"
     */
    {
	/*
	 * we use default upper Ants.pheromone trail limit for MMAS and hence we
	 * do not have to worry regarding keeping the upper limit
	 */

	int iteration_best_ant;

	// TRACE ( System.out.println("MAX-MIN Ant System Ants.pheromone deposit\n"); );

	if (InOut.iteration % Ants.u_gb == 0) {
	    iteration_best_ant = Ants.find_best();
	    Ants.global_update_pheromone(Ants.ant[iteration_best_ant]);
	} else {
	    if (Ants.u_gb == 1 && (InOut.iteration - InOut.restart_found_best > 50))
		Ants.global_update_pheromone(Ants.best_so_far_ant);
	    else
		Ants.global_update_pheromone(Ants.restart_best_ant);
	}

	if (LocalSearch.ls_flag != 0) {
	    /*
	     * implement the schedule for Ants.u_gb as defined in the
	     * Future Generation Computer Systems article or in Stuetzle's PhD thesis.
	     * This schedule is only applied if local search is used.
	     */
	    if ((InOut.iteration - InOut.restart_iteration) < 25)
		Ants.u_gb = 25;
	    else if ((InOut.iteration - InOut.restart_iteration) < 75)
		Ants.u_gb = 5;
	    else if ((InOut.iteration - InOut.restart_iteration) < 125)
		Ants.u_gb = 3;
	    else if ((InOut.iteration - InOut.restart_iteration) < 250)
		Ants.u_gb = 2;
	    else
		Ants.u_gb = 1;
	} else
	    Ants.u_gb = 25;

    }

    static void bwas_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for Best-Worst Ant System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: either the iteration-best or the best-so-far ant deposit Ants.pheromone
     * on matrix "Ants.pheromone"
     */
    {
	int iteration_worst_ant, distance_best_worst;

	// TRACE ( System.out.println("Best-worst Ant System Ants.pheromone deposit\n"); );

	Ants.global_update_pheromone(Ants.best_so_far_ant);
	iteration_worst_ant = Ants.find_worst();
	Ants.bwas_worst_ant_update(Ants.ant[iteration_worst_ant], Ants.best_so_far_ant);
	distance_best_worst = Ants.distance_between_ants(Ants.best_so_far_ant, Ants.ant[iteration_worst_ant]);
	/*
	 * System.out.println("distance_best_worst %ld, tour length worst %ld\n",distance_best_worst,ant[iteration_worst_ant
	 * ].tour_length);
	 */
	if (distance_best_worst < (int) (0.05 * (double) Tsp.n)) {
	    Ants.restart_best_ant.tour_length = Integer.MAX_VALUE;
	    Ants.init_pheromone_trails(Ants.trail_0);
	    InOut.restart_iteration = InOut.iteration;
	    InOut.restart_time = Timer.elapsed_time();
	    System.out.println("init Ants.pheromone trails with " + Ants.trail_0 + ", iteration " + InOut.iteration);
	} else
	    Ants.bwas_pheromone_mutation();
    }

    static void acs_global_update()
    /*
     * FUNCTION: manage global Ants.pheromone deposit for Ant Colony System
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: the best-so-far ant deposits Ants.pheromone on matrix "Ants.pheromone"
     * COMMENTS: global Ants.pheromone deposit in ACS is done per default using
     * the best-so-far ant; Gambardella & Dorigo examined also iteration-best
     * update (see their IEEE Trans. on Evolutionary Computation article),
     * but did not use it for the published computational results.
     */
    {
	// TRACE ( System.out.println("Ant colony System global Ants.pheromone deposit\n"); );

	Ants.global_acs_pheromone_update(Ants.best_so_far_ant);
    }

    static void pheromone_trail_update()
    /*
     * FUNCTION: manage global Ants.pheromone trail update for the ACO algorithms
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: Ants.pheromone trails are evaporated and Ants.pheromones are deposited
     * according to the rules defined by the various ACO algorithms.
     */
    {
	/*
	 * Simulate the Ants.pheromone evaporation of all Ants.pheromones; this is not necessary
	 * for ACS (see also ACO Book)
	 */
	if (Ants.as_flag || Ants.eas_flag || Ants.ras_flag || Ants.bwas_flag || Ants.mmas_flag) {
	    if (LocalSearch.ls_flag != 0) {
		if (Ants.mmas_flag)
		    Ants.mmas_evaporation_nn_list();
		else
		    Ants.evaporation_nn_list();
		/*
		 * evaporate only Ants.pheromones on arcs of candidate list to make the
		 * Ants.pheromone evaporation faster for being able to tackle large TSP
		 * Tsp.instances. For MMAS additionally check lower Ants.pheromone trail limits.
		 */
	    } else {
		/* if no local search is used, evaporate all Ants.pheromone trails */
		Ants.evaporation();
	    }
	}

	/* Next, apply the Ants.pheromone deposit for the various ACO algorithms */
	if (Ants.as_flag)
	    as_update();
	else if (Ants.eas_flag)
	    eas_update();
	else if (Ants.ras_flag)
	    ras_update();
	else if (Ants.mmas_flag)
	    mmas_update();
	else if (Ants.bwas_flag)
	    bwas_update();
	else if (Ants.acs_flag)
	    acs_global_update();

	/*
	 * check Ants.pheromone trail limits for MMAS; not necessary if local
	 * search is used, because in the local search case lower Ants.pheromone trail
	 * limits are checked in procedure mmas_evaporation_nn_list
	 */
	if (Ants.mmas_flag && LocalSearch.ls_flag == 0)
	    Ants.check_pheromone_trail_limits();

	/*
	 * Compute combined information Ants.pheromone times heuristic info after
	 * the Ants.pheromone update for all ACO algorithms except ACS; in the ACS case
	 * this is already done in the Ants.pheromone update procedures of ACS
	 */
	if (Ants.as_flag || Ants.eas_flag || Ants.ras_flag || Ants.mmas_flag || Ants.bwas_flag) {
	    if (LocalSearch.ls_flag != 0) {
		Ants.compute_nn_list_total_information();
	    } else {
		Ants.compute_total_information();
	    }
	}
    }

    /* --- main program ------------------------------------------------------ */

    public static void main(String[] args) {
	/*
	 * FUNCTION: main control for running the ACO algorithms
	 * INPUT: none
	 * OUTPUT: none
	 * (SIDE)EFFECTS: none
	 * COMMENTS: this function controls the run of "max_tries" independent trials
	 */
	for (String argument : args) {
	    System.out.println(argument);
	}
	Timer.start_timers();

	InOut.init_program(args);

	Tsp.instance.nn_list = Tsp.compute_nn_lists();
	Ants.pheromone = Utilities.generate_double_matrix(Tsp.n, Tsp.n);
	Ants.total = Utilities.generate_double_matrix(Tsp.n, Tsp.n);

	InOut.time_used = Timer.elapsed_time();
	System.out.println("Initialization took " + InOut.time_used + " seconds\n");

	for (InOut.n_try = 0; InOut.n_try < InOut.max_tries; InOut.n_try++) {

	    init_try(InOut.n_try);

	    while (!termination_condition()) {

		construct_solutions();

		if (LocalSearch.ls_flag > 0)
		    local_search();

		update_statistics();

		pheromone_trail_update();

		search_control_and_statistics();

		InOut.iteration++;
	    }
	    InOut.exit_try(InOut.n_try);
	}
	InOut.exit_program();

	// Added by AW
	int aw_best_tour_length = Utilities.best_of_vector(InOut.best_in_try, InOut.max_tries);
	String aw_best_tour = InOut.aw_best_tour_in_try[Utilities.aw_best_tour_index()];
	try {
	    Writer w = new OutputStreamWriter(new FileOutputStream("tour." + Tsp.instance.name), "UTF8");
	    BufferedWriter out = new BufferedWriter(w);
	    out.write(aw_best_tour_length + "\n");
	    out.write(aw_best_tour);
	    out.close();
	} catch (IOException e) {
	    System.err.print("Could not write file tour." + Tsp.instance.name + " " + e.getMessage());
	    System.exit(1);
	}
	System.out.println();
	System.out.println("Best tour:");
	System.out.println(aw_best_tour_length);
	System.out.println(aw_best_tour);
    }
}
