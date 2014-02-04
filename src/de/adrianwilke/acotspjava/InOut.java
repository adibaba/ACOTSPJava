package de.adrianwilke.acotspjava;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import de.adrianwilke.acotspjava.Tsp.problem;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class InOut {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: InOut.c
     * Author: Thomas Stuetzle
     * Purpose: mainly input / output / statistic routines
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

    enum Distance_type {
	EUC_2D, CEIL_2D, GEO, ATT
    };

    static Distance_type distance_type;

    public static final String PROG_ID_STR = "ACO algorithms for the TSP";

    static int[] best_in_try;
    static int[] best_found_at;
    static double[] time_best_found;
    static double[] time_total_run;
    static String[] aw_best_tour_in_try;

    static int n_try; /* try counter */
    static int n_tours; /* counter of number constructed tours */
    static int iteration; /* iteration counter */
    static int restart_iteration; /* remember iteration when restart was done if any */
    static double restart_time; /* remember time when restart was done if any */
    static int max_tries; /* maximum number of independent tries */
    static int max_tours; /* maximum number of tour constructions in one try */

    static double lambda; /* Parameter to determine branching factor */
    static double branch_fac; /* If branching factor < branch_fac => update trails */

    static double max_time; /* maximal allowed run time of a try */
    static double time_used; /* time used until some given event */
    static double time_passed; /* time passed until some moment */
    static int optimal; /* optimal solution or bound to find */

    static double mean_ants; /* average tour length */
    static double stddev_ants; /* stddev of tour lengths */
    static double branching_factor; /* average node branching factor when searching */
    static double found_branching; /* branching factor when best solution is found */

    static int found_best; /* iteration in which best solution is found */
    static int restart_found_best;/* iteration in which restart-best solution is found */

    /* ------------------------------------------------------------------------ */

    static File report, comp_report, stat_report;
    private static Map<String, BufferedWriter> writer = new HashMap<String, BufferedWriter>();

    static String name_buf;
    static int opt;
    static boolean quiet_flag; /* --quiet was given in the command-line. */

    static Tsp.point[] read_etsp(String tsp_file_name) throws IOException
    /*
     * FUNCTION: parse and read Tsp.instance file
     * INPUT: Tsp.instance name
     * OUTPUT: list of coordinates for all nodes
     * COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
     */
    {
	String buf;
	int i;
	Tsp.point[] nodeptr = null;

	if (tsp_file_name == null) {
	    System.err.println("No instance file specified, abort");
	    System.exit(1);
	}

	if (!new File(tsp_file_name).canRead()) {
	    System.err.println("Can not read file " + tsp_file_name);
	    System.exit(1);
	}

	System.out.println("\nreading tsp-file " + tsp_file_name + " ... ");

	i = 0;
	boolean found_coord_section = false;
	Reader reader = new InputStreamReader(new FileInputStream(tsp_file_name), "UTF8");
	BufferedReader bufferedReader = new BufferedReader(reader);
	String line = bufferedReader.readLine();
	while (line != null) {

	    if (line.startsWith("EOF")) {
		break;
	    }

	    if (!found_coord_section) {
		if (line.startsWith("NAME")) {
		    Tsp.instance.name = line.split(":")[1].trim();
		} else if (line.startsWith("COMMENT")) {
		} else if (line.startsWith("TYPE") && !line.contains("TSP")) {
		    System.err.println("Not a TSP Tsp.instance in TSPLIB format !!");
		    System.exit(1);
		} else if (line.startsWith("DIMENSION")) {
		    Tsp.n = Integer.parseInt(line.split(":")[1].trim());
		    Tsp.instance.n = Tsp.n;
		    nodeptr = new Tsp.point[Tsp.n];
		    assert (Tsp.n > 2 && Tsp.n < 6000);
		} else if (line.startsWith("DISPLAY_DATA_TYPE")) {
		} else if (line.startsWith("EDGE_WEIGHT_TYPE")) {
		    buf = line.split(":")[1].trim();
		    if (buf.equals("EUC_2D")) {
			distance_type = Distance_type.EUC_2D;
		    } else if (buf.equals("CEIL_2D")) {
			distance_type = Distance_type.CEIL_2D;
		    } else if (buf.equals("GEO")) {
			distance_type = Distance_type.GEO;
		    } else if (buf.equals("ATT")) {
			distance_type = Distance_type.ATT;
		    } else {
			System.err.println("EDGE_WEIGHT_TYPE " + buf + " not implemented");
			System.exit(1);
		    }
		}
	    } else {
		String[] city_info = line.split(" ");
		nodeptr[i] = new Tsp.point();
		nodeptr[i].x = Double.parseDouble(city_info[1]);
		nodeptr[i].y = Double.parseDouble(city_info[2]);
		i++;
	    }

	    if (line.startsWith("NODE_COORD_SECTION")) {
		found_coord_section = true;
	    }

	    line = bufferedReader.readLine();
	}

	if (!found_coord_section) {
	    System.err.println("Some error ocurred finding start of coordinates from tsp file !!");
	    System.exit(1);
	}

	bufferedReader.close();

	// TRACE ( System.out.println("number of cities is %ld\n",Tsp.n); )
	// TRACE ( System.out.println("\n... done\n"); )
	System.out.println();

	return (nodeptr);
    }

    static void write_report()
    /*
     * FUNCTION: output some info about trial (best-so-far solution quality, time)
     * INPUT: none
     * OUTPUT: none
     * COMMENTS: none
     */
    {
	System.out.println("best " + Ants.best_so_far_ant.tour_length + ", iteration: " + iteration + ", time "
		+ Timer.elapsed_time());
	if (comp_report != null)
	    printToFile(comp_report, "best " + Ants.best_so_far_ant.tour_length + "\t iteration " + iteration
		    + "\t tours " + n_tours + "\t time " + time_used);
    }

    static void fprintf_parameters(File file) {
	printToFile(file, "max_tries\t\t " + max_tries);
	printToFile(file, "max_tours\t\t " + max_tours);
	printToFile(file, "max_time\t\t " + max_time);
	printToFile(file, "Utilities.seed\t\t " + Utilities.seed);
	printToFile(file, "optimum\t\t\t " + optimal);
	printToFile(file, "n_ants\t\t\t " + Ants.n_ants);
	printToFile(file, "Ants.nn_ants\t\t " + Ants.nn_ants);
	printToFile(file, "Ants.alpha\t\t " + Ants.alpha);
	printToFile(file, "Ants.beta\t\t " + Ants.beta);
	printToFile(file, "Ants.rho\t\t " + Ants.rho);
	printToFile(file, "Ants.q_0\t\t " + Ants.q_0);
	printToFile(file, "Ants.elitist_ants\t " + Ants.elitist_ants);
	printToFile(file, "Ants.ras_ranks\t\t " + Ants.ras_ranks);
	printToFile(file, "LocalSearch.ls_flag\t " + LocalSearch.ls_flag);
	printToFile(file, "LocalSearch.nn_ls\t " + LocalSearch.nn_ls);
	printToFile(file, "LocalSearch.dlb_flag\t " + LocalSearch.dlb_flag);
	printToFile(file, "Ants.as_flag\t\t " + Ants.as_flag);
	printToFile(file, "Ants.eAnts.as_flag\t " + Ants.eas_flag);
	printToFile(file, "rAnts.as_flag\t\t " + Ants.ras_flag);
	printToFile(file, "mmAnts.as_flag\t\t " + Ants.mmas_flag);
	printToFile(file, "Ants.bwAnts.as_flag\t " + Ants.bwas_flag);
	printToFile(file, "Ants.acs_flag\t\t " + Ants.acs_flag);
    }

    static void print_default_parameters()
    /*
     * FUNCTION: output default parameter settings
     * INPUT: none
     * OUTPUT: none
     * COMMENTS: none
     */
    {
	System.err.println("\nDefault parameter settings are:\n\n");
	fprintf_parameters(null);
    }

    static void set_default_as_parameters() {
	assert (Ants.as_flag);
	Ants.n_ants = -1; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.5;
	Ants.q_0 = 0.0;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = 0;
    }

    static void set_default_eas_parameters() {
	assert (Ants.eas_flag);
	Ants.n_ants = -1; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.5;
	Ants.q_0 = 0.0;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = Ants.n_ants;
    }

    static void set_default_ras_parameters() {
	assert (Ants.ras_flag);
	Ants.n_ants = -1; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.1;
	Ants.q_0 = 0.0;
	Ants.ras_ranks = 6;
	Ants.elitist_ants = 0;
    }

    static void set_default_bwas_parameters() {
	assert (Ants.bwas_flag);
	Ants.n_ants = -1; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.1;
	Ants.q_0 = 0.0;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = 0;
    }

    static void set_default_mmas_parameters() {
	assert (Ants.mmas_flag);
	Ants.n_ants = -1; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.02;
	Ants.q_0 = 0.0;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = 0;
    }

    static void set_default_acs_parameters() {
	assert (Ants.acs_flag);
	Ants.n_ants = 10; /* number of ants (-1 means Tsp.instance size) */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.1;
	Ants.q_0 = 0.9;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = 0;
    }

    static void set_default_ls_parameters() {
	assert (LocalSearch.ls_flag != 0);
	LocalSearch.dlb_flag = true; /* apply don't look bits in local search */
	LocalSearch.nn_ls = 20; /* use fixed radius search in the 20 nearest neighbours */

	Ants.n_ants = 25; /* number of ants */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.5;
	Ants.q_0 = 0.0;

	if (Ants.mmas_flag) {
	    Ants.n_ants = 25;
	    Ants.rho = 0.2;
	    Ants.q_0 = 0.00;
	} else if (Ants.acs_flag) {
	    Ants.n_ants = 10;
	    Ants.rho = 0.1;
	    Ants.q_0 = 0.98;
	} else if (Ants.eas_flag) {
	    Ants.elitist_ants = Ants.n_ants;
	}
    }

    static void set_default_parameters()
    /*
     * FUNCTION: set default parameter settings
     * INPUT: none
     * OUTPUT: none
     * COMMENTS: none
     */
    {
	LocalSearch.ls_flag = 3; /* per default run 3-opt */
	LocalSearch.dlb_flag = true; /* apply don't look bits in local search */
	LocalSearch.nn_ls = 20; /* use fixed radius search in the 20 nearest neighbours */
	Ants.n_ants = 25; /* number of ants */
	Ants.nn_ants = 20; /* number of nearest neighbours in tour construction */
	Ants.alpha = 1.0;
	Ants.beta = 2.0;
	Ants.rho = 0.5;
	Ants.q_0 = 0.0;
	max_tries = 10;
	max_tours = 0;
	Utilities.seed = (int) System.currentTimeMillis();
	max_time = 10.0;
	optimal = 1;
	branch_fac = 1.00001;
	Ants.u_gb = Integer.MAX_VALUE;
	Ants.as_flag = false;
	Ants.eas_flag = false;
	Ants.ras_flag = false;
	Ants.mmas_flag = true;
	Ants.bwas_flag = false;
	Ants.acs_flag = false;
	Ants.ras_ranks = 0;
	Ants.elitist_ants = 0;
    }

    static void population_statistics()
    /*
     * FUNCTION: compute some population statistics like average tour length,
     * standard deviations, average distance, branching-factor and
     * output to a file gathering statistics
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: none
     */
    {
	int j, k;
	int[] l;
	double pop_mean, pop_stddev, avg_distance = 0.0;

	l = new int[Ants.n_ants];
	for (k = 0; k < Ants.n_ants; k++) {
	    l[k] = Ants.ant[k].tour_length;
	}

	pop_mean = Utilities.mean(l, Ants.n_ants);
	pop_stddev = Utilities.std_deviation(l, Ants.n_ants, pop_mean);
	branching_factor = node_branching(lambda);

	for (k = 0; k < Ants.n_ants - 1; k++)
	    for (j = k + 1; j < Ants.n_ants; j++) {
		avg_distance += (double) Ants.distance_between_ants(Ants.ant[k], Ants.ant[j]);
	    }
	avg_distance /= ((double) Ants.n_ants * (double) (Ants.n_ants - 1) / 2.);

	if (stat_report != null) {
	    printToFile(stat_report, iteration + "\t" + pop_mean + "\t" + pop_stddev + "\t" + (pop_stddev / pop_mean)
		    + "\t" + branching_factor + "\t" + ((branching_factor - 1.) * (double) Tsp.n) + "\t" + avg_distance
		    + "\t" + (avg_distance / (double) Tsp.n));
	}
    }

    static double node_branching(double l)
    /*
     * FUNCTION: compute the average node lambda-branching factor
     * INPUT: lambda value
     * OUTPUT: average node branching factor
     * (SIDE)EFFECTS: none
     * COMMENTS: see the ACO book for a definition of the average node
     * lambda-branching factor
     */
    {
	int i, m;
	double min, max, cutoff;
	double avg;
	double[] num_branches;

	num_branches = new double[Tsp.n];

	for (m = 0; m < Tsp.n; m++) {
	    /* determine max, min to calculate the cutoff value */
	    min = Ants.pheromone[m][Tsp.instance.nn_list[m][1]];
	    max = Ants.pheromone[m][Tsp.instance.nn_list[m][1]];
	    for (i = 1; i < Ants.nn_ants; i++) {
		if (Ants.pheromone[m][Tsp.instance.nn_list[m][i]] > max)
		    max = Ants.pheromone[m][Tsp.instance.nn_list[m][i]];
		if (Ants.pheromone[m][Tsp.instance.nn_list[m][i]] < min)
		    min = Ants.pheromone[m][Tsp.instance.nn_list[m][i]];
	    }
	    cutoff = min + l * (max - min);

	    for (i = 0; i < Ants.nn_ants; i++) {
		if (Ants.pheromone[m][Tsp.instance.nn_list[m][i]] > cutoff)
		    num_branches[m] += 1.;
	    }
	}
	avg = 0.;
	for (m = 0; m < Tsp.n; m++) {
	    avg += num_branches[m];
	}
	/* Norm branching factor to minimal value 1 */
	return (avg / (double) (Tsp.n * 2));
    }

    static void output_solution()
    /*
     * FUNCTION: output a solution together with node coordinates
     * INPUT: none
     * OUTPUT: none
     * COMMENTS: not used in the default implementation but may be useful anyway
     */
    {

	int i;
	if (stat_report != null) {
	    for (i = 0; i < Tsp.n; i++) {
		printToFile(stat_report, Ants.best_so_far_ant.tour[i] + " "
			+ Tsp.instance.nodeptr[Ants.best_so_far_ant.tour[i]].x + " "
			+ Tsp.instance.nodeptr[Ants.best_so_far_ant.tour[i]].y);
	    }

	}
    }

    static void exit_try(int ntry)
    /*
     * FUNCTION: save some statistical information on a trial once it finishes
     * INPUT: trial number
     * OUTPUT: none
     * COMMENTS:
     */
    {
	checkTour(Ants.best_so_far_ant.tour);
	/* printTourFile( best_so_far_ant.tour ); */

	System.out.println("Best Solution in try " + ntry + " is " + Ants.best_so_far_ant.tour_length);

	if (report != null)
	    printToFile(report, "Best: " + Ants.best_so_far_ant.tour_length + "\t Iterations: " + found_best
		    + "\t B-Fac " + found_branching + "\t Time " + time_used + "\t Tot.time " + Timer.elapsed_time());
	System.out.println(" Best Solution was found after " + found_best + " iterations\n");

	best_in_try[ntry] = Ants.best_so_far_ant.tour_length;
	best_found_at[ntry] = found_best;
	time_best_found[ntry] = time_used;
	time_total_run[ntry] = Timer.elapsed_time();
	aw_best_tour_in_try[ntry] = Arrays.toString(Ants.best_so_far_ant.tour);

	System.out.println("\ntry " + ntry + ", Best " + best_in_try[ntry] + ", found at iteration "
		+ best_found_at[ntry] + ", found at time " + time_best_found[ntry] + "\n");

	if (comp_report != null)
	    printToFile(comp_report, "end try " + ntry + "\n");
	if (stat_report != null)
	    printToFile(stat_report, "end try " + ntry + "\n");
	// TRACE (output_solution();)

    }

    static void exit_program()
    /*
     * FUNCTION: save some final statistical information on a trial once it finishes
     * INPUT: none
     * OUTPUT: none
     * COMMENTS:
     */
    {
	int best_tour_length, worst_tour_length;
	double t_avgbest, t_stdbest, t_avgtotal, t_stdtotal;
	double avg_sol_quality = 0., avg_cyc_to_bst = 0., stddev_best, stddev_iterations;

	best_tour_length = Utilities.best_of_vector(best_in_try, max_tries);
	worst_tour_length = Utilities.worst_of_vector(best_in_try, max_tries);

	avg_cyc_to_bst = Utilities.mean(best_found_at, max_tries);
	stddev_iterations = Utilities.std_deviation(best_found_at, max_tries, avg_cyc_to_bst);

	avg_sol_quality = Utilities.mean(best_in_try, max_tries);
	stddev_best = Utilities.std_deviation(best_in_try, max_tries, avg_sol_quality);

	t_avgbest = Utilities.meanr(time_best_found, max_tries);
	System.out.println(" t_avgbest = " + t_avgbest);
	t_stdbest = Utilities.std_deviationr(time_best_found, max_tries, t_avgbest);

	t_avgtotal = Utilities.meanr(time_total_run, max_tries);
	System.out.println(" t_avgtotal = " + t_avgtotal);
	t_stdtotal = Utilities.std_deviationr(time_total_run, max_tries, t_avgtotal);

	if (report != null) {
	    printToFile(report, "\nAverage-Best: " + avg_sol_quality + "\t Average-Iterations: " + avg_cyc_to_bst);
	    printToFile(report, "Stddev-Best: " + stddev_best + " \t Stddev Iterations: " + stddev_iterations);
	    printToFile(report, "Best try: " + best_tour_length + "\t\t Worst try: " + worst_tour_length);
	    printToFile(report, "\nAvg.time-best: " + t_avgbest + " stddev.time-best: " + t_stdbest);
	    printToFile(report, "\nAvg.time-Ants.total: " + t_avgtotal + " stddev.time-Ants.total: " + t_stdtotal);

	    if (optimal > 0) {
		printToFile(report, " excess best = " + ((double) (best_tour_length - optimal) / (double) optimal)
			+ ", excess average = " + ((double) (avg_sol_quality - optimal) / (double) optimal) + ","
			+ " excess worst = " + ((double) (worst_tour_length - optimal) / (double) optimal));
	    }
	}

	if (comp_report != null)
	    printToFile(comp_report, "end problem " + Tsp.instance.name);

	for (String key : writer.keySet()) {
	    try {
		writer.get(key).close();
	    } catch (IOException e) {
		System.err.println("Could not close file " + key + " " + e.getMessage());
	    }
	}
    }

    static void init_program(String[] args)
    /*
     * FUNCTION: initialize the program,
     * INPUT: program arguments, needed for parsing commandline
     * OUTPUT: none
     * COMMENTS:
     */
    {
	Tsp.instance = new problem();

	String temp_buffer;

	System.out.println(InOut.PROG_ID_STR);
	set_default_parameters();
	Parse.parse_commandline(args);

	assert (max_tries <= Utilities.MAXIMUM_NO_TRIES);

	best_in_try = new int[max_tries];
	best_found_at = new int[max_tries];
	time_best_found = new double[max_tries];
	time_total_run = new double[max_tries];

	aw_best_tour_in_try = new String[max_tries];

	// TRACE ( System.out.println("read problem data  ..\n\n"); )

	try {
	    Tsp.instance.nodeptr = read_etsp(name_buf);
	} catch (IOException e) {
	    System.err.println("Could not read input file. " + e.getMessage());
	    System.exit(1);
	}

	// TRACE ( System.out.println("\n .. done\n\n"); )

	if (Ants.n_ants < 0)
	    Ants.n_ants = Tsp.n;
	/*
	 * default setting for Ants.elitist_ants is 0; if EAS is applied and
	 * option Ants.elitist_ants is not used, we set the default to
	 * Ants.elitist_ants = n
	 */
	if (Ants.eas_flag && Ants.elitist_ants <= 0)
	    Ants.elitist_ants = Tsp.n;

	LocalSearch.nn_ls = Math.min(Tsp.n - 1, LocalSearch.nn_ls);

	assert (Ants.n_ants < Ants.MAX_ANTS - 1);
	assert (Ants.nn_ants < Ants.MAX_NEIGHBOURS);
	assert (Ants.nn_ants > 0);
	assert (LocalSearch.nn_ls > 0);

	if (!quiet_flag) {
	    Writer w;
	    try {
		temp_buffer = "best." + Tsp.instance.name;
		// // TRACE ( System.out.println("%s\n",temp_buffer); )
		report = new File(temp_buffer);
		w = new OutputStreamWriter(new FileOutputStream(temp_buffer), "UTF8");
		writer.put(report.getName(), new BufferedWriter(w));

		temp_buffer = "cmp." + Tsp.instance.name;
		// // TRACE ( System.out.println("%s\n",temp_buffer); )
		comp_report = new File(temp_buffer);
		w = new OutputStreamWriter(new FileOutputStream(temp_buffer), "UTF8");
		writer.put(comp_report.getName(), new BufferedWriter(w));

		temp_buffer = "stat." + Tsp.instance.name;
		// // TRACE ( System.out.println("%s\n",temp_buffer); )
		stat_report = new File(temp_buffer);
		w = new OutputStreamWriter(new FileOutputStream(temp_buffer), "UTF8");
		writer.put(stat_report.getName(), new BufferedWriter(w));
	    } catch (IOException e) {
		System.err.println("Could not write file. " + e.getMessage());
		System.exit(1);
	    }
	} else {
	    report = null;
	    comp_report = null;
	    stat_report = null;
	}

	System.out.println("calculating distance matrix ..");
	Tsp.instance.distance = Tsp.compute_distances();
	System.out.println(" .. done\n");
	write_params();

	if (comp_report != null)
	    printToFile(comp_report, "begin problem " + name_buf);
	System.out.println("allocate ants' memory ..");
	Ants.allocate_ants();
	System.out.println(" .. done\n");
    }

    static void printDist()
    /*
     * FUNCTION: print distance matrix
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j;

	System.out.println("Distance Matrix:\n");
	for (i = 0; i < Tsp.n; i++) {
	    System.out.println("From " + i);
	    for (j = 0; j < Tsp.n - 1; j++) {
		System.out.println(" " + Tsp.instance.distance[i][j]);
	    }
	    System.out.println(" " + Tsp.instance.distance[i][Tsp.n - 1]);
	    System.out.println("\n");
	}
	System.out.println("\n");
    }

    static void printHeur()
    /*
     * FUNCTION: print heuristic information
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j;

	System.out.println("Heuristic information:\n");
	for (i = 0; i < Tsp.n; i++) {
	    System.out.println("From " + i + ":  ");
	    for (j = 0; j < Tsp.n - 1; j++) {
		System.out.println(" " + Ants.HEURISTIC(i, j));
	    }
	    System.out.println(" " + Ants.HEURISTIC(i, j));
	    System.out.println("\n");
	}
	System.out.println("\n");
    }

    static void printTrail()
    /*
     * FUNCTION: print Ants.pheromone trail values
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j;

	System.out.println("Ants.pheromone Trail matrix, iteration: " + iteration + "\n");
	for (i = 0; i < Tsp.n; i++) {
	    System.out.println("From " + i + ": ");
	    for (j = 0; j < Tsp.n; j++) {
		System.out.println(" " + Ants.pheromone[i][j] + " ");
		if (Ants.pheromone[i][j] > 1.0)
		    System.out.println("XXXXX\n");
	    }
	    System.out.println("\n");
	}
	System.out.println("\n");
    }

    static void printTotal()
    /*
     * FUNCTION: print values of Ants.pheromone times heuristic information
     * INPUT: none
     * OUTPUT: none
     */
    {
	int i, j;

	System.out.println("combined Ants.pheromone and heuristic info\n\n");
	for (i = 0; i < Tsp.n; i++) {
	    for (j = 0; j < Tsp.n - 1; j++) {
		System.out.println(" " + Ants.total[i][j]);
		if (Ants.total[i][j] > 1.0)
		    System.out.println("XXXXX\n");
	    }
	    System.out.println(" " + Ants.total[i][Tsp.n - 1]);
	    if (Ants.total[i][Tsp.n - 1] > 1.0)
		System.out.println("XXXXX\n");
	}
	System.out.println("\n");
    }

    static void printProbabilities()
    /*
     * FUNCTION: prints the selection probabilities as encountered by an ant
     * INPUT: none
     * OUTPUT: none
     * COMMENTS: this computation assumes that no choice has been made yet.
     */
    {
	int i, j;
	double p[];
	double sum_prob;

	System.out.println("Selection Probabilities, iteration: " + iteration);
	p = new double[Tsp.n];

	for (i = 0; i < Tsp.n; i++) {
	    System.out.println("From " + i);
	    sum_prob = 0.;
	    for (j = 0; j < Tsp.n; j++) {
		if (i == j)
		    p[j] = 0.;
		else
		    p[j] = Ants.total[i][j];
		sum_prob += p[j];
	    }
	    for (j = 0; j < Tsp.n; j++) {
		p[j] = p[j] / sum_prob;
	    }
	    for (j = 0; j < Tsp.n - 1; j++) {
		System.out.println(" " + p[j] + " ");
	    }
	    System.out.println(" " + p[Tsp.n - 1]);
	    if ((j % 26) == 0) {
		System.out.println("\n");
	    }
	    System.out.println("\n");
	}
	System.out.println("\n");
    }

    static void printTour(int[] t)
    /*
     * FUNCTION: print the tour *t
     * INPUT: pointer to a tour
     * OUTPUT: none
     */
    {
	int i;

	System.out.println("\n");
	for (i = 0; i <= Tsp.n; i++) {
	    if (i % 25 == 0)
		System.out.println("\n");
	    System.out.println(t[i]);
	}
	System.out.println("\n");
	System.out.println("Tour length = " + Tsp.compute_tour_length(t));
    }

    static void printTourFile(int[] t)
    /*
     * FUNCTION: print the tour *t to cmp.tsplibfile
     * INPUT: pointer to a tour
     * OUTPUT: none
     */
    {
	int i;
	if (comp_report == null)
	    return;

	printToFile(comp_report, "begin solution\n");
	for (i = 0; i < Tsp.n; i++) {
	    printToFile(comp_report, t[i] + " ");
	}
	printToFile(comp_report, "\n");
	printToFile(comp_report, "Tour length " + Tsp.compute_tour_length(t));
	printToFile(comp_report, "end solution\n");
    }

    static void checkTour(int[] t)
    /*
     * FUNCTION: make a simple check whether tour *t can be feasible
     * INPUT: pointer to a tour
     * OUTPUT: none
     */
    {
	int i, sum = 0;

	for (i = 0; i < Tsp.n; i++) {
	    sum += t[i];
	}
	if (sum != (Tsp.n - 1) * Tsp.n / 2) {
	    System.err.println("Next tour must be flawed !!\n");
	    printTour(t);
	    System.exit(1);
	}
    }

    static void write_params()
    /*
     * FUNCTION: writes chosen parameter settings in standard output and in
     * report files
     * INPUT: none
     * OUTPUT: none
     */
    {
	System.out.println("Parameter-settings:\n");
	fprintf_parameters(null);
	System.out.println("\n");

	if (report != null) {
	    printToFile(report, "Parameter-settings: \n\n");
	    fprintf_parameters(report);
	    printToFile(report, "\n");
	}

	if (comp_report != null) {
	    printToFile(comp_report, PROG_ID_STR);
	    printToFile(comp_report, "Parameter-settings: \n\n");
	    fprintf_parameters(comp_report);
	    printToFile(comp_report, "\n");
	}
    }

    static void printToFile(File file, String string) {
	if (file == null) {
	    System.out.println(string);
	} else {
	    try {
		writer.get(file.getName()).write(string + "\n");
	    } catch (IOException e) {
		System.err.print("Could not write file " + file.getName() + " " + e.getMessage());
		System.exit(1);
	    }
	}
    }
}
