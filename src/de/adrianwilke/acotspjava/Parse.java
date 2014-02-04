package de.adrianwilke.acotspjava;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class Parse {
    /***************************************************************************
     * Program's name: ACOTSPJava
     * 
     * Command line parser for 'ACO algorithms for the TSP'
     * 
     * Copyright (C) 2014 Adrian Wilke
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
     * You should have received a copy of the GNU General Public License along
     * with this program; if not, write to the Free Software Foundation, Inc.,
     * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
     ***************************************************************************/

    static class OptComparator implements Comparator<Option> {

	Map<String, Integer> opt = new HashMap<String, Integer>();

	public OptComparator() {
	    int i = 0;
	    opt.put("r", i++);
	    opt.put("s", i++);
	    opt.put("t", i++);
	    opt.put("seed", i++);
	    opt.put("i", i++);
	    opt.put("o", i++);
	    opt.put("g", i++);
	    opt.put("a", i++);
	    opt.put("b", i++);
	    opt.put("e", i++);
	    opt.put("q", i++);
	    opt.put("c", i++);
	    opt.put("f", i++);
	    opt.put("k", i++);
	    opt.put("l", i++);
	    opt.put("d", i++);
	    opt.put("u", i++);
	    opt.put("v", i++);
	    opt.put("w", i++);
	    opt.put("x", i++);
	    opt.put("quiet", i++);
	    opt.put("h", i++);
	}

	@Override
	public int compare(Option o1, Option o2) {
	    if (o1.getValue() == null || o2.getValue() == null)
		return 0;
	    else
		return (opt.get(o1.getOpt()) - opt.get(o2.getOpt()));
	}
    }

    static int parse_commandline(String args[]) {

	// TODO range check

	if (args.length == 0) {
	    System.err.println("No options are specified.");
	    System.err.println("Try `--help' for more information.");
	    System.exit(1);
	}

	Options options = new Options();
	options.addOption("r", "tries", true, "# number of independent trials");
	options.addOption("s", "tours", true, "# number of steps in each trial");
	options.addOption("t", "time", true, "# maximum time for each trial");
	options.addOption("seed", true, "# seed for the random number generator");
	options.addOption("i", "tsplibfile", true, "f inputfile (TSPLIB format necessary)");
	options.addOption("o", "optimum", true, "# stop if tour better or equal optimum is found");
	options.addOption("m", "ants", true, "# number of ants");
	options.addOption("g", "nnants", true, "# nearest neighbours in tour construction");
	options.addOption("a", "alpha", true, "# alpha (influence of pheromone trails)");
	options.addOption("b", "beta", true, "# beta (influence of heuristic information)");
	options.addOption("e", "rho", true, "# rho: pheromone trail evaporation");
	options.addOption("q", "q0", true, "# q_0: prob. of best choice in tour construction");
	options.addOption("c", "elitistants", true, "# number of elitist ants");
	options.addOption("f", "rasranks", true, "# number of ranks in rank-based Ant System");
	options.addOption("k", "nnls", true, "# No. of nearest neighbors for local search");
	options.addOption("l", "localsearch", true, "0:no local search  1:2-opt  2:2.5-opt  3:3-opt");
	options.addOption("d", "dlb", false, "1 use don't look bits in local search");
	options.addOption("u", "as", false, "apply basic Ant System");
	options.addOption("v", "eas", false, "apply elitist Ant System");
	options.addOption("w", "ras", false, "apply rank-based version of Ant System");
	options.addOption("x", "mmas", false, "apply MAX-MIN ant_colony system");
	options.addOption("y", "bwas", false, "apply best-worst ant_colony system");
	options.addOption("z", "acs", false, "apply ant_colony colony system");
	options.addOption("quiet", false, "reduce output to a minimum, no extra files");
	options.addOption("h", "help", false, "display this help text and exit");

	CommandLine cmd = null;
	CommandLineParser parser = new BasicParser();
	try {
	    cmd = parser.parse(options, args);
	} catch (ParseException e) {
	    System.err.println("Error: " + e.getMessage());
	    System.exit(1);
	}

	if (cmd.hasOption("h")) {
	    HelpFormatter formatter = new HelpFormatter();
	    formatter.setSyntaxPrefix("Usage: ");
	    formatter.setOptionComparator(new OptComparator());
	    formatter.printHelp(InOut.PROG_ID_STR + " [OPTION]... [ARGUMENT]...", "Options:", options, "");
	    System.exit(0);
	}

	System.out.println("OPTIONS:");

	if (cmd.hasOption("quiet")) {
	    InOut.quiet_flag = true;
	    System.out.println("-quiet Quiet mode is set");
	}

	if (cmd.hasOption("t")) {
	    InOut.max_time = Float.parseFloat(cmd.getOptionValue("t"));
	    System.out.println("-t/time Time limit with argument " + InOut.max_time);
	} else {
	    System.out.println("Note: Time limit is set to default " + InOut.max_time + " seconds");
	}

	if (cmd.hasOption("r")) {
	    InOut.max_tries = Integer.parseInt(cmd.getOptionValue("r"));
	    System.out.println("-r/tries Number of tries with argument " + InOut.max_tries);
	} else {
	    System.out.println("Note: Number of tries is set to default " + InOut.max_tries);
	}

	if (cmd.hasOption("s")) {
	    InOut.max_tours = Integer.parseInt(cmd.getOptionValue("s"));
	    System.out.println("-s/tours Maximum number tours with argument " + InOut.max_tours);
	} else {
	    System.out.println("Note: Maximum number tours is set to default " + InOut.max_tours);
	}

	if (cmd.hasOption("seed")) {
	    Utilities.seed = Integer.parseInt(cmd.getOptionValue("seed"));
	    System.out.println("-seed with argument " + Utilities.seed);
	} else {
	    System.out.println("Note: A seed was generated as " + Utilities.seed);
	}

	if (cmd.hasOption("o")) {
	    InOut.optimal = Integer.parseInt(cmd.getOptionValue("o"));
	    System.out.println("-o/optimum Optimal solution with argument " + InOut.optimal);
	} else {
	    System.out.println("Note: Optimal solution value is set to default " + InOut.optimal);
	}

	if (cmd.hasOption("i")) {
	    InOut.name_buf = cmd.getOptionValue("i");
	    System.out.println("-i/tsplibfile File with argument " + InOut.name_buf);
	} else {
	    System.err.println("Error: No input file given");
	    System.exit(1);
	}

	// Choice of ONE algorithm
	int algorithmCount = 0;
	if (cmd.hasOption("u")) {
	    algorithmCount++;
	}
	if (cmd.hasOption("w")) {
	    algorithmCount++;
	}
	if (cmd.hasOption("x")) {
	    algorithmCount++;
	}
	if (cmd.hasOption("v")) {
	    algorithmCount++;
	}
	if (cmd.hasOption("y")) {
	    algorithmCount++;
	}
	if (cmd.hasOption("z")) {
	    algorithmCount++;
	}
	if (algorithmCount > 1) {
	    System.err.println("Error: More than one ACO algorithm enabled in the command line.");
	    System.exit(1);
	} else if (algorithmCount == 1) {
	    Ants.as_flag = false;
	    Ants.eas_flag = false;
	    Ants.ras_flag = false;
	    Ants.mmas_flag = false;
	    Ants.bwas_flag = false;
	    Ants.acs_flag = false;
	}

	if (cmd.hasOption("u")) {
	    Ants.as_flag = true;
	    InOut.set_default_as_parameters();
	    System.out.println("-u/as is set, run basic Ant System");
	}
	if (cmd.hasOption("v")) {
	    Ants.eas_flag = true;
	    InOut.set_default_eas_parameters();
	    System.out.println("-v/eas is set, run Elitist Ant System");
	}
	if (cmd.hasOption("w")) {
	    Ants.ras_flag = true;
	    InOut.set_default_ras_parameters();
	    System.out.println("-w/ras is set, run Rank-based Ant System");
	}
	if (cmd.hasOption("x") || algorithmCount == 0) {
	    Ants.mmas_flag = true;
	    InOut.set_default_mmas_parameters();
	    System.out.println("-x/mmas is set, run MAX-MIN Ant System");
	}
	if (cmd.hasOption("y")) {
	    Ants.bwas_flag = true;
	    InOut.set_default_bwas_parameters();
	    System.out.println("-y/bwas is set, run Best-Worst Ant System");
	}
	if (cmd.hasOption("z")) {
	    Ants.acs_flag = true;
	    InOut.set_default_acs_parameters();
	    System.out.println("-z/acs is set, run Ant Colony System");
	}

	// Local search
	if (cmd.hasOption("l")) {
	    LocalSearch.ls_flag = Integer.parseInt(cmd.getOptionValue("l"));

	    switch (LocalSearch.ls_flag) {
	    case 0:
		System.out.println("Note: local search flag is set to default 0 (disabled)");
		break;
	    case 1:
		System.out.println("Note: local search flag is set to default 1 (2-opt)");
		break;
	    case 2:
		System.out.println("Note: local search flag is set to default 2 (2.5-opt)");
		break;
	    case 3:
		System.out.println("Note: local search flag is set to default 3 (3-opt)");
		break;
	    default:
		System.out.println("-l/localsearch with argument " + LocalSearch.ls_flag);
		break;
	    }
	}
	if (LocalSearch.ls_flag != 0) {
	    InOut.set_default_ls_parameters();
	}

	if (cmd.hasOption("m")) {
	    Ants.n_ants = Integer.parseInt(cmd.getOptionValue("m"));
	    System.out.println("-m/ants Number of ants with argument " + Ants.n_ants);
	} else {
	    System.out.println("Note: Number of ants is set to default " + Ants.n_ants);
	}

	if (cmd.hasOption("a")) {
	    Ants.alpha = Float.parseFloat(cmd.getOptionValue("a"));
	    System.out.println("-a/alpha with argument " + Ants.alpha);
	} else {
	    System.out.println("Note: Alpha is set to default " + Ants.alpha);
	}

	if (cmd.hasOption("b")) {
	    Ants.beta = Float.parseFloat(cmd.getOptionValue("b"));
	    System.out.println("-b/beta with argument " + Ants.beta);
	} else {
	    System.out.println("Note: Beta is set to default " + Ants.beta);
	}

	if (cmd.hasOption("e")) {
	    Ants.rho = Float.parseFloat(cmd.getOptionValue("e"));
	    System.out.println("-e/rho with argument " + Ants.rho);
	} else {
	    System.out.println("Note: Rho is set to default " + Ants.rho);
	}

	if (cmd.hasOption("q")) {
	    Ants.q_0 = Float.parseFloat(cmd.getOptionValue("q"));
	    System.out.println("-q/q0 with argument " + Ants.q_0);
	} else {
	    System.out.println("Note: q0 is set to default " + Ants.q_0);
	}

	if (cmd.hasOption("c")) {
	    Ants.elitist_ants = Integer.parseInt(cmd.getOptionValue("c"));
	    System.out.println("-c/elitistants Number of elitist ants with argument " + Ants.elitist_ants);
	} else {
	    System.out.println("Note: Number of elitist ants is set to default " + Ants.elitist_ants);
	}

	if (cmd.hasOption("f")) {
	    Ants.ras_ranks = Integer.parseInt(cmd.getOptionValue("f"));
	    System.out.println("-f/rasranks Number of ranks with argument " + Ants.ras_ranks);
	} else {
	    System.out.println("Note: Number of ranks is set to default " + Ants.ras_ranks);
	}

	if (cmd.hasOption("k")) {
	    LocalSearch.nn_ls = Integer.parseInt(cmd.getOptionValue("k"));
	    System.out.println("-k/nnls Number nearest neighbours with argument " + LocalSearch.nn_ls);
	} else {
	    System.out
		    .println("Note: Number nearest neighbours in local search is set to default " + LocalSearch.nn_ls);
	}

	if (cmd.hasOption("d")) {
	    LocalSearch.dlb_flag = true;
	    System.out.println("-d/dlb Don't-look-bits flag with argument " + LocalSearch.dlb_flag);
	} else {
	    System.out.println("Note: Don't-look-bits flag is set to default " + LocalSearch.dlb_flag);
	}

	return 0;
    }
}
