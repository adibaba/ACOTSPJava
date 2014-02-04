package de.adrianwilke.acotspjava.test;

import java.util.LinkedList;
import java.util.List;

import de.adrianwilke.acotspjava.AcoTsp;

public class Test {

    private static final String TEST_FILE = "tsp/d1291.tsp";

    private static List<String> argList = new LinkedList<String>();

    public static void main(String[] args) {

	// System.out.println(new File("").getAbsolutePath());

	argList.add("-h");

	// putTestParameters();
	// argList.add("-l");
	// argList.add("0");

	// ./AcoTsp --ants 20 --time 2 --tours 2 --tries 2 -i ../../d1291.tsp -l 0 --as

	// Comments: ACOTSP, ACOTSPJava

	// argList.add("--as");
	// -l 0: 62571,67053
	// -l 1: 52842,53925
	// -l 2: 52547,53464
	// -l 3: 51759,52609

	// argList.add("--eas");
	// -l 0: 84079,90539
	// -l 3: 51524,52335

	// argList.add("--ras");
	// -l 0: 68201,78482
	// -l 3: 51401,51890

	// argList.add("--mmas");
	// -l 0: 97732,98702
	// -l 3: 52179,52698

	// argList.add("--bwas");
	// -l 0:92378,94083
	// -l 3:51076,52068

	// argList.add("--acs");
	// -l 0:61881,63012
	// -l 3:50962,51401

	AcoTsp.main(argList.toArray(new String[0]));
    }

    private static void putTestParameters() {
	argList.add("--ants");
	argList.add("20");

	argList.add("--time");
	argList.add("2");

	argList.add("--tours");
	argList.add("2");

	argList.add("--tries");
	argList.add("2");

	argList.add("-i");
	argList.add(TEST_FILE);
    }

}
