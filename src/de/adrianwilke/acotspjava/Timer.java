package de.adrianwilke.acotspjava;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 */
public class Timer {
    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: unix_timer.c
     * Author: Thomas Stuetzle
     * Purpose: routines for measuring elapsed time (CPU or real)
     * Check: README.txt and legal.txt
     */

    private static long startTime;

    static void start_timers()
    /*
     * FUNCTION: virtual and real time of day are computed and stored to
     * allow at later time the computation of the elapsed time
     * (virtual or real)
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: virtual and real time are computed
     */
    {
	startTime = System.currentTimeMillis();
    }

    static double elapsed_time()
    /*
     * FUNCTION: return the time used in seconds (virtual or real, depending on type)
     * INPUT: TUtilities.IMER_TYPE (virtual or real time)
     * OUTPUT: seconds since last call to start_timers (virtual or real)
     * (SIDE)EFFECTS: none
     */
    {
	return (System.currentTimeMillis() - startTime) / 1000.0;
    }

}
