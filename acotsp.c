/*

  AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
 AA  AA  CC     OO  OO    TT    SS      PP  PP
 AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
 AA  AA  CC     OO  OO    TT        SS  PP
 AA  AA   CCCC   OOOO     TT    SSSSS   PP
 
 
 MM   MM  PPPPP   IIIIII
 MMMMMMM  PP  PP    II
 MM M MM  PPPPP     II
 MM   MM  PP        II
 MM   MM  PP      IIIIII
	
 
  OOOO    MM   MM  PPPPP
 OO  OO   MMMMMMM  PP  PP
 OO  OO   MM M MM  PPPPP
 OO  OO   MM   MM  PP
  OOOO    MM   MM  PP
 
 
 
 ######################################################
 ##########    MPI+OMP PARALLEL              ##########
 ##########    ACO algorithms for the TSP    ##########
 ######################################################

      Version: 1.0
      File:    acotsp.c
      Purpose: main routines and control for the ACO algorithms
      Check:   README and gpl.txt
*/


#include <mpi.h>
//#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>

#include "parallel.h"
#include "ants.h"
#include "utilities.h"
#include "InOut.h"
#include "TSP.h"
#include "timer.h"
#include "ls.h"
#include "ft_aco.h"



long int termination_condition( void )
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
    stopColonies = 0;

    printf("[termination_condition] Process %d / %d -> iterations = %d\n", mpi_id, NPROC, iteration);
    if (NPROC == 1) 
    	best_global_tour_length = best_so_far_ant -> tour_length;

    return ( elapsed_time( REAL ) >= max_time ||
            iteration > max_iters || best_global_tour_length <= optimal );
    //return best_global_tour_length <= optimal;
    }



void construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    TRACE ( printf("construct solutions for all ants\n"); );

    /***  OMP PARALLEL LOOP ****/
    for ( k = 0 ; k < n_ants ; k++) {
        step=0;
        /* Mark all cities as unvisited */
        ant_empty_memory( &ant[k] );
        /* Place the ants on same initial city */
        place_ant( &ant[k], step);

        while ( step < n-1 ) {
	    step++;
	    neighbour_choose_and_move_to_next( &ant[k], step);

        }

        step = n;
        ant[k].tour[n] = ant[k].tour[0];
        ant[k].tour_length = compute_tour_length( ant[k].tour );

    }

    n_tours += n_ants;
}



void init_try( long int ntry ) 
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
{

    TRACE ( printf("INITIALIZE TRIAL\n"); );

  //  MPI_Barrier(comm);

    start_timers();
    time_used = elapsed_time( REAL );
    time_passed = time_used;

    if (comp_report) {
        fprintf(comp_report,"seed %ld\n",seed);
        fflush(comp_report);
    }
    
    /* Initialize variables concerning statistics etc. */
  
    n_tours      = 0;
    iteration    = 1;
    restart_iteration = 1;
    lambda       = 0.05;            
    best_so_far_ant->tour_length = INFTY;
    best_global_tour_length = INFTY;
    found_best   = 0;
  
    /* Initialize the Pheromone trails */
    if ( !(mmas_flag || bwas_flag) ) {
	trail_0 = 1. / ( (rho) * nn_tour() );
	/* in the original papers on Ant System, Elitist Ant System, and
	   Rank-based Ant System it is not exactly defined what the
	   initial value of the pheromones is. Here we set it to some
	   small constant, analogously as done in MAX-MIN Ant System.  
	*/
	init_pheromone_trails( trail_0 );
    } 
    if ( bwas_flag ) {
	trail_0 = 1. / ( (double) n * (double) nn_tour()) ;
	init_pheromone_trails( trail_0 );
    } 
    if ( mmas_flag ) {
	trail_max = 1. / ( (rho) * nn_tour() );
	trail_min = trail_max / ( 2. * n );
	init_pheromone_trails( trail_max );   
    }
  
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();
    
    if (comp_report) fprintf(comp_report,"\nbegin try %li ...\n",ntry);
    if (cc_report) fprintf(cc_report,"\nbegin try %li ...\n",ntry);
    if (stat_report) fprintf(stat_report,"\nbegin try %li ...\n",ntry);
}



void local_search( void )
/*    
      FUNCTION:       manage the local search phase; apply local search to ALL ants; in 
                      dependence of ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
		      is chosen.
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants of the colony have locally optimal tours
      COMMENTS:       typically, best performance is obtained by applying local search 
                      to all ants. It is known that some improvements (e.g. convergence 
		      speed towards high quality solutions) may be obtained for some 
		      ACO algorithms by applying local search to only some of the ants.
		      Overall best performance is typcially obtained by using 3-opt.
*/
{
    long int k;

    TRACE ( printf("apply local search to all ants\n"); );

    /****   OMP PARALLEL LOOP ****/
    for ( k = 0 ; k < n_ants ; k++ ) {
	switch (ls_flag) {
        case 1:
	    two_opt_first( ant[k].tour );    /* 2-opt local search */
            break;
        case 2:
	    two_h_opt_first( ant[k].tour );  /* 2.5-opt local search */
            break;
        case 3:
	    three_opt_first( ant[k].tour );  /* 3-opt local search */
            break;
        default:
	    fprintf(stderr,"type of local search procedure not correctly specified\n");
	    exit(1);
	}
	ant[k].tour_length = compute_tour_length( ant[k].tour );
 }
 
}



void update_statistics(MPI_Comm comm)
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    long int iteration_best_ant;
    double p_x; /* only used by MMAS */

    iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */

    if ( ant[iteration_best_ant].tour_length < best_so_far_ant->tour_length ) {

	time_used = elapsed_time( REAL ); /* best sol found after time_used */
	copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
	copy_from_to( &ant[iteration_best_ant], restart_best_ant );

    /*** Asynchronous communication to other colonies ***/
    if (best_so_far_ant->tour_length < best_global_tour_length) {
        sendBestSolutionToColonies(comm);
        printf("[update_statistics] %d / %d: sent best solution\n", mpi_id, NPROC);
    }

	found_best = iteration;
	restart_found_best = iteration;
	found_branching = node_branching(lambda);
	branching_factor = found_branching;
	if ( mmas_flag ) {
	    if ( !ls_flag ) {
		p_x = exp(log(0.05)/n); 
		trail_min = 1. * (1. - p_x) / (p_x * (double)((nn_ants + 1) / 2));
		trail_max = 1. / ( (rho) * best_so_far_ant->tour_length );
		trail_0 = trail_max;
		trail_min = trail_max * trail_min; 
	    } else {
		trail_max = 1. / ( (rho) * best_so_far_ant->tour_length );
		trail_min = trail_max / ( 2. * n );
		trail_0 = trail_max;
	    }
	}
	write_report();
    }
    
    if ( ant[iteration_best_ant].tour_length < restart_best_ant->tour_length ) {
	copy_from_to( &ant[iteration_best_ant], restart_best_ant );
	restart_found_best = iteration;
	printf("restart best: %ld, restart_found_best %ld, time %.2f\n",restart_best_ant->tour_length, restart_found_best, elapsed_time ( REAL )); 
    }

    /*** Listen from other colonies ***/
    printf("[update_statistics] %d / %d: before 'listenTours'\n", mpi_id, NPROC);
	listenTours(comm);
    printf("[update_statistics] %d / %d: after 'listenTours'\n", mpi_id, NPROC);

 
}


void search_control_and_statistics( void )
/*    
      FUNCTION:       occasionally compute some statistics and check whether algorithm 
                      is converged 
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{
    TRACE ( printf("SEARCH CONTROL AND STATISTICS\n"); );

    if (!(iteration % 100)) {
	population_statistics();
	branching_factor = node_branching(lambda);
        
    }

	if ( mmas_flag && (branching_factor < branch_fac) && (iteration - restart_found_best > 250) ) {
	    /* MAX-MIN Ant System was the first ACO algorithm to use
	       pheromone trail re-initialisation as implemented
	       here. Other ACO algorithms may also profit from this mechanism.
	    */
	    printf("INIT TRAILS!!!\n"); restart_best_ant->tour_length = INFTY; 
	    init_pheromone_trails( trail_max );
	    compute_total_information();
	    restart_iteration = iteration;
	    best_global_tour_length = INFTY;
	    restart_time = elapsed_time( REAL );
	}
}




void as_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
    
    if (NPROC > 1)
            foreign_solution_update_pheromone( best_global_tour );
    

}



void eas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Elitist Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants plus elitist ant deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Elitist Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
    global_update_pheromone_weighted( best_so_far_ant, elitist_ants );
    
    if (NPROC > 1) {
        if (best_global_tour_length < best_so_far_ant->tour_length) {
            foreign_solution_update_pheromone_weighted( best_global_tour, elitist_ants);
        } else {
            foreign_solution_update_pheromone( best_global_tour );
        }
    }

}



void ras_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    long int i, k, b, target;
    long int *help_b;

    TRACE ( printf("Rank-based Ant System pheromone deposit\n"); );

    help_b = malloc( n_ants  * sizeof(long int) );
    for ( k = 0 ; k < n_ants ; k++ )
	help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ ) {
	b = help_b[0]; target = 0;
	for ( k = 0 ; k < n_ants ; k++ ) {
	    if ( help_b[k] < b ) {
		b = help_b[k]; target = k;
	    }
	}
	help_b[target] = LONG_MAX;
	global_update_pheromone_weighted( &ant[target], ras_ranks-i-1 );
    }
    global_update_pheromone_weighted( best_so_far_ant, ras_ranks );
    
    if (NPROC > 1)
            foreign_solution_update_pheromone_weighted( best_global_tour, ras_ranks );
    
    free ( help_b );
}



void mmas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for MAX-MIN Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    /* we use default upper pheromone trail limit for MMAS and hence we
       do not have to worry regarding keeping the upper limit */

    long int iteration_best_ant;

    TRACE ( printf("MAX-MIN Ant System pheromone deposit\n"); );

    if ( iteration % u_gb ) {
	iteration_best_ant = find_best();
	global_update_pheromone( &ant[iteration_best_ant] );
    }
    else {
        if ( u_gb == 1 && (iteration - restart_found_best > 50))
	    global_update_pheromone( best_so_far_ant );
        else 
	    global_update_pheromone( restart_best_ant );
    }

    if (NPROC > 1 && (best_global_tour_length < best_so_far_ant->tour_length))
            foreign_solution_update_pheromone( best_global_tour );
    

    
    if ( ls_flag ) {
	/* implement the schedule for u_gb as defined in the 
	   Future Generation Computer Systems article or in Stuetzle's PhD thesis.
	   This schedule is only applied if local search is used.
	*/
	if ( ( iteration - restart_iteration ) < 25 )
	    u_gb = 25;
	else if ( (iteration - restart_iteration) < 75 )
	    u_gb = 5;
	else if ( (iteration - restart_iteration) < 125 )
	    u_gb = 3;
	else if ( (iteration - restart_iteration) < 250 )
	    u_gb = 2;
	else 
	    u_gb = 1;
    } else
	u_gb = 25;
  
}



void bwas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Best-Worst Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    long int   iteration_worst_ant, distance_best_worst;

    TRACE ( printf("Best-worst Ant System pheromone deposit\n"); );

    global_update_pheromone( best_so_far_ant );
    if (NPROC > 1)
            foreign_solution_update_pheromone( best_global_tour );

    iteration_worst_ant = find_worst();
    bwas_worst_ant_update( &ant[iteration_worst_ant], best_so_far_ant );
    distance_best_worst = distance_between_ants( best_so_far_ant, &ant[iteration_worst_ant]);
/*    printf("distance_best_worst %ld, tour length worst %ld\n",distance_best_worst,ant[iteration_worst_ant].tour_length); */
    if ( distance_best_worst < (long int) (0.05 * (double) n) ) {
	restart_best_ant->tour_length = INFTY;
	init_pheromone_trails( trail_0 );
	restart_iteration = iteration;    
	restart_time = elapsed_time( REAL );
	printf("init pheromone trails with %.15f, iteration %ld\n",trail_0,iteration);
    }
    else 
	bwas_pheromone_mutation();
}


void pheromone_trail_update( void )  
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
                      according to the rules defined by the various ACO algorithms.
*/
{
    /* Simulate the pheromone evaporation of all pheromones */
    if ( as_flag || eas_flag || ras_flag || bwas_flag || mmas_flag ) {
	if ( ls_flag ) {
	    if ( mmas_flag )
		mmas_evaporation_nn_list();
	    else
		evaporation_nn_list();
	    /* evaporate only pheromones on arcs of candidate list to make the 
	       pheromone evaporation faster for being able to tackle large TSP 
	       instances. For MMAS additionally check lower pheromone trail limits.
	    */
	} else {
	    /* if no local search is used, evaporate all pheromone trails */
	    evaporation();
	}
    }

    /* Next, apply the pheromone deposit for the various ACO algorithms */
    if ( as_flag )
	as_update(); 
    else if ( eas_flag )
	eas_update();
    else if ( ras_flag )
	ras_update();
    else if ( mmas_flag )
	mmas_update();
    else if ( bwas_flag )
	bwas_update();
   

  /* check pheromone trail limits for MMAS; not necessary if local
     search is used, because in the local search case lower pheromone trail
     limits are checked in procedure mmas_evaporation_nn_list */
    if ( mmas_flag && !ls_flag )
	check_pheromone_trail_limits();

  /* Compute combined information pheromone times heuristic info after
     the pheromone update for ACO algorithms */
    if ( as_flag || eas_flag || ras_flag || mmas_flag || bwas_flag ) {
	if ( ls_flag ) {
	    compute_nn_list_total_information();
	} else {
	    compute_total_information();
	}
    }
}

#if defined FT_ACO && defined FT_ABORT_ON_FAILURE
void cleanup(void * ptr)
{
    printf("CLEANUP FUNCTION %d\n", *((int *) ptr));
}
#endif

/* --- main program ------------------------------------------------------ */

int main(int argc, char *argv[]) {
/*    
      FUNCTION:       main control for running the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
      COMMENTS:       this function controls the run of "max_tries" independent trials

*/

    MPI_Comm comm;
    #if  defined FT_ACO && (defined FT_ERRORS_ARE_FATAL || defined FT_ERRORS_RETURN || defined FT_ABORT_ON_FAILURE  || defined FT_IGNORE_ON_FAILURE)
    MPI_Errhandler error_handler;
    #ifdef FT_ABORT_ON_FAILURE
    int ft_parameter = 666;
    #endif 
    #endif
    long int i = 1;
    start_timers();

    int provided;
    int isNewSpawnee = 0;
    int killFlag = 0;
    /** MPI Initialization **/

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided == MPI_THREAD_MULTIPLE) {
        printf("Executing in MPI_THREAD_MULTIPLE mode\n");
    }

	FT_init();

    //MPI_Init(&argc, &argv);
    // handles creation by checking parent 
    MPI_Comm parent;
    MPI_Comm_get_parent(&parent);
    if (parent == MPI_COMM_NULL) {
        printf("[main loop] original process\n");
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    } else {
        printf("[main loop new spawnee] new spawnee after failure\n");
        printf("[main loop new spawnee] Trying to attach to quasi-global communicator...\n");
        FT_set_respawn_data(argv, &mpi_id, &NPROC);
        int stt = attach_to_comm(&comm);
        printf("[main loop new spawnee] status: %d\n", stt);
        printf("[main loop new spawnee] new comm\n");
        if (MPI_COMM_NULL == comm) {
            printf("[main loop new spawnee] MPI_COMM_NULL\n");
        } else {
            printf("[main loop new spawnee] not  MPI_COMM_NULL\n");
        }
        printf("[main loop spawnee] before comm test\n");
        MPI_Comm_rank(comm, &mpi_id);
        MPI_Comm_size(comm, &NPROC);
        printf("[main loop spawnee] after comm test, %d / %d\n", mpi_id, NPROC);
        MPI_Comm_set_errhandler(comm, FT_RESPAWN_ON_FAILURE_HANDLER);
        printf("[main loop spawnee]: Using FT_RESPAWN_ON_FAILURE_HANDLER\n");
        isNewSpawnee = 1;
    }
    if (isNewSpawnee == 1)
        printf("[main loop new spawnee] computing initial size\n");
    MPI_Comm_size(comm, &NPROC);
    if (isNewSpawnee == 1)
        printf("[main loop new spawnee] computing initial rank\n");
    MPI_Comm_rank(comm, &mpi_id);  

    if (isNewSpawnee == 1)
        printf("[main loop new spawnee] %d / %d: initial parameters computed\n", mpi_id, NPROC);


    #ifdef FT_ACO
    #ifdef FT_ERRORS_ARE_FATAL
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_ARE_FATAL);
    printf("Using FT_ERRORS_ARE_FATAL_HANDLER\n");
    #endif
    #ifdef FT_ERRORS_RETURN
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    printf("Using FT_ERRORS_RETURN_HANLDER\n");
    #endif
    #ifdef FT_ABORT_ON_FAILURE
    MPI_Comm_set_errhandler(comm, FT_ABORT_ON_FAILURE_HANDLER);
    FT_set_cleanup_function(cleanup, (void *) &ft_parameter);
    printf("Using FT_ABORT_ON_FAILURE_HANDLER\n");
    #endif
    #ifdef FT_IGNORE_ON_FAILURE
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    MPI_Comm_set_errhandler(comm, FT_IGNORE_ON_FAILURE_HANDLER);
    printf("[%s]: Using FT_IGNORE_ON_FAILURE_HANDLER\n", get_current_time(time));
    FT_set_respawn_data(argv, &mpi_id, &NPROC);
    printf("[%s] Respawn data shared\n", get_current_time(time));
    #endif
    #ifdef FT_RESPAWN_ON_FAILURE
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    MPI_Comm_set_errhandler(comm, FT_RESPAWN_ON_FAILURE_HANDLER);
    printf("[%s]: Using FT_RESPAWN_ON_FAILURE_HANDLER\n", get_current_time(time));
    FT_set_respawn_data(argv, &mpi_id, &NPROC);
    printf("[%s] Respawn data shared\n", get_current_time(time));
    #endif
    #ifdef FT_ELASTIC_RESPAWN_ON_FAILURE
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    MPI_Comm_set_errhandler(comm, FT_ELASTIC_RESPAWN_ON_FAILURE_HANDLER);
    printf("[%s]: Using FT_ELASTIC_RESPAWN_ON_FAILURE_HANDLER\n", get_current_time(time));
    FT_set_elastic_respawn_data(argv, &mpi_id, &NPROC, 1, 4);
    printf("[%s] Respawn data shared\n", get_current_time(time));
    #endif

    //MPI_Barrier(comm);
    printf("Process %d / %d: reporting alive\n", mpi_id, NPROC);
    #endif

    

    init_program(argc, argv);
    
    instance.nn_list = compute_nn_lists();
    
    pheromone = generate_double_matrix( n, n );
    total = generate_double_matrix( n, n );

    time_used = elapsed_time( REAL );
    printf("Initialization took %.10f seconds\n",time_used);

    for ( n_try = 0 ; n_try < max_tries ; n_try++ ) {

        init_try(n_try);

        if ( NPROC > 1 ) {
            startCommColoniesTour(comm); /*prepare buffer for communications from Colonies */
        }

        while ( !termination_condition() ) {
            
            //printf("Entering while loop\n");
            construct_solutions();
            
            printf("[main loop] %d / %d: constructing solutions\n", mpi_id, NPROC);
            if ( ls_flag > 0 )
                local_search();

            printf("[main loop] %d / %d: updating statistics\n", mpi_id, NPROC);
            update_statistics(comm);

            printf("[main loop] %d / %d: pheromone trail update\n", mpi_id, NPROC);
            pheromone_trail_update();  

            printf("[main loop] %d / %d: search control & statistics\n", mpi_id, NPROC);
            search_control_and_statistics();


            #ifdef FT_ACO_ASDFOJ

            printf("[main loop] %d / %d TRY: %d\n", mpi_id, NPROC, n_try);
            if (killFlag == 0 && NPROC > 1 && mpi_id == NPROC-1) {
                printf("[main loop] %d / %d: old MÁTOME AQUÍ\n", mpi_id, NPROC);
                raise(SIGKILL);
            }
            killFlag = 0;


            //MPI_Bcast(&NPROC, 1, MPI_INT, 0, comm);
            //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
            //MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
            //printf("AFTER FAILURE REW PROCESS MAPPING: %d / %d\n", mpi_id, NPROC);
            #endif
            #ifdef KILL_AT_50
            if (iteration == 5 && mpi_id == NPROC - 1) {
                printf("[main loop] %d / %d: MÁTOME AQUÍ\n", mpi_id, NPROC);
                raise(SIGKILL);
            }
            #endif
            #ifdef KILL_AT_25_AND_50
            if (iteration == 250 && mpi_id == 4) {
                printf("[main loop] 25% %d / %d: MÁTOME AQUÍ\n", mpi_id, NPROC);
                raise(SIGKILL);
            }            
            if (iteration == 500 && mpi_id == 3) {
                printf("[main loop] 50% %d / %d: MÁTOME AQUÍ\n", mpi_id, NPROC);
                raise(SIGKILL);
            }            
            #endif

            #ifdef KILL_AT_25_AND_50
            /*
            printf("[main loop] kill at 25 and 50\n");
            if (iteration == 250 && mpi_id == 15) {
                printf("[main loop] %d / %d: MÁTOME AQUÍ na iteración %d; try = %d\n", mpi_id, NPROC, iteration, n_try);
                raise(SIGKILL);
            }
            if (iteration == 500 && mpi_id == 14) {
                printf("[main loop] %d / %d: MÁTOME AQUÍ na iteración %d; try = %d\n", mpi_id, NPROC, iteration, n_try);
                raise(SIGKILL);
            }
            */
            #endif
            #ifdef KILL_ALL_BUT_ONE
            printf("kill all but one, id = %d\n", mpi_id);
            if (iteration % 16 == 0) {
                if (mpi_id != 0 && mpi_id == NPROC - i) {
                    //i++;
                    printf("[main loop] %d / %d: MÁTOME AQUÍ iteration = %d\n", mpi_id, NPROC, iteration);
                    raise(SIGKILL);
                }
                i++;
            }
            #endif
            iteration++;
        }

        printf("[main loop] about to exit try\n");
        exit_try(comm, n_try);
        printf("[main loop] try exited\n");
            
    }
    
    exit_program();

    free( instance.distance );
    free( instance.nn_list );
    free( pheromone );
    free( total );
    free( best_in_try );
    free( best_found_at );
    free( time_best_found );
    free( time_total_run );
    for ( i = 0 ; i < n_ants ; i++ ) {
	free( ant[i].tour );
	free( ant[i].visited );
    }
    free( ant );
    free( best_so_far_ant->tour );
    free( best_so_far_ant->visited );
    free( prob_of_selection );

    /** MPI finalize **/
    fflush(parallel);
    fflush(cc_report);
    fflush(report);
    fflush(comp_report);
    fflush(stat_report);
    
    #ifdef FT_ACO
    //MPI_Errhandler_free(&error_handler);
    FT_finalize();
    #endif


    //MPI_Barrier(comm);
    MPI_Finalize();

    return(0);
}
