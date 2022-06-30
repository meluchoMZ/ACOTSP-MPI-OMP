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
 File:    parallel.c
 Purpose: routines added in the parallel implementation
 */


#include <signal.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "TSP.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"
#include "timer.h"
#include "parallel.h"

MPI_Request request[50][100]; //request for a maximum of 50 colonies and 100 tries
MPI_Request SRrequest;
MPI_Request Srequest;
MPI_Status status;
MPI_Status Sstatus;

int mpi_id;
int NPROC;
int stopColonies;
long int ss;

void parallel_init(MPI_Comm comm, MPI_Errhandler error_handler)
{
    MPI_Init(NULL,NULL);
 
    MPI_Comm_rank(comm, &mpi_id);
    MPI_Comm_size(comm, &NPROC);
    
    // avoid errors to be fatal, depends on MPI framework
    MPI_Comm_set_errhandler(comm, error_handler);
}


/********** Asynchronous communication protocol *****************/

/**
 Routine to start recv. of communications from Colonies
 **/

void startCommColoniesTour (MPI_Comm comm)
{
    /*
    int i;
    
    /* Prepare to receive the best solutions found in
     the rest of the colonies */
    /*
    for (i=0 ; i<NPROC ; i++)
        if (i!=mpi_id)
   	    	MPI_Irecv(&global_tour[i][0], n+1, MPI_LONG, i,
        	      3000, comm, &request[i][n_try]);
    */
}


/**
 Routine to send the best solution found to the rest of Colonies
 **/
void sendBestSolutionToColonies (MPI_Comm comm)
{
    int     i,j;
    int rank, size;
    int tag = 3000 + n_try;

    /* Send best solution found (tour) to the rest of the colonies */
    printf("[sendBestSolutionToColonies] %d / %d entering method\n", mpi_id, NPROC);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    printf("[sendBestSolutionToColonies] global: %d / %d; local: %d / %d\n", mpi_id, NPROC, rank, size);
    for( i=0 ; i<NPROC ; i++ )
       if ( i != mpi_id){
           MPI_Isend(best_so_far_ant->tour, n+1, MPI_LONG, i,
                     3000, comm, &Srequest);
           printf("[sendBestSolutionToColonies] %d / %d sending best solution to peer %d / %d\n", mpi_id, NPROC, i, NPROC);
        }

    for ( j=0 ; j<n+1 ; j++ )
            best_global_tour[j] = best_so_far_ant->tour[j];
    best_global_tour_length = best_so_far_ant->tour_length;
    
    if (mpi_id==0 && cc_report) fprintf(cc_report,"%ld \t %f\n",best_global_tour_length,elapsed_time(REAL));
    printf("[sendBestSolutionToColonies] Process %d / %d: sent best solution to colonies\n", rank, size);

}


/**
 Routine to check if there are pending Tour messages from Colonies
 **/
void listenTours(MPI_Comm comm)
{
    int     i, j;
    int     flag;
    long int    recv_tour_length;
    int tag = 3000 + n_try;
    
    /* Loop to listen best solutions from colonies */
    //for(i=0; i<NPROC ; i++) {
 
      //if(i!=mpi_id){

        flag = 1;
        printf("[listenTours] %d / %d: entering\n", mpi_id, NPROC);

        while ( flag == 1 ) {

	        //MPI_Test(&request[i][n_try], &flag, &status);
            printf("[listenTours] %d / %d before i-probe 1\n", mpi_id, NPROC);
            printf("[listenTours] %d / %d before i-probe 2\n", mpi_id, NPROC);
            printf("[listenTours] %d / %d before i-probe 3\n", mpi_id, NPROC);
            printf("[listenTours] %d / %d before i-probe 4\n", mpi_id, NPROC);
            printf("[listenTours] %d / %d before i-probe 5\n", mpi_id, NPROC);
            MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &flag, &status);
            printf("[listenTours] %d / %d i-probe with flag = %d received from %d / %d with tag=%d\n", mpi_id, NPROC, flag, status.MPI_SOURCE, NPROC, status.MPI_TAG);
        
            if (flag==1) {
                MPI_Recv(&global_tour[status.MPI_SOURCE][0], n + 1, MPI_LONG, status.MPI_SOURCE, tag, comm, &status);
                printf("[listenTours] %d / %d received message from %d / %d with tag=%d", mpi_id, NPROC, status.MPI_SOURCE, NPROC, status.MPI_TAG);

              recv_tour_length = compute_tour_length( &global_tour[status.MPI_SOURCE][0] );
            
              if ( recv_tour_length < best_global_tour_length) {
                best_global_tour_length = recv_tour_length;
                for ( j=0 ; j<n+1 ; j++ )
                    best_global_tour[j] = global_tour[status.MPI_SOURCE][j];
                  
                if (mpi_id==0 && cc_report) fprintf(cc_report,"%ld \t %f\n",best_global_tour_length,elapsed_time(REAL));
                write_report();
              }
                
            /* Prepare to receive more solutions */
               MPI_Irecv(&global_tour[status.MPI_SOURCE][0], n+1, MPI_LONG,
                      status.MPI_SOURCE, 3000, comm, &request[i][n_try]);
            
            }
     // }
    // }
   }
    
}


/********************* Update pheromone **********************************/

/**
 Routines to update pheromone using the foreigner tour
 **/
void foreign_solution_update_pheromone( long int *ftour )
/*
 FUNCTION:      reinforces edges used in foreigner solution
 INPUT:         pointer to foreigner tour that updates the pheromone trail
 OUTPUT:        none
 (SIDE)EFFECTS: pheromones of arcs in foreigner tour are increased
 */
{
    long int i, j, h, ftour_length;
    double   d_tau;
    
    TRACE ( printf("global pheromone update with foreigner solutions\n"); );
    
    ftour_length = compute_tour_length( ftour );
    d_tau = 1.0 / (double) ftour_length;
    for ( i = 0 ; i < n ; i++ ) {
        j = ftour[i];
        h = ftour[i+1];
        pheromone[j][h] += d_tau;
        pheromone[h][j] = pheromone[j][h];
    }
}

void foreign_solution_update_pheromone_weighted( long int *ftour, long int weight )
/*
 FUNCTION:      reinforces edges used in foreigner solution with weight
 INPUT:         pointer to foreigner tour that updates the pheromone trail and its weight
 OUTPUT:        none
 (SIDE)EFFECTS: pheromones of arcs in foreigner tour are increased
 */
{
    long int      i, j, h, ftour_length;
    double        d_tau;
    
    TRACE ( printf("global pheromone update with foreigner solutions weighted\n"); );
    
    ftour_length = compute_tour_length( ftour );
    d_tau = (double) weight / (double) ftour_length;
    for ( i = 0 ; i < n ; i++ ) {
        j = ftour[i];
        h = ftour[i+1];
        pheromone[j][h] += d_tau;
        pheromone[h][j] = pheromone[j][h];
    }       
}


/***************************  Final MPI-Report ************************************/

/**
 Routine to write one summarize mpi report
 **/
void write_mpi_report (MPI_Comm comm)
{
    int     i, com_id;
    long int    best_com_tour, com_tour_length, com_iter, com_found_best;
    double  com_found_branching, com_b_fac, com_time, com_best_time;
    
    best_com_tour = best_so_far_ant->tour_length;
    com_found_best = found_best;
    com_found_branching = found_branching;
    com_best_time = time_used;
    com_id = mpi_id;
    
    if( mpi_id == 0 ) {
        for( i=1 ; i<NPROC ; i++ ) {
              MPI_Recv(&com_tour_length, 1, MPI_LONG, i, 2000, comm, &status);
              MPI_Recv(&com_iter, 1, MPI_LONG, i, 2000, comm, &status);
              MPI_Recv(&com_b_fac, 1, MPI_DOUBLE, i, 2000, comm, &status);
              MPI_Recv(&com_time, 1, MPI_DOUBLE, i, 2000, comm, &status);
            
              if ( com_tour_length < best_com_tour ) {
                best_com_tour = com_tour_length;
                com_found_best = com_iter;
                com_found_branching = com_b_fac;
                com_best_time = com_time;
                com_id=i;
              }
	    
        }
        fprintf(summary, "Try: %ld\t Best: %ld\t Iterations: %6ld\t B-Fac %.5f\t Time %.2f\t Tot.time %.2f\t PROC: %d\n",
                n_try,best_com_tour, com_found_best, com_found_branching,
                com_best_time, elapsed_time( REAL ), com_id);
        fflush(summary);
        fflush(cc_report);
    } else {
        MPI_Isend(&best_com_tour, 1, MPI_LONG, 0, 2000, comm,&SRrequest);
        MPI_Isend(&com_found_best, 1, MPI_LONG, 0, 2000, comm,&SRrequest);
        MPI_Isend(&com_found_branching, 1, MPI_DOUBLE, 0, 2000, comm,&SRrequest);
        MPI_Isend(&com_best_time, 1, MPI_DOUBLE, 0, 2000, comm,&SRrequest);
    }
    

}

