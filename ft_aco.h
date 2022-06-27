/*
 * Fault tolerant ACO
 * Author: Miguel Blanco Godón
 * Computer Engineering, 2022
 */

#ifndef __FT_ACO_H
#define __FT_ACO_H

#include <mpi.h>
#include <mpi-ext.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "parallel.h"

#define FT_FAILURE 1
#define FT_SUCCESS 0
#define SHARE_COMM_TAG 0x00FEDCBA

// PREDEFINED ERROR HANDLERS

extern MPI_Errhandler ft_abort_on_failure_error_handler;
extern MPI_Errhandler ft_ignore_on_failure_error_handler;

/*
 * FT_ERRORS_ARE_FATAL_ON_FAILURE_HANDLER
 * This handler is a wrapper to the default MPI_ERRORS_ARE_FATAL handler.
 * The usage of this handler causes abnormal termination of the MPI program 
 * whenever an error is detected.
 */
#define FT_ERRORS_ARE_FATAL_ON_FAILURE_HANDLER MPI_ERRORS_ARE_FATAL

/*
 * FT_ERRORS_RETURN_ON_FAILURE_HANDLER
 * This handler is a wrapper to the default MPI_ERRORS_RETURN handler.
 * The usage of this handler causes MPI functions to return an error code
 * when an error is detected, and the execution continues without any 
 * further error control. From this call on, deterministic behavior is not 
 * guaranteed.
 */
#define FT_ERRORS_RETURN_ON_FAILURE_HANDLER MPI_ERRORS_RETURN

/*
 * FT_ABORT_ON_FAILURE_HANDLER
 * When it is called, this handler tries to finalize the MPI program in the 
 * most harmless way. A "clean up" function can be defined 
 * (see FT_set_cleanup_funciton for further details) in order to finalize 
 * the MPI program storing some data and restoring dynamic memory.
 * At the end of the routine, MPI_Finalize() is called.
 */
#define FT_ABORT_ON_FAILURE_HANDLER ft_abort_on_failure_error_handler

/*
 * FT_IGNORE_ON_FAILURE_HANDLER
 * When it is called, this handler tries to continue the MPI program by 
 * ignoring the MPI errors. In order to do so, it remplaces the current 
 * communicatior with a new communicator involving the processes who 
 * are reported alive whenever the handler is invoked.
 */
#define FT_IGNORE_ON_FAILURE_HANDLER ft_ignore_on_failure_error_handler


/*
 * Initializes predefined error handlers.
 * INPUT: none.
 * OUTPUT: none.
 * PRECONDITIONS: none.
 * SIDE EFFECTS: none.
 */
void FT_init(void);

/*
 * Destroys predefined error handlers. Should be called only at the end of the program.
 * INPUT: none.
 * OUTPUT: none.
 * PRECONDITIONS: previous call to FT_init. Otherwise, behavior is not defined.
 * SIDE EFFECTS: all predefined error handlers are set to null.
 */
void FT_finalize(void);

/*
 * Asociates the error handler 'error_handler' to the communicator 'comm'.
 * INPUT: a comunicatior and an error handler.
 * OUTPUT: none.
 * PRECONDITIONS: none.
 * SIDE_EFFECTS: none.
 */
void FT_set_error_handler(MPI_Comm comm, MPI_Errhandler error_handler);

/*
 * Sets the value of the cleanup function and parameters.
 * INPUT: a pointer to de defined function and a pointer to the function parameters.
 * OUTPUT: none.
 * PRECONDITIONS: none.
 * SIDE EFFECTS: none.
 */
void FT_set_cleanup_function(void * function, void * parameters);

/*
 * Sets the value of the cleanup function parameters.
 * INPUT: a pointer to the struct containing the paramters.
 * OUTPUT: FT_SUCCESS on success, FT_FAILURE on error.
 * PRECONDITIONS: the cleanup function cannot be null.
 * SIDE EFFECTS: none.
 */
int FT_set_cleanup_params(void * parameters);

/*
 * Finalizes all processes safely.
 * If there is a "cleanup" function defined in the scope, it will be called before 
 * the program is finalized.
 * MPI_Finalize is called at the end of the routine.
 * INPUT: communicator to the processes and an error set.
 * OUTPUT: none.
 * PRECONDITIONS: none
 * SIDE EFFECTS: finalizes MPI execution.
 */
void FT_abort_on_failure(MPI_Comm * comm, int * err, ...);

/*
 * Catches an error and continues without further control. 
 * Should be used only if each worker task is completely independent from others,
 * thus, not using further communication.
 * INPUT: communicator to the processes and an error set.
 * OUTPUT: none.
 * PRECONDITIONS: none.
 * SIDE EFFECTS: none.
 */
void FT_ignore_on_failure(MPI_Comm * comm, int * err, ...);

void FT_set_respawn_data(char ** argv, int * id, int * nprocs);


char * get_current_time(char * output);
 
#endif
