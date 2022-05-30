/*
 * Fault tolerant ACO
 * Author: Miguel Blanco Godón
 * Computer Engineering, 2022
 */

#include "ft_aco.h"

/*
 * Predefined error handlers. They are initialized on FT_init procedure
 * and destroyed on FT_finalize procedure.
 */
MPI_Errhandler ft_abort_on_failure_error_handler = NULL;
MPI_Errhandler ft_ignore_on_failure_error_handler = NULL;

/*
 * Cleanup function.
 * This function is called before finilizing MPI on an error.
 * The value of this function can be set with the primitive FT_Set_Cleanup_Funcion.
 * INPUT: a struct containing the necessary information needed to end the program
 * without too many collateral damage (freeing memory and saving data).
 * OUTPUT: none.
 * PRECONDITIONS: none.
 * CONSECUENCES: potential program state could be lost.
 */
static void (* ft_cleanup_function) (void *) = NULL;

/*
 * Stores data needed for cleanup function to work properly.
 * This value can be set when calling to FT_set_cleanup_function or when calling 
 * to FT_set_cleanup_params.
 */
static void * ft_cleanup_params = NULL;

void FT_init(void)
{
    MPI_Comm_create_errhandler(&FT_abort_on_failure, &ft_abort_on_failure_error_handler);
    MPI_Comm_create_errhandler(&FT_ignore_on_failure ,&ft_ignore_on_failure_error_handler);
}

void FT_finalize(void)
{
    MPI_Errhandler_free(&ft_abort_on_failure_error_handler);
    MPI_Errhandler_free(&ft_ignore_on_failure_error_handler);
}

void FT_set_error_handler(MPI_Comm comm, MPI_Errhandler error_handler)
{
    MPI_Comm_set_errhandler(comm, error_handler);
}

void FT_set_cleanup_function(void * function, void * parameters)
{
    ft_cleanup_function = function; 
    ft_cleanup_params = parameters;
}

int FT_set_cleanup_params(void * parameters)
{
    if (ft_cleanup_function == NULL) {
        return FT_FAILURE;
    }
    ft_cleanup_params = parameters;
    return FT_SUCCESS;
}

void FT_abort_on_failure(MPI_Comm * comm, int * err, ...)
{
    char error_info[MPI_MAX_ERROR_STRING];
    int error_info_length, rank, size;

    MPI_Comm_rank(*comm, &rank);
    MPI_Comm_size(*comm, &size);
    // Retrieves error cause
    MPI_Error_string(*err, error_info, &error_info_length);
    printf("Process %d / %d received error with exception: %s\n", rank, size, error_info);
    
    if (ft_cleanup_function != NULL) {
        printf("Process %d / %d: cleaning up environment...\n", rank, size);
        (* ft_cleanup_function)(ft_cleanup_params);
    }

    printf("Process %d / %d: finalizing MPI execution", rank, size);

    MPI_Finalize();
    exit(EXIT_FAILURE);
}

void FT_ignore_on_failure(MPI_Comm * comm, int * err, ...)
{
    if (comm != NULL) {
        printf("Comm is not null, error: %i\n", *err);
    }
    /*
    char error_info[MPI_MAX_ERROR_STRING];
    int rank, size, error_info_length;
    MPI_Comm * new_comm = malloc(sizeof(MPI_Comm));
    int new_rank, new_size;
    int mpi_call_status = -34132514;
    MPI_Comm_rank(*comm, &rank);
    MPI_Comm_size(*comm, &size);
    MPI_Error_string(*err, error_info, &error_info_length);
    MPIX_Comm_failure_ack(*comm);
    MPIX_Comm_shrink(*comm, new_comm);
    MPI_Comm_rank(*new_comm, &new_rank);
    MPI_Comm_size(*new_comm, &new_size);
    printf("NEW COMMUNICATOR: process map updated (%d -> %d) / (%d -> %d)\n", rank, size, new_rank, new_size);
    printf("Trying to revoke old communicator %d / %d\n", rank, size);
    mpi_call_status = MPIX_Comm_revoke(*comm);
    printf("Big Ol' communicator revoked %d / %d, SUCCESS: %s\n", rank, size, mpi_call_status == MPI_SUCCESS ? "yes" : "no");
    printf("Before old comm free %d / %d\n", rank, size);
    printf("After old comm free %d / %d; call status = %d; MPI_SUCCESS = %d", rank, size, mpi_call_status, MPI_SUCCESS);
    *comm = *new_comm;
    printf("Shrinked down old communicator. Replaced by new shrinked one\n");
    int test_rank, test_size;
    MPI_Comm_rank(*comm, &test_rank);
    MPI_Comm_size(*comm, &test_size);
    printf("New parameters : %d / %d\n", test_rank, test_size);
    MPIX_Comm_failure_ack(*comm);
    */
}