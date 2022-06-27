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

char ** gargv;
int * global_rank, *global_procs;

int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm) {
    int verbose = 1;
    MPI_Comm icomm, /* the intercomm between the spawnees and the old (shrinked) world */
             scomm, /* the local comm for each sides of icomm */
             mcomm; /* the intracomm, merged from icomm */
    MPI_Errhandler backup;
    MPI_Comm_get_errhandler(comm, &backup);
    int rc, flag, rflag, nc, ns, nd, rank;
    printf("Replace: entering replace\n");

redo:
    if( comm == MPI_COMM_NULL ) { /* am I a new process? */
        /* I am a new spawnee, waiting for my new rank assignment
         * it will be sent by rank 0 in the old world */
        MPI_Comm_get_parent(&icomm);
        scomm = MPI_COMM_WORLD;
    }
    else {
        /* I am a survivor: Spawn the appropriate number
         * of replacement processes (we check that this operation worked
         * before we procees further) */
        /* First: remove dead processes */
        MPIX_Comm_shrink(comm, &scomm);
        MPI_Comm_size(scomm, &ns);
        MPI_Comm_size(comm, &nc);
        nd = nc-ns; /* number of deads */
        if( 0 == nd ) {
            /* Nobody was dead to start with. We are done here */
            MPI_Comm_free(&scomm);
            *newcomm = comm;
            return MPI_SUCCESS;
        }
        /* We handle failures during this function ourselves... */
        MPI_Comm_set_errhandler( scomm, MPI_ERRORS_RETURN );

        rc = MPI_Comm_spawn(gargv[0], &gargv[1], nd, MPI_INFO_NULL,
                           0, scomm, &icomm, MPI_ERRCODES_IGNORE);
        flag = (MPI_SUCCESS == rc);
        MPIX_Comm_agree(scomm, &flag);
        if( !flag ) {
            if( MPI_SUCCESS == rc ) {
                MPIX_Comm_revoke(icomm);
                MPI_Comm_free(&icomm);
            }
            MPI_Comm_free(&scomm);
            if( verbose ) fprintf(stderr, "%04d: comm_spawn failed, redo\n", rank);
            goto redo;
        }
    }

#if 0
    /* Move this failure around to see what happens 
    if( 0 == rank ) {
        fprintf(stderr, "%04d: injecting another failure!\n", rank);
        raise(SIGKILL);
    }
    */
#endif
    /* Merge the intercomm, to reconstruct an intracomm (we check
     * that this operation worked before we proceed further) */
    printf("Replace: process spawned\n");
    rc = MPI_Intercomm_merge(icomm, 1, &mcomm);
    rflag = flag = (MPI_SUCCESS==rc);
    printf("Replace: before agreeing\n");
    MPIX_Comm_agree(scomm, &flag);
    printf("Replace: after agreeing\n");
    if( MPI_COMM_WORLD != scomm ) MPI_Comm_free(&scomm);
    MPIX_Comm_agree(icomm, &rflag);
    MPI_Comm_free(&icomm);
    if( !(flag && rflag) ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free(&mcomm);
        }
        if( verbose ) fprintf(stderr, "%04d: Intercomm_merge failed, redo\n", rank);
        goto redo;
    }

    /* restore the error handler */
    if( MPI_COMM_NULL != comm ) {
        MPI_Errhandler errh;
        MPI_Comm_get_errhandler( comm, &errh );
        MPI_Comm_set_errhandler( mcomm, errh );
    }
    *newcomm = mcomm;

    MPI_Comm_rank(mcomm, &rank);
    return MPI_SUCCESS;
}

void repair(MPI_Comm * comm) {
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    MPI_Comm * scomm;
    int ns, nc;
    if ( *comm != MPI_COMM_NULL) {
        scomm = (MPI_Comm *) malloc(sizeof(MPI_Comm));
        MPIX_Comm_shrink(*comm, scomm);
        printf("[%s] Repair: can shrink the communicator\n", get_current_time(time));
        MPI_Comm_size(*scomm, &ns);
        printf("[%s] Repair: can do first MPI_Comm_size\n", get_current_time(time));
        MPI_Comm_size(*comm, &nc);
        printf("[%s] Repair: found %d dead processes\n", get_current_time(time), nc - ns);
        //MPIX_Comm_revoke(*comm);
        //printf("[%s] Repair: communicator revoked\n", get_current_time(time));
        if (MPI_COMM_WORLD != *comm) {
            printf("[%s] Repair: communicator is Not THE WORLD\n", get_current_time(time));
            int rc = MPI_Comm_free(comm);
            if (rc == MPI_SUCCESS) {
                printf("[%s] Repair: communicator destroyed\n", get_current_time(time));
            } else {
                char error_info[MPI_MAX_ERROR_STRING];
                int error_info_length;
                MPI_Error_string(rc, error_info, &error_info_length);
                printf("[%s] Repair: communicator not destroyed: %s", get_current_time(time), error_info);
            }
        }
        if (*scomm == MPI_COMM_NULL) {
            printf("[%s] Repair: communicator is NULL\n", get_current_time(time));
        }
        *comm = *scomm;
        printf("[%s] Repair: communicator reasigned\n", get_current_time(time));
    }
}

void FT_set_respawn_data(char ** argv, int * mpi_id, int * NPROC) {
    gargv = argv;
    global_rank = (int *) mpi_id;
	global_procs = (int *) NPROC;
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    printf("[%s] command: %s; args: %s\n", get_current_time(time), gargv[0], gargv[1]);
}

void FT_ignore_on_failure(MPI_Comm * comm, int * err, ...)
{
    char time[12] = {'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};
    char error_info[MPI_MAX_ERROR_STRING];
    int error_info_length, size, rank, number_of_dead;
    MPI_Group group_f, group_c;
    int ranks_gf, ranks_gc;

    MPI_Error_string(*err, error_info, &error_info_length);
    MPI_Comm_size(*comm, &size);
    MPI_Comm_rank(*comm, &rank);

    MPIX_Comm_failure_ack(*comm);
    MPIX_Comm_failure_get_acked(*comm, &group_f);
    MPI_Group_size(group_f, &number_of_dead);

    printf("[%s] found error: %s\n", get_current_time(time), error_info);
    printf("[%s] Error handler: entering\n", get_current_time(time));
    MPI_Comm new_comm;
    printf("[%s] Error handler: before replace\n", get_current_time(time));
    //MPIX_Comm_replace(*comm, &new_comm);
    repair(comm);
    printf("[%s] Error handler: after replace\n", get_current_time(time));
    //MPI_Comm_free(comm);
    //*comm = new_comm;
    printf("[%s] Error handler: exiting handler\n", get_current_time(time));
    if (*comm == MPI_COMM_NULL) {
        printf("[%s] Error handler: communicator is NULL\n", get_current_time(time));
    }
    MPI_Comm_group(*comm, &group_c);
    MPI_Group_translate_ranks(group_f, number_of_dead, &ranks_gf, group_c, &ranks_gc);
	MPI_Comm_rank(*comm, &ranks_gc);
	printf("[handler - %s] re-ranking processes. New map %d / %d -> %d / %d\n", get_current_time(time), rank, size, ranks_gc, size-number_of_dead);
	*global_rank = (int) ranks_gc;
	*global_procs = (int) size - number_of_dead;
    /*
    MPI_Comm * new_comm;
    new_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    char error_info[MPI_MAX_ERROR_STRING];
    //int rank, size, error_info_length;
    int error_info_length;
    MPI_Comm * new_comm;

    if (*err == MPIX_ERR_REVOKED) {
        printf("Process %d / %d has found a REVOKED communicator\n", mpi_id, NPROC);
    } else {
        MPI_Error_string(*err, error_info, &error_info_length);
        printf("Process %d / %d has detected an MPI failure: %s\n", mpi_id, NPROC, error_info);
    }
    MPI_Group group;
    MPIX_Comm_failure_ack(*comm);
    MPIX_Comm_failure_get_acked(*comm, &group);
    MPIX_Comm_revoke(*comm);
    printf("Process %d / %d: shrinking communicator\n", mpi_id, NPROC);
    new_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    MPIX_Comm_shrink(*comm, new_comm);
    MPI_Comm_free(comm);
    comm = new_comm;
    //MPI_Comm_rank(*comm, &mpi_id);
    //MPI_Comm_size(*comm, &NPROC);
    //printf("Process %d / %d remapped to %d / %d\n", rank, size, mpi_id, NPROC);
    //int new_rank, new_size, mask;
    MPI_Group group;
    MPI_Comm * new_comm = NULL;
    MPI_Comm_rank(*comm, &rank);
    MPI_Comm_size(*comm, &size);
    MPI_Error_string(*err, error_info, &error_info_length);
    if (*err == MPI_ERR_REVOKED) {
        MPI_Irecv(new_comm, 1, MPI_PACKED, MPI_ANY_SOURCE, SHARE_COMM_TAG, *comm, NULL);
        if (new_comm == NULL) {
            // no received comm yet
            // will re enter the handler until the buffer has data
            return;
        }
    } else {
        // the first process that finds the failure
        // creates a new communicator and sends it
        printf("Process %d / %d has encountered an MPI error: %s\n", rank, size, error_info);
        new_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm));
        MPIX_Comm_shrink(*comm, new_comm);
        MPIX_Comm_agree(*new_comm, &mask)
        MPIX_Comm_revoke(*comm);
        for (i = 0; i < size; i++) {
            if (rank != i) {
                printf("Sendindg communicator to process %d / %d\n", i, size);
                MPI_Isend(new_comm, 1, MPI_PACKED, i, SHARE_COMM_TAG, *comm, NULL);
            }
        }
    }
    free(comm);
    comm = new_comm;
    MPIX_Comm_failure_ack(*comm);
    MPIX_Comm_agree(*comm, &mask);
    new_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    MPIX_Comm_failure_ack(*comm);
    MPIX_Comm_shrink(*comm, new_comm);
    MPIX_Comm_agree(*new_comm, &mask);
    MPI_Comm_free(comm);
    comm = new_comm;
    MPI_Comm_rank(*comm, &new_rank);
    MPI_Comm_size(*comm, &new_size);
    printf("Process %d / %d has encountered an MPI error: %s\n", rank, size, error_info);
    printf("Process %d / %d: using new communicator -> new process ID: %d / %d\n", rank, size, new_rank, new_size);
    
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


char * get_current_time(char * output) {
    time_t t;
    struct tm * time_s;
    t = time(NULL);
    time_s = localtime(&t);
    sprintf(output, "%02d:%02d:%03d", time_s->tm_hour, time_s->tm_min, time_s->tm_sec);
    return output;
}
