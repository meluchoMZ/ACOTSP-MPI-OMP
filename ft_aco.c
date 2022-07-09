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
MPI_Errhandler ft_respawn_on_failure_error_handler = NULL;
int spawn_threshold = 0;

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
    MPI_Comm_create_errhandler(&FT_respawn_on_failure, &ft_respawn_on_failure_error_handler);
}

void FT_finalize(void)
{
    MPI_Errhandler_free(&ft_abort_on_failure_error_handler);
    MPI_Errhandler_free(&ft_ignore_on_failure_error_handler);
    MPI_Errhandler_free(&ft_respawn_on_failure_error_handler);
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
        printf("[%s] % d / %d Repair: can shrink the communicator\n", get_current_time(time), *global_rank, *global_procs);
        MPI_Comm_size(*scomm, &ns);
        printf("[%s] % d / %d Repair: can do first MPI_Comm_size\n", get_current_time(time), *global_rank, *global_procs);
        MPI_Comm_size(*comm, &nc);
        printf("[%s] %d / %d Repair: found %d dead processes\n", get_current_time(time), *global_rank, *global_procs, nc - ns);
        //MPIX_Comm_revoke(*comm);
        //printf("[%s] Repair: communicator revoked\n", get_current_time(time));
        if (MPI_COMM_WORLD != *comm) {
            printf("[%s] %d / %d Repair: communicator is Not THE WORLD\n", get_current_time(time), *global_rank, *global_procs);
            int rc = MPI_Comm_free(comm);
            if (rc == MPI_SUCCESS) {
                printf("[%s] %d / %d Repair: communicator destroyed\n", get_current_time(time), *global_rank, *global_procs);
            } else {
                char error_info[MPI_MAX_ERROR_STRING];
                int error_info_length;
                MPI_Error_string(rc, error_info, &error_info_length);
                printf("[%s] %d / %d Repair: communicator not destroyed: %s", get_current_time(time), *global_rank, *global_procs, error_info);
            }
        }
        if (*scomm == MPI_COMM_NULL) {
            printf("[%s] %d / %d Repair: communicator is NULL\n", get_current_time(time), *global_rank, *global_procs);
        }
        *comm = *scomm;
        printf("[%s] %d / %d Repair: communicator reasigned\n", get_current_time(time), *global_rank, *global_procs);
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

    printf("[%s] %d / %d found error: %s\n", get_current_time(time), *global_rank, *global_procs, error_info);
    printf("[%s] %d / %d Error handler: entering\n", get_current_time(time), *global_rank, *global_procs);
    MPI_Comm new_comm;
    printf("[%s] %d / %d Error handler: before replace\n", get_current_time(time), *global_rank, *global_procs);
    MPIX_Comm_replace(*comm, &new_comm);
    //MPI_Comm_free(comm);
    //*comm = new_comm;
    printf("[%s] %d / %d Error handler: exiting handler\n", get_current_time(time), *global_rank, *global_procs);
    if (*comm == MPI_COMM_NULL) {
        printf("[%s] %d / %d Error handler: communicator is NULL\n", get_current_time(time), *global_rank, *global_procs);
    }
    MPI_Comm_group(*comm, &group_c);
    //MPI_Group_translate_ranks(group_f, number_of_dead, &ranks_gf, group_c, &ranks_gc);
    int shrinked_rank, shrinked_size;
	MPI_Comm_rank(*comm, &shrinked_rank);
    printf("[handler - %s] %d / %d: new rank = %d\n", get_current_time(time), rank, size, shrinked_rank);
    MPI_Comm_size(*comm, &shrinked_size);
    printf("[handler - %s] %d / %d: new size = %d\n", get_current_time(time), rank, size, shrinked_size);
	printf("[handler - %s] %d / %d re-ranking processes. New map %d / %d -> %d / %d\n", get_current_time(time), *global_rank, *global_procs, rank, size, shrinked_rank, shrinked_size);
	*global_rank = (int) shrinked_rank;
	*global_procs = (int) shrinked_size;
    printf("[handler - %s] %d / %d: using new global rank and size values\n", get_current_time(time), *global_rank, *global_procs);
}


char * get_current_time(char * output) {
    time_t t;
    struct tm * time_s;
    t = time(NULL);
    time_s = localtime(&t);
    sprintf(output, "%02d:%02d:%03d", time_s->tm_hour, time_s->tm_min, time_s->tm_sec);
    return output;
}

int attach_to_comm(MPI_Comm * comm) {
    printf("[attach_to_comm] entering\n");
    respawn(MPI_COMM_NULL, comm);
    printf("[attach_to_comm] leaving\n");
}

void respawn(MPI_Comm comm, MPI_Comm * newcomm)
{
    // respawns 1 process
    MPI_Comm icomm, /* the intercomm between the spawnees and the old (shrinked) world */
             scomm, /* the local comm for each sides of icomm */
             mcomm; /* the intracomm, merged from icomm */
    int verbose = 1;
    MPI_Group cgrp, sgrp, dgrp;
    int rc, flag, rflag, i, nc, ns, nd, crank, srank, drank;
    int isNewSpawnee = 1;

    printf("[respawn] entering respawn function\n");
redo:
    printf("[respawn] redo:\n");
    if( comm == MPI_COMM_NULL ) { /* am I a new process? */
        /* I am a new spawnee, waiting for my new rank assignment
         * it will be sent by rank 0 in the old world */
        printf("[respawn new spawnee] I am a new spawnee -> waiting for my new rank assignment\n");
        MPI_Comm_get_parent(&icomm);
        printf("[respanw new spawnee] new spawnee: parent obtained\n");
        scomm = MPI_COMM_WORLD;
        printf("[respawn new spawnee] new spawnee: scomm = WORLD\n");
        isNewSpawnee = 0;
    }
    else {
        /* I am a survivor: Spawn the appropriate number
         * of replacement processes (we check that this operation worked
         * before we procees further) */
        /* First: remove dead processes */
        printf("[respawn] I am a survivor :)\n");
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
            if( verbose ) fprintf(stderr, "%04d: comm_spawn failed, redo\n", *global_rank);
            goto redo;
        }
        printf("[respawn] updating sizes and ranks\n");
        MPI_Comm_rank(comm, global_rank);
        MPI_Comm_size(comm, global_procs);
        printf("[respawn] leaving respawn function\n");
    }

#if 0
    /* Move this failure around to see what happens */
    if( 0 == rank ) {
        fprintf(stderr, "%04d: injecting another failure!\n", rank);
        raise(SIGKILL);
    }
#endif

    /* Merge the intercomm, to reconstruct an intracomm (we check
     * that this operation worked before we proceed further) */
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] before intercomm merge\n");
    rc = MPI_Intercomm_merge(icomm, 1, &mcomm);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] after intercomm merge\n");
    rflag = flag = (MPI_SUCCESS==rc);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] before MPI_COMM_AGREE\n");
    MPIX_Comm_agree(scomm, &flag);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] freeing communicator\n");
    if( MPI_COMM_WORLD != scomm ) MPI_Comm_free(&scomm);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] communicator FREED\n");
    MPIX_Comm_agree(icomm, &rflag);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] after agree again\n");
    MPI_Comm_free(&icomm);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] intercomm freeing\n");
    if( !(flag && rflag) ) {
        if( MPI_SUCCESS == rc ) {
            if (isNewSpawnee == 0) {
                printf("[respawn new spawnee] mcomm freeing\n");
            }
            MPI_Comm_free(&mcomm);
            if (isNewSpawnee == 0) {
                printf("[respawn new spawnee] mcomm FREED\n");
            }
        }
        if( verbose ) fprintf(stderr, "%04d: Intercomm_merge failed, redo\n", *global_rank);
        goto redo;
    }

    if (isNewSpawnee == 0) {
        printf("[respawn new spawnee] error handler restoration 1\n");
        printf("[respawn new spawnee] error handler restoration 2\n");
        printf("[respawn new spawnee] error handler restoration 3\n");
        printf("[respawn new spawnee] error handler restoration 4\n");
        printf("[respawn new spawnee] error handler restoration 5\n");
        printf("[respawn new spawnee] error handler restoration 6\n");
    }
    /* restore the error handler */
    if( MPI_COMM_NULL != comm ) {
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] COMM IS NOT NULL!!!!!!\n");
        MPI_Errhandler errh;
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] before GETTER\n");
        MPI_Comm_get_errhandler( comm, &errh );
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] after GETTER\n");
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] before SETTER\n");
        MPI_Comm_set_errhandler( mcomm, errh );
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] after SETTER\n");
        if (isNewSpawnee == 0)
            printf("[respawn new spawnee] error handler restored\n");
    } else {
        if (isNewSpawnee == 0) {
            printf("[respawn new spawnee] COMM IS NULL!!!!!!\n");
            //MPI_Comm_set_errhandler(mcomm, FT_RESPAWN_ON_FAILURE_HANDLER);
            printf("[respawn new spawnee] error handler set to FT_RESPAWN ON FAILURE\n");
        }
    }
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] newcomm = mcomm\n");
    *newcomm = mcomm;
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] after newcomm assignment\n");

    MPI_Comm_rank(mcomm, global_rank);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] global rank updated\n");
    MPI_Comm_size(mcomm, global_procs);
    if (isNewSpawnee == 0)
        printf("[respawn new spawnee] global size updated\n");
    return MPI_SUCCESS;
}

void FT_respawn_on_failure(MPI_Comm * comm, int * err, ...)
{
    // shrinks communicator, spawns process and obtains a new communicator including the new spawnee
    // changes handler to ERR_RETURN to avoid stagnation from spawn errors
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

    printf("[respawn handler - %s] %d / %d entering FT_RESPAWN_ON_FAILURE_HANDLER\n", get_current_time(time), *global_rank, *global_procs);
    printf("[respawn handler - %s] %d / %d found error: %s\n", get_current_time(time), *global_rank, *global_procs, error_info);
    printf("[respawn handler - %s] %d / %d Error handler: entering\n", get_current_time(time), *global_rank, *global_procs);
    MPI_Comm new_comm;
    printf("[respawn handler - %s] %d / %d Error handler: before replace\n", get_current_time(time), *global_rank, *global_procs);
    //MPIX_Comm_replace(*comm, &new_comm);
    respawn(*comm, comm);
    printf("[respawn handler - %s] %d / %d Error handler: after replace\n", get_current_time(time), *global_rank, *global_procs);
    //MPI_Comm_free(comm);
    //*comm = new_comm;
    printf("[respawn handler - %s] %d / %d Error handler: exiting handler\n", get_current_time(time), *global_rank, *global_procs);
    if (*comm == MPI_COMM_NULL) {
        printf("[respawn handler - %s] %d / %d Error handler: communicator is NULL\n", get_current_time(time), *global_rank, *global_procs);
    }
    MPI_Comm_group(*comm, &group_c);
    //MPI_Group_translate_ranks(group_f, number_of_dead, &ranks_gf, group_c, &ranks_gc);
    int shrinked_rank, shrinked_size;
	MPI_Comm_rank(*comm, &shrinked_rank);
    printf("[respawn handler - %s] %d / %d: new rank = %d\n", get_current_time(time), rank, size, shrinked_rank);
    MPI_Comm_size(*comm, &shrinked_size);
    printf("[respawn handler - %s] %d / %d: new size = %d\n", get_current_time(time), rank, size, shrinked_size);
	printf("[respawn handler - %s] %d / %d re-ranking processes. New map %d / %d -> %d / %d\n", get_current_time(time), *global_rank, *global_procs, rank, size, shrinked_rank, shrinked_size);
	*global_rank = (int) shrinked_rank;
	*global_procs = (int) shrinked_size;
    printf("[respawn handler - %s] %d / %d: using new global rank and size values\n", get_current_time(time), *global_rank, *global_procs);
}

void set_spawn_threshold(int threshold) {
    spawn_threshold = threshold;
}

void FT_elastic_respawn_on_failure(MPI_Comm * comm, int * err, ...)
{
    // same as previous one, but spawns multiple processes from 
    // a lower threshold to an upper one.
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

    printf("[elastic handler - %s] %d / %d entering FT_ELASTIC_RESPAWN_ON_FAILURE_HANDLER\n", get_current_time(time), *global_rank, *global_procs);
    printf("[elastic handler - %s] %d / %d found error: %s\n", get_current_time(time), *global_rank, *global_procs, error_info);
    printf("[elastic handler - %s] %d / %d Error handler: entering\n", get_current_time(time), *global_rank, *global_procs);
    MPI_Comm new_comm;
    printf("[elastic handler - %s] %d / %d Error handler: before replace\n", get_current_time(time), *global_rank, *global_procs);
    //MPIX_Comm_replace(*comm, &new_comm);
    // TODO: test only right number of processes are spawned
    if (spawn_threshold > (size - number_of_dead)) {
        repair(comm);
    } else {
        respawn(*comm, comm);
    }
    printf("[elastic handler - %s] %d / %d Error handler: after replace\n", get_current_time(time), *global_rank, *global_procs);
    //MPI_Comm_free(comm);
    //*comm = new_comm;
    printf("[elastic handler - %s] %d / %d Error handler: exiting handler\n", get_current_time(time), *global_rank, *global_procs);
    if (*comm == MPI_COMM_NULL) {
        printf("[elastic handler %s] %d / %d Error handler: communicator is NULL\n", get_current_time(time), *global_rank, *global_procs);
    }
    MPI_Comm_group(*comm, &group_c);
    //MPI_Group_translate_ranks(group_f, number_of_dead, &ranks_gf, group_c, &ranks_gc);
    int shrinked_rank, shrinked_size;
	MPI_Comm_rank(*comm, &shrinked_rank);
    printf("[elastic handler - %s] %d / %d: new rank = %d\n", get_current_time(time), rank, size, shrinked_rank);
    MPI_Comm_size(*comm, &shrinked_size);
    printf("[elastic handler - %s] %d / %d: new size = %d\n", get_current_time(time), rank, size, shrinked_size);
	printf("[elastic handler - %s] %d / %d re-ranking processes. New map %d / %d -> %d / %d\n", get_current_time(time), *global_rank, *global_procs, rank, size, shrinked_rank, shrinked_size);
	*global_rank = (int) shrinked_rank;
	*global_procs = (int) shrinked_size;
    printf("[elastic handler - %s] %d / %d: using new global rank and size values\n", get_current_time(time), *global_rank, *global_procs);
}
