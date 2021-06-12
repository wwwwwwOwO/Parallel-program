/*
 * Compile:  mpicc -g -Wall -o MPI_solver MPI_solver.c -lm
 * Run:      mpiexec -n <comm_sz> ./MPI_solver
 *
 * Input:    nbody.txt
 *
 * Output:   N*dT时间后各个body的位置、速度
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define SCALE 1024
#define dT 0.005
#define G 1
#define N 20

#define DIM 3
#define X 0
#define Y 1
#define Z 2

#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}

void Get_input(double masses[], double loc_pos[], double loc_vel[], 
    int loc_n, int my_rank, int comm_sz, MPI_Comm comm);

void Check_for_error(int local_ok, char fname[], char message[],
    MPI_Comm  comm);

void Print(double masses[],double pos[], double vel[]);

int main(void) {
    int i, q, k, phase, owner;
    double x_diff, y_diff, z_diff, dist, dist_cubed;
    double start, stop;
    int my_rank, comm_sz;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);            //Initialize
    comm = MPI_COMM_WORLD;           //Communicator
    MPI_Comm_size(comm, &comm_sz);   //comm_sz stores the number of processes
    MPI_Comm_rank(comm, &my_rank);   //the ID of the running process of communicator

    int loc_n = SCALE/comm_sz;
    int source = (my_rank+1)%comm_sz;
    int dest = (my_rank-1+comm_sz)%comm_sz;

    double* masses;
    double *loc_pos, *tmp_pos;
    double *loc_forces, *tmp_forces;
    double *loc_vel;
    double *pos, *vel;

        
    double force_qk[DIM];
    
    masses = malloc(SCALE*sizeof(double));
    loc_pos = malloc(loc_n*DIM*sizeof(double));
    tmp_pos = malloc(loc_n*DIM*sizeof(double));
    loc_vel = malloc(loc_n*DIM*sizeof(double));
    loc_forces = malloc(loc_n*DIM*sizeof(double));
    tmp_forces = malloc(loc_n*DIM*sizeof(double));

    Get_input(masses, loc_pos, loc_vel, loc_n, my_rank, comm_sz, comm);

    GET_TIME(start);
    for(i=0; i<N; i++){
        // copy loc_pos into tmp_pos
        memcpy(tmp_pos, loc_pos, loc_n*DIM*sizeof(double));
        memset(loc_forces, 0, loc_n*DIM*sizeof(double));
        memset(tmp_forces, 0, loc_n*DIM*sizeof(double));
        owner = my_rank;

        // Compute forces due to interactions among local particles
        for(q=0; q<loc_n; q++){
            for(k=0; k<loc_n; k++){
                    if(loc_n*my_rank+q<loc_n*owner+k){
                        x_diff = loc_pos[q*DIM+X] - tmp_pos[k*DIM+X];
                        y_diff = loc_pos[q*DIM+Y] - tmp_pos[k*DIM+Y];
                        z_diff = loc_pos[q*DIM+Z] - tmp_pos[k*DIM+Z];
                        dist = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
                        dist_cubed = dist*dist*dist;

                        force_qk[X] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*x_diff;
                        force_qk[Y] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*y_diff;
                        force_qk[Z] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*z_diff;

                        loc_forces[q*DIM+X] -= force_qk[X];
                        loc_forces[q*DIM+Y] -= force_qk[Y];
                        loc_forces[q*DIM+Z] -= force_qk[Z];
                        tmp_forces[k*DIM+X] += force_qk[X];
                        tmp_forces[k*DIM+Y] += force_qk[Y];
                        tmp_forces[k*DIM+Z] += force_qk[Z];
                    }
            }
        }

        for(phase=1; phase<comm_sz; phase++){
            
            MPI_Send(tmp_forces, loc_n*DIM, MPI_DOUBLE, dest, 0, comm);
            MPI_Send(tmp_pos, loc_n*DIM, MPI_DOUBLE, dest, 0, comm);
            MPI_Recv(tmp_forces, loc_n*DIM, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
            MPI_Recv(tmp_pos, loc_n*DIM, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);

            owner = (my_rank+phase)%comm_sz;

            for(q=0; q<loc_n; q++){
                for(k=0; k<loc_n; k++){
                    if(loc_n*my_rank+q<loc_n*owner+k){
                        x_diff = loc_pos[q*DIM+X] - tmp_pos[k*DIM+X];
                        y_diff = loc_pos[q*DIM+Y] - tmp_pos[k*DIM+Y];
                        z_diff = loc_pos[q*DIM+Z] - tmp_pos[k*DIM+Z];
                        dist = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
                        dist_cubed = dist*dist*dist;

                        force_qk[X] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*x_diff;
                        force_qk[Y] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*y_diff;
                        force_qk[Z] = G*masses[loc_n*my_rank+q]*masses[loc_n*owner+k]/dist_cubed*z_diff;

                        loc_forces[q*DIM+X] -= force_qk[X];
                        loc_forces[q*DIM+Y] -= force_qk[Y];
                        loc_forces[q*DIM+Z] -= force_qk[Z];
                        tmp_forces[k*DIM+X] += force_qk[X];
                        tmp_forces[k*DIM+Y] += force_qk[Y];
                        tmp_forces[k*DIM+Z] += force_qk[Z];
                    }
                }
            } 
        }
        MPI_Send(tmp_forces, loc_n*DIM, MPI_DOUBLE, dest, 0, comm);
        MPI_Send(tmp_pos, loc_n*DIM, MPI_DOUBLE, dest, 0, comm);
        MPI_Recv(tmp_forces, loc_n*DIM, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
        MPI_Recv(tmp_pos, loc_n*DIM, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);


        for(q=0; q<loc_n; q++){
            loc_forces[q*DIM+X] += tmp_forces[q*DIM+X];
            loc_forces[q*DIM+Y] += tmp_forces[q*DIM+Y];
            loc_forces[q*DIM+Z] += tmp_forces[q*DIM+Z];
        }

        for(q=0; q<loc_n; q++){
            loc_vel[q*DIM+X] += dT/masses[loc_n*my_rank+q]*loc_forces[q*DIM+X];
            loc_vel[q*DIM+Y] += dT/masses[loc_n*my_rank+q]*loc_forces[q*DIM+Y];
            loc_vel[q*DIM+Z] += dT/masses[loc_n*my_rank+q]*loc_forces[q*DIM+Z];

            loc_pos[q*DIM+X] += loc_vel[q*DIM+X]*dT;
            loc_pos[q*DIM+Y] += loc_vel[q*DIM+Y]*dT;
            loc_pos[q*DIM+Z] += loc_vel[q*DIM+Z]*dT;
        }
    }
    GET_TIME(stop);

    if (my_rank  == 0){
        pos = malloc(SCALE*DIM*sizeof(double));
        vel = malloc(SCALE*DIM*sizeof(double));
    }
    MPI_Gather(loc_pos, loc_n*DIM, MPI_DOUBLE, pos, loc_n*DIM, MPI_DOUBLE, 0, comm);
    MPI_Gather(loc_vel, loc_n*DIM, MPI_DOUBLE, vel,loc_n*DIM, MPI_DOUBLE, 0, comm);
    if(my_rank == 0){
        Print(masses, pos, vel);
        printf("MPI version run time: %e\n", stop-start);   
        free(pos);
        free(vel);
    }

    free(masses);
    free(loc_pos);
    free(tmp_pos);
    free(loc_vel);
    free(loc_forces);
    free(tmp_forces);

    MPI_Finalize();      //Finalize
    return 0;

} /* main */
void Get_input(double masses[], double loc_pos[], double loc_vel[], 
    int loc_n, int my_rank, int comm_sz, MPI_Comm comm){
    int i;
    int local_ok = 1;
    double* pos = NULL;
    double* vel = NULL;

    if(my_rank==0){
        FILE* data= fopen("nbody.txt", "rb+");
        if(data==NULL)
            fprintf(stderr, "Can not open data file.\n");

        pos = malloc(SCALE*DIM*sizeof(double));
        vel = malloc(SCALE*DIM*sizeof(double));
        for(i=0; i<SCALE; i++){
            fscanf(data, "%lf %lf %lf %lf %lf %lf %lf", &masses[i], &pos[i*DIM+X], &pos[i*DIM+Y], &pos[i*DIM+Z], &vel[i*DIM+X], &vel[i*DIM+Y], &vel[i*DIM+Z]);
        }
        fclose(data);

        MPI_Scatter(pos, loc_n*DIM, MPI_DOUBLE, 
        loc_pos, loc_n*DIM, MPI_DOUBLE, 0, comm);

        MPI_Scatter(vel, loc_n*DIM, MPI_DOUBLE, 
        loc_vel, loc_n*DIM, MPI_DOUBLE, 0, comm);

        free(pos);
        free(vel);
    }
    else{
        MPI_Scatter(pos, loc_n*DIM, MPI_DOUBLE, 
        loc_pos, loc_n*DIM, MPI_DOUBLE, 0, comm);

        MPI_Scatter(vel, loc_n*DIM, MPI_DOUBLE, 
        loc_vel, loc_n*DIM, MPI_DOUBLE, 0, comm);
    }

    MPI_Bcast(masses, SCALE, MPI_DOUBLE, 0, comm);

    if(SCALE % comm_sz != 0) local_ok = 0;

    Check_for_error(local_ok, "Get_input",
                "data_count must be evenly divisible by comm_sz",
                comm);
}


void Check_for_error(
     int       local_ok   /* in */,
     char      fname[]    /* in */,
     char      message[]  /* in */,
     MPI_Comm  comm       /* in */) {
   int ok;

   MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
   if (ok == 0) {
      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
         fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname,
               message);
         fflush(stderr);
      }
      MPI_Finalize();
      exit(-1);
   }
}  /* Check_for_error */

void Print(double masses[], double pos[], double vel[]){
   int q;
   FILE*fp = fopen ("nbody_last.txt", "w+");
   for(q=0; q<SCALE; q++){
      fprintf(fp, "%15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n",masses[q],  pos[q][X], pos[q][Y], pos[q][Z], vel[q][X], vel[q][Y], vel[q][Z]);
   }
}  /* Print*/
