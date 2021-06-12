/*
 * Compile:  gcc -g -Wall -o Pthread_solver Pthread_solver.c -lm -lpthread
 * Run:      ./Pthread_solver <thread_count>
 *
 * Input:    nbody.txt
 *
 * Output:   N*dT时间后各个body的位置、速度
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

// #define DEBUG
#define MAX_THREADS 32
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

int thread_count;

double masses[SCALE];
double pos[SCALE][DIM];
double vel[SCALE][DIM];

double force_qk[DIM];
double forces[SCALE][DIM];
double loc_forces[MAX_THREADS][SCALE][DIM];

void Usage(char prog_name[]);
void Get_args(char* argv[], int* thread_count_p);
void* thread_initForces(void* args_p);
void* thread_compForce(void* args_p);
void* thread_reducForce(void* args_p);
void* thread_compPosAndVel(void* args_p);
void Print();

int main(int argc, char* argv[]) {
   int i;
   int q;
   long thread;
   double start, stop;
   if (argc != 2) Usage(argv[0]);
   Get_args(argv, &thread_count);
   

   pthread_t *thread_handles=malloc(thread_count*sizeof(pthread_t));
   
   FILE* data= fopen("nbody.txt", "rb+");
   if(data==NULL)
      fprintf(stderr, "Can not open data file.\n");

   for(i=0; i<SCALE; i++){
      fscanf(data, "%lf %lf %lf %lf %lf %lf %lf", &masses[i], &pos[i][X], &pos[i][Y], &pos[i][Z], &vel[i][X], &vel[i][Y], &vel[i][Z]);
   }
   fclose(data);

   GET_TIME(start);
   for(i=0; i<N; i++){
        for(thread=0; thread<thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, thread_initForces, (void*) thread);
        
        for(thread=0; thread<thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);
        
        
        for(thread=0; thread<thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, thread_compForce, (void*) thread);

        for(thread=0; thread<thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);
        
        for(thread=0; thread<thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, thread_reducForce, (void*) thread);

        for(thread=0; thread<thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);
        
        for(thread=0; thread<thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, thread_compPosAndVel, (void*) thread);

        for(thread=0; thread<thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);        
   }
   GET_TIME(stop);
   Print();

   printf("Pthread run time: %e\n", stop-start);   

   return 0;
}  /* main */
/*---------------------------------------------------------------------
 * Function:  Usage 
 * Purpose:   Print a message showing how to run the program and quit
 * In arg:    prog_name:  the name of the program from the command line
 */
void Usage(char prog_name[]) {
   fprintf(stderr, "usage: %s <thread_count> (thread_count <= 32)\n", prog_name);
   exit(0);
}  /* Usage */

/*---------------------------------------------------------------------
 * Function:  Get_args
 * Purpose:   Get the command line arguments
 * In arg:    argv:  strings from command line
 * Out args:  thread_count_p: number of threads
 *            n_p: number of elements
 */
void Get_args(char* argv[], int* thread_count_p) {
   *thread_count_p = strtol(argv[1], NULL, 10);
   if(*thread_count_p>32)
    Usage(argv[0]);
}  /* Get_args */

void* thread_initForces(void* args_p){
    long my_rank = (long)args_p;
    int q;
    for(q=0; q<SCALE; q++)
        loc_forces[my_rank][q][X] = loc_forces[my_rank][q][Y] = loc_forces[my_rank][q][Z] = 0;
}

void* thread_compForce(void* args_p){
    int q, k;
    long my_rank = (long)args_p;
    int my_n = SCALE/thread_count;
    int my_first_q = my_rank*my_n, my_last_q = my_first_q+my_n;
    double x_diff, y_diff, z_diff, dist, dist_cubed;
    double force_qk[DIM];    

    for(q=my_first_q; q<my_last_q; q++){
        for(k=q+1; k<SCALE; k++){
            x_diff = pos[q][X] - pos[k][X];
            y_diff = pos[q][Y] - pos[k][Y];
            z_diff = pos[q][Z] - pos[k][Z];
            dist = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            dist_cubed = dist*dist*dist;

            force_qk[X] = G*masses[q]*masses[k]/dist_cubed*x_diff;
            force_qk[Y] = G*masses[q]*masses[k]/dist_cubed*y_diff;
            force_qk[Z] = G*masses[q]*masses[k]/dist_cubed*z_diff;

            loc_forces[my_rank][q][X] -= force_qk[X];
            loc_forces[my_rank][q][Y] -= force_qk[Y];
            loc_forces[my_rank][q][Z] -= force_qk[Z];
            loc_forces[my_rank][k][X] += force_qk[X];
            loc_forces[my_rank][k][Y] += force_qk[Y];
            loc_forces[my_rank][k][Z] += force_qk[Z];
        }
    }
}
void* thread_reducForce(void* args_p){
    long my_rank = (long)args_p;
    int my_n = SCALE/thread_count;
    int my_first_q = my_rank*my_n, my_last_q = my_first_q+my_n;
    int q, thread;
    for(q=my_first_q; q<my_last_q; q++){
        forces[q][X] = forces[q][Y] = forces[q][Z] = 0;
        for(thread=0; thread<thread_count; thread++){
            forces[q][X] += loc_forces[thread][q][X];
            forces[q][Y] += loc_forces[thread][q][Y];
            forces[q][Z] += loc_forces[thread][q][Z];
        }
    }
}
void* thread_compPosAndVel(void* args_p){
    long my_rank = (long)args_p;
    int my_n = SCALE/thread_count;
    int my_first_q = my_rank*my_n, my_last_q = my_first_q+my_n;
    int q;
    for(q=my_first_q; q<my_last_q; q++){
        vel[q][X] += dT/masses[q]*forces[q][X];
        vel[q][Y] += dT/masses[q]*forces[q][Y];
        vel[q][Z] += dT/masses[q]*forces[q][Z];
        pos[q][X] += vel[q][X]*dT;
        pos[q][Y] += vel[q][Y]*dT;
        pos[q][Z] += vel[q][Z]*dT;
    }
}


void Print(){
   int q;
   FILE*fp = fopen ("nbody_last.txt", "w+");
   for(q=0; q<SCALE; q++){
      fprintf(fp, "%15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n",masses[q], pos[q][X], pos[q][Y], pos[q][Z], vel[q][X], vel[q][Y], vel[q][Z]);
   }
}  /* Print*/

