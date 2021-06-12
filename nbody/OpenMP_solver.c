/*
 * Compile:  gcc -g -Wall -o OpenMP_solver OpenMP_solver.c -lm -fopenmp
 * Run:      ./OpenMP_solver <thread_count>
 *
 * Input:    nbody.txt
 *
 * Output:   N*dT时间后各个body的位置、速度
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
// #define DEBUG

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

double masses[SCALE];
double pos[SCALE][DIM];
double vel[SCALE][DIM];

double forces[SCALE][DIM];

void Usage(char prog_name[]);
void Get_args(char* argv[], int* thread_count_p);
void Print();


int main(int argc, char* argv[]) {
   int i;
   int q;
   int thread_count, thread;
   double start, stop;
   if (argc != 2) Usage(argv[0]);
   Get_args(argv, &thread_count);
   double loc_forces[thread_count][SCALE][DIM];

   FILE* data= fopen("nbody.txt", "rb+");
   if(data==NULL)
      fprintf(stderr, "Can not open data file.\n");

   for(i=0; i<SCALE; i++){
      fscanf(data, "%lf %lf %lf %lf %lf %lf %lf", &masses[i], &pos[i][X], &pos[i][Y], &pos[i][Z], &vel[i][X], &vel[i][Y], &vel[i][Z]);
   }
   fclose(data);

   GET_TIME(start);
   for(i=0; i<N; i++){
      # pragma omp for  
         for(q=0; q<SCALE; q++){
            for(thread=0; thread<thread_count; thread++)
               loc_forces[thread][q][X] = loc_forces[thread][q][Y] = loc_forces[thread][q][Z] = 0;
         }


      # pragma omp parallel num_threads(thread_count)
      {
         int q, k;
         int my_rank = omp_get_thread_num();
         double x_diff, y_diff, z_diff, dist, dist_cubed;
         double force_qk[DIM];

         # pragma omp for
         for(q=0; q<SCALE; q++){
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
      # pragma omp for  
      for(q=0; q<SCALE; q++){
         forces[q][X] = forces[q][Y] = forces[q][Z] = 0;
         for(thread=0; thread<thread_count; thread++){
            forces[q][X] += loc_forces[thread][q][X];
            forces[q][Y] += loc_forces[thread][q][Y];
            forces[q][Z] += loc_forces[thread][q][Z];
         }
      }
       
      # pragma omp for  
      for(q=0; q<SCALE; q++){
         vel[q][X] += dT/masses[q]*forces[q][X];
         vel[q][Y] += dT/masses[q]*forces[q][Y];
         vel[q][Z] += dT/masses[q]*forces[q][Z];

         pos[q][X] += vel[q][X]*dT;
         pos[q][Y] += vel[q][Y]*dT;
         pos[q][Z] += vel[q][Z]*dT;

      }
   }
   GET_TIME(stop);
   Print();

   printf("OpenMP run time: %e\n", stop-start);   

   return 0;
}  /* main */

/*---------------------------------------------------------------------
 * Function:  Usage 
 * Purpose:   Print a message showing how to run the program and quit
 * In arg:    prog_name:  the name of the program from the command line
 */
void Usage(char prog_name[]) {
   fprintf(stderr, "usage: %s <thread_count> \n", prog_name);
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
}  /* Get_args */

void Print(){
   int q;
   FILE*fp = fopen ("nbody_last.txt", "w+");
   for(q=0; q<SCALE; q++){
      fprintf(fp, "%15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n",masses[q], pos[q][X], pos[q][Y], pos[q][Z], vel[q][X], vel[q][Y], vel[q][Z]);
   }
}  /* Print*/

