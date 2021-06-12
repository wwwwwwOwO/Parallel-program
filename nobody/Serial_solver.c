/*
 * Compile:  gcc -g -Wall -o Serial_solver Serial_solver.c -lm
 * Run:      ./Serial_solver
 *
 * Input:    nbody.txt
 *
 * Output:   N*dT时间后各个body的位置、速度
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

double masses[SCALE];
double pos[SCALE][DIM];
double vel[SCALE][DIM];

double force_qk[DIM];
double forces[SCALE][DIM];

void Print();

int main(int argc, char* argv[]) {
   int i;
   int q, k;
   double x_diff, y_diff, z_diff, dist, dist_cubed;
   double x_a, y_a, z_a;
   double start, stop;


   FILE* data= fopen("nbody.txt", "rb+");
   if(data==NULL)
      fprintf(stderr, "Can not open data file.\n");

   for(i=0; i<SCALE; i++){
      fscanf(data, "%lf %lf %lf %lf %lf %lf %lf", &masses[i], &pos[i][X], &pos[i][Y], &pos[i][Z], &vel[i][X], &vel[i][Y], &vel[i][Z]);
   }
   fclose(data);

   GET_TIME(start);
   for(i=0; i<N; i++){
      /* 将每个body受到的力初始化为0 */
      for(q=0; q<SCALE; q++){
         forces[q][X] = forces[q][Y] = forces[q][Z] = 0;
      }

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

            forces[q][X] -= force_qk[X];
            forces[q][Y] -= force_qk[Y];
            forces[q][Z] -= force_qk[Z];
            forces[k][X] += force_qk[X];
            forces[k][Y] += force_qk[Y];
            forces[k][Z] += force_qk[Z];
         }
      }
   
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
   printf("Serial run time: %e\n", stop-start);   

   return 0;
}  /* main */

void Print(){
   int q;
   FILE*fp = fopen ("nbody_last.txt", "w+");
   for(q=0; q<SCALE; q++){
      fprintf(fp, "%15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n",masses[q], pos[q][X], pos[q][Y], pos[q][Z], vel[q][X], vel[q][Y], vel[q][Z]);
   }
}  /* Print*/

