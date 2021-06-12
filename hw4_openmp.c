/* File:      hw4_openmp.c
 * Purpose:   Implementation of parallel count sort
 *
 * Compile:   gcc -fopenmp -o hw4_openmp hw4_openmp.c
 *
 * Run:      ./hw4_openmp <thread_count> <n>
 *
 * Input:   none
 * Output:  runtime of serial count sort
 *          runtime of parallel count sort
 *          runtime of qsort library function
 *
 */

 
void Usage(char prog_name[]);
void Get_args(char* argv[], int* thread_count_p, int* n_p);
void Gen_data(int a[], int n);
void Count_sort_serial(int a[], int n);
void Count_sort_parallel(int a[], int n, int thread_count);
void Library_qsort(int a[], int n);
int  My_compare(const void* a, const void* b);
void Print_data(int a[], int n, char msg[]);
int  Check_sort(int a[], int n);

int main(int argc, char* argv[]) {
   int n, thread_count;
   int *a, *copy;
   double start, stop;
   
   /* please choose terms 'n', and the threads 'thread_count' here. */
   n = 10;
   thread_count = 4;

   /* You can also get number of threads from command line */
   if (argc != 3) Usage(argv[0]);
   Get_args(argv, &thread_count, &n);
   
   /* Allocate storage and generate data for a */
   a = malloc(n*sizeof(int));
   Gen_data(a, n);
   
   /* Allocate storage for copy */
   copy = malloc(n*sizeof(int));
   
   /* Serial count sort */
   memcpy(copy, a, n*sizeof(int));
#  ifdef DEBUG   
   Print_data(copy, n, "Original: Serial sort a");
#  endif
   GET_TIME(start);
   Count_sort_serial(copy, n);
   GET_TIME(stop);
#  ifdef DEBUG   
   Print_data(copy, n, "Sorted: Serial sort a");
#  endif
   if (!Check_sort(copy, n))
      printf("Serial sort failed\n");
   printf("Serial run time: %e\n\n", stop-start);

   /* Parallel count sort */
   memcpy(copy, a, n*sizeof(int));
#  ifdef DEBUG   
   Print_data(copy, n, "Original: Parallel qsort a");
#  endif
   GET_TIME(start);
   Count_sort_parallel(copy, n, thread_count);
   GET_TIME(stop);
#  ifdef DEBUG   
   Print_data(copy, n, "Sorted: Parallel sort a");
#  endif
   if (!Check_sort(copy, n))
      printf("Parallel sort failed\n");
   printf("Parallel run time: %e\n\n", stop-start);   
   
   /* qsort library */
   memcpy(copy, a, n*sizeof(int));
#  ifdef DEBUG   
   Print_data(copy, n, "Original: Library qsort a");
#  endif   
   GET_TIME(start);
   Library_qsort(copy, n);
   GET_TIME(stop);
#  ifdef DEBUG   
   Print_data(copy, n, "Sorted: Library qsort a");
#  endif
   if (!Check_sort(copy, n))
      printf("Library sort failed\n");
   printf("qsort run time: %e\n", stop-start);   

   free(a);
   free(copy);
   
   return 0;
}  /* main */

/*---------------------------------------------------------------------
 * Function:  Usage 
 * Purpose:   Print a message showing how to run the program and quit
 * In arg:    prog_name:  the name of the program from the command line
 */
void Usage(char prog_name[]) {
   fprintf(stderr, "usage: %s <thread_count> <n>\n", prog_name);
   exit(0);
}  /* Usage */


/*---------------------------------------------------------------------
 * Function:  Get_args
 * Purpose:   Get the command line arguments
 * In arg:    argv:  strings from command line
 * Out args:  thread_count_p: number of threads
 *            n_p: number of elements
 */
void Get_args(char* argv[], int* thread_count_p, int* n_p) {
   *thread_count_p = strtol(argv[1], NULL, 10);
   *n_p = strtol(argv[2], NULL, 10);
}  /* Get_args */


/*---------------------------------------------------------------------
 * Function:  Gen_data
 * Purpose:   Generate random ints in the range 1 to n
 * In args:   n: number of elements
 * Out arg:   a: array of elements
 */

void Gen_data(int a[], int n) {
   int i;
   
   for (i = 0; i < n; i++)
      a[i] = random() % n + 1; // (double) RAND_MAX;

#  ifdef DEBUG
   Print_data(a, n, "a");
#  endif
}  /* Gen_data */


/*---------------------------------------------------------------------
 * Function:     Count_sort_serial
 * Purpose:      sort elements in an array using count sort
 * In args:      n: number of elements
 * In/out arg:   a: array of elements
 */

void Count_sort_serial(int a[], int n) {
   int i, j, count;
   int* temp = malloc(n*sizeof(int));
   
   for (i = 0; i < n; i++) {
      count = 0;
      for (j = 0; j < n; j++) 
         if (a[j] < a[i])
            count++;
         else if (a[j] == a[i] && j < i)
            count++;
      temp[count] = a[i];
   }
   
   memcpy(a, temp, n*sizeof(int));
   free(temp);
}  /* Count_sort_serial */


void Count_sort_parallel(int a[], int n, int thread_count) {
   int *temp = malloc(n * sizeof(int));

   # pragma omp parallel num_threads(thread_count)
   {
      int my_n = n / thread_count;
      int my_rank = omp_get_thread_num();
      int my_first_i = my_n * my_rank;
      int my_last_i = my_first_i + my_n;
      int my_count;

      for(int i = my_first_i; i < my_last_i; i++){
         my_count = 0;
         for (int j = 0; j < n; j++){
            if (a[j] < a[i] || (a[j] == a[i] && j < i)){
               my_count++;
            }
         }
         temp[my_count] = a[i];
      }
   }

   memcpy(a, temp, n * sizeof(int));
   free(temp);
}
/*---------------------------------------------------------------------
 * Function:     Library_qsort
 * Purpose:      sort elements in an array using qsort library function
 * In args:      n: number of elements
 * In/out arg:   a: array of elements
 */

void Library_qsort(int a[], int n) {
   qsort(a, n, sizeof(int), My_compare);
}  /* Library_qsort */

/*---------------------------------------------------------------------
 * Function:     My_compare
 * Purpose:      compare integer elements for use with qsort function
 * In args:      element a, element b
 * Return val:   positive if a > b, negative if b > a, 0 if equal
 */
int My_compare(const void* a, const void* b) {
   const int* int_a = (const int*) a;
   const int* int_b = (const int*) b;
   
   return (*int_a - *int_b);
}  /* My_compare */


/*---------------------------------------------------------------------
 * Function:  Print_data
 * Purpose:   print an array
 * In args:   a: array of elements
 *            n: number of elements
 *            msg: name of array
 */

void Print_data(int a[], int n, char msg[]) {
   int i;

   printf("%s = ", msg);
   for (i = 0; i < n; i++)
      printf("%d ", a[i]);
   printf("\n");
}  /* Print_data */


/*---------------------------------------------------------------------
 * Function:  Check_sort
 * Purpose:   Determine whether an array is sorted
 * In args:   a: array of elements
 *            n: number of elements
 * Ret val:   true if sorted, false if not sorted
 */

int  Check_sort(int a[], int n) {
   int i;

   for (i = 1; i < n; i++)
      if (a[i-1] > a[i]) return 0;
   return 1;
}  /* Check_sort */
