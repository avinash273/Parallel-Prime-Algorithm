/*
Parallel Prime Sieve of Eratosthenes Algorithm

Done By:
Avinash Shanker
ID: 1001668570

Implemented on Stampede2 (Texas Advanced Computing Center)


•	To compile: mpicc prime_p.c
•	Run the code in 2 modes:
    Mode1(Prime count and Time): ibrun -np 2 ./a.out 1000000
    Mode2(Prime count, Prime numbers and Time): ibrun -np 2 ./a.out 100000 PRIMES
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define first_prime 2
#define MASTER_START_INT 1
int debug = 0;

void distributeWork(int*, int*, unsigned long*,
                    unsigned int*, unsigned int*, unsigned int*,
                    unsigned int*, unsigned int*, unsigned int*);

int main(int argc, char** argv) {

   unsigned long n = 4294967295;
   double timer;
   int p_count;
   int p_id;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;

   unsigned int n_master;
   unsigned int n_worker;
   unsigned int remainder;
   unsigned int p_work;
   unsigned int p_first;
   unsigned int p_last;
   unsigned int prime;
   int j;
   int final_count;
   char *mark_table;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    timer = -MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    if(argc >= 2)
       n = strtoul(argv[1], NULL, 0);
    if(argc == 3 && !strcmp(argv[2],"PRIMES"))
    {
       debug = 1;
    }

    distributeWork(&p_id, &p_count,
                   &n,
                   &p_last, &p_first, &p_work,
                   &n_master, &n_worker, &remainder);

    if(debug)
    {
       printf("Proc %d: Assigned %d (%d - %d)\n", p_id, p_work, p_first, p_last);
       MPI_Barrier(MPI_COMM_WORLD);
    }


    mark_table = (char*) malloc( sizeof(char)*p_work );
    if(mark_table == NULL)
    {
       printf("Failed to allocate memory for process %d\n", p_id);
       MPI_Finalize();
       exit(1);
    }

    for(j = 0; j < p_work; j++)
       mark_table[j] = '0';

    if(!p_id) mark_table[0] = '1';

    int marker;
    prime = first_prime;
    do
    {

       if(prime < p_first)
       {
          int mod = p_first % prime;
          if(mod)
             marker = prime - mod;

          else
             marker = 0;
       }
       else
       {
          marker = 2*prime - p_first;
       }


       for(j = marker; j < p_work; j += prime)
       {
          mark_table[j] = '1';
       }


       if(!p_id)
       {
          int next_index = prime - MASTER_START_INT;
          do{
             if(++next_index >= n_master)
             {
                next_index == n_master - MASTER_START_INT;
                break;
             }
          }while(mark_table[next_index] != '0');
          prime = next_index + MASTER_START_INT;
       }


       MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);


    }while(prime <= n_master );

    if(debug)
    {

       if(p_id)
       {
          MPI_Recv(&prime, 1, MPI_INT, (p_id-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       }


       printf("Proc %d: ", p_id);
       for(j = 0; j < p_work; j++)
       {
          if(mark_table[j] == '0')
             printf("%d, ", j + p_first);
       }
       printf("\n");
       fflush(stdout);


       if(p_id != (p_count - 1))
       {

          MPI_Send(&prime, 1, MPI_INT, (p_id+1), 1, MPI_COMM_WORLD);

       }
    }


    prime = 0;

    #pragma omp parallel for reduction(+:prime)
    for(j = 0; j < p_work; j++)
    {
       if(mark_table[j] == '0')
          prime++;
    }
    MPI_Reduce(&prime, &final_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    timer += MPI_Wtime();


    if(!p_id)
    {
     if(debug) printf("Count of primes up to %d: %d\n", n, final_count);
     printf("N:%lu,P:%d,Time:%10.6f\n",n,p_count,timer);
    }


    MPI_Finalize();
}



void distributeWork(int *p_id,
                    int *p_count,
                    unsigned long *n,
                    unsigned int *p_last,
                    unsigned int *p_first,
                    unsigned int *p_work,
                    unsigned int *n_master,
                    unsigned int *n_worker,
                    unsigned int *remainder)
{
    *n_master = (unsigned int) ceil(sqrt((double) *n));
    *n_worker = (*n - *n_master)/(*p_count - 1);
    *remainder = (*n - *n_master)%(*p_count - 1);


    if(*n_worker > *n_master)
    {
       *n_master = *n_worker = *n / *p_count;
       *remainder = *n % *p_count;
    }

    if(!*p_id)
    {
       *p_last = *p_work = *n_master;
       *p_first = MASTER_START_INT;
    }
    else
    {
       *p_work = *n_worker;

       if(*p_id <= *remainder)
          (*p_work)++;

       *p_last = *n_master + (*n_worker * *p_id);

       if(*p_id <= *remainder)
          *p_last += *p_id;
       else
          *p_last += *remainder;
       *p_first = *p_last - *p_work + 1;
    }
}