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

#define MASTER_START_INT 1
#define start_prime 2

int print_prime = 0;

void Distribute_Work(int*, int*, unsigned long*,
                    unsigned int*, unsigned int*, unsigned int*,
                    unsigned int*, unsigned int*, unsigned int*);

int main(int argc, char** argv) {

   double program_timer;
   int processor_count;
   int processor_ID;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   int name_len;
   unsigned long n = 4294967295;

   unsigned int master_node;
   unsigned int slave_node;
   unsigned int remainder;
   unsigned int processor_work;
   unsigned int processor_first;
   unsigned int processor_last;
   unsigned int prime;
   int j;
   int final_count;
   char *mark_index;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    program_timer = -MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &processor_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &processor_ID);

    if(argc >= 2)
       n = strtoul(argv[1], NULL, 0);
    if(argc == 3 && !strcmp(argv[2],"PRIMES"))
    {
       print_prime = 1;
    }

    Distribute_Work(&processor_ID, &processor_count,
                   &n,
                   &processor_last, &processor_first, &processor_work,
                   &master_node, &slave_node, &remainder);

    if(print_prime)
    {
       printf("Proc %d: Assigned %d (%d - %d)\n", processor_ID, processor_work, processor_first, processor_last);
       MPI_Barrier(MPI_COMM_WORLD);
    }


    mark_index = (char*) malloc( sizeof(char)*processor_work );
    if(mark_index == NULL)
    {
       printf("Failed to allocate memory for process %d\n", processor_ID);
       MPI_Finalize();
       exit(1);
    }

    for(j = 0; j < processor_work; j++)
       mark_index[j] = '0';

    if(!processor_ID) mark_index[0] = '1';

    int marker;
    prime = start_prime;
    do
    {

       if(prime < processor_first)
       {
          int check_reminder = processor_first % prime;
          if(check_reminder)
             marker = prime - check_reminder;

          else
             marker = 0;
       }
       else
       {
          marker = 2*prime - processor_first;
       }


       for(j = marker; j < processor_work; j += prime)
       {
          mark_index[j] = '1';
       }


       if(!processor_ID)
       {
          int next_index = prime - MASTER_START_INT;
          do{
             if(++next_index >= master_node)
             {
                next_index == master_node - MASTER_START_INT;
                break;
             }
          }while(mark_index[next_index] != '0');
          prime = next_index + MASTER_START_INT;
       }


       MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);


    }while(prime <= master_node );

    if(print_prime)
    {

       if(processor_ID)
       {
          MPI_Recv(&prime, 1, MPI_INT, (processor_ID-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       }


       printf("Processor %d: ", processor_ID);
       for(j = 0; j < processor_work; j++)
       {
          if(mark_index[j] == '0')
             printf("%d, ", j + processor_first);
       }
       printf("\n");
       fflush(stdout);


       if(processor_ID != (processor_count - 1))
       {

          MPI_Send(&prime, 1, MPI_INT, (processor_ID+1), 1, MPI_COMM_WORLD);

       }
    }


    prime = 0;

    #pragma omp parallel for reduction(+:prime)
    for(j = 0; j < processor_work; j++)
    {
       if(mark_index[j] == '0')
          prime++;
    }
    MPI_Reduce(&prime, &final_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    program_timer += MPI_Wtime();


    if(!processor_ID)
    {
     if(print_prime) printf("Get Primes till %d: %d\n", n, final_count);
     printf("Number:%lu Processors:%d  Time consumed:%10.6f\n",n,processor_count,program_timer);
    }


    MPI_Finalize();
}



void Distribute_Work(int *processor_ID,
                    int *processor_count,
                    unsigned long *n,
                    unsigned int *processor_last,
                    unsigned int *processor_first,
                    unsigned int *processor_work,
                    unsigned int *master_node,
                    unsigned int *slave_node,
                    unsigned int *remainder)
{
    *master_node = (unsigned int) ceil(sqrt((double) *n));
    *slave_node = (*n - *master_node)/(*processor_count - 1);
    *remainder = (*n - *master_node)%(*processor_count - 1);


    if(*slave_node > *master_node)
    {
       *master_node = *slave_node = *n / *processor_count;
       *remainder = *n % *processor_count;
    }

    if(!*processor_ID)
    {
       *processor_last = *processor_work = *master_node;
       *processor_first = MASTER_START_INT;
    }
    else
    {
       *processor_work = *slave_node;

       if(*processor_ID <= *remainder)
          (*processor_work)++;

       *processor_last = *master_node + (*slave_node * *processor_ID);

       if(*processor_ID <= *remainder)
          *processor_last += *processor_ID;
       else
          *processor_last += *remainder;
       *processor_first = *processor_last - *processor_work + 1;
    }
}