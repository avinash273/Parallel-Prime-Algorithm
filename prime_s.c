/*
Serialized Prime Sieve of Eratosthenes Algorithm

Done By:
Avinash Shanker
ID: 1001668570

Implemented on Stampede2 (Texas Advanced Computing Center)

•	To compile: mpicc prime_s.c
•	Run: ibrun -np 1 ./a.out 1000000
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>

int main(int argc, char** argv){
    unsigned long long int i,j;
    struct timeval  tv1, tv2;
    int *primes;
    int z = 1;
    int limit = 0;
    double time_spent = 0.0;
    gettimeofday(&tv1, NULL);
    limit = atoi(argv[1]);
    primes = malloc(sizeof(int) * limit);

    for (i = 2;i < limit; i++)
        primes[i] = 1;

    for (i = 2;i < limit; i++)
        if (primes[i])
            for (j = i;i * j < limit; j++)
                primes[i * j] = 0;

    printf("\nPrime numbers are:\n");
    for (i = 2;i < limit; i++)
        if (primes[i])
            printf("%d\n", i);

    gettimeofday(&tv2, NULL);
    printf ("Number N:%d,Time:%f\n",limit,
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));

return 0;
}