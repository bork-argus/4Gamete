#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <error.h>
#include <string.h>

typedef struct {
	char 			locusName[128];
	unsigned long 	sample;
} Locus;

/*
 * Function to compare two 24 bit numbers, finding unique pairs of
 * bits.
 */
int uniquePairs(unsigned long p1, unsigned long p2)
{
	unsigned char uniques[4] = { 0, 0, 0, 0 };

	register unsigned long bitMask = 1UL;
	for (unsigned short i = 0; i < 24; i++)
	{
		register unsigned long p1_bit = p1 & bitMask;
		register unsigned long p2_bit = p2 & bitMask;

		/* 
		 * if both p1 bit and p2 bit are set, increment unique[0] 
		 * if p1 bit is set, but p2 isn't, increment unique[1]
		 * if p2 bit is set, but p1 isn't, increment unique[2]
		 * if neither p1 bit nor p2 bit is set, increment unique[3]
		 */
		if (p1_bit && p2_bit)
			uniques[0] = 1;
		else if (p1_bit)
			uniques[1] = 1;
		else if (p2_bit)
			uniques[2] = 1;
		else
			uniques[3] = 1;

		/* shift bitMask to the left to process the next bit */
		bitMask <<= 1;
	}

	return uniques[0] + uniques[1] + uniques[2] + uniques[3];
}

/*
 *  Set this to the maximum # of threads to run to process the data.
 *  This should be set to the number of cores available to the system.
 */
#define MAX_THREADS 11

/*
 *  Holds the samples read from the input file.
 *  Dynamically allocated.
 */
static Sample *samples;

/*
 *  Holds the next currently processed row
 */
static unsigned long current_row;

/*
 *  Holds the total number of rows in the samples array.
 */
static unsigned long total_rows;

/*
 *  Tracks total matches
 */
static unsigned long total_matches = 0;

/*
 *  Controls access to current_row and total_matches
 */
static pthread_mutex_t current_row_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t total_matches_mutex = PTHREAD_MUTEX_INITIALIZER;

/*
 *  Matches array
 */
typedef struct {
	unsigned long m1;
	unsigned long m2;
} RowMatch;

typedef struct {
	RowMatch		*matches;
	unsigned long	total;
} RowMatches;


/*
 *  Threading function
 *  This function represents a chunk of work. It processes
 *  100 rows of the sample array starting at current_row.
 *  When it is done, it starts on the next 100 until there
 *  are no more to process.
 */
void *processRows(void *arg)
{
	const unsigned long process_rows = 100;
	const unsigned long match_chunk = 1000000;
	unsigned long matches = 0;

	/* Initialize the "local" row match array */
	RowMatches *rm = (RowMatches *) arg;
	rm->total = 0;
	rm->matches = NULL;

	while (1)
	{
		pthread_mutex_lock(&current_row_mutex);
		unsigned long start_row = current_row;
		current_row += process_rows;
		pthread_mutex_unlock(&current_row_mutex);

		/* if the thread detects we are done ... exit */
		if (start_row >= total_rows - 1)
			break;

		unsigned long end_row = start_row + process_rows;
		if (end_row > total_rows - 1)
			end_row = total_rows - 1;

		fprintf(stderr, ".");
		for (unsigned long row = start_row; row < end_row; row++)
		{
			for (unsigned long cmp_row = row + 1; 
								 cmp_row < total_rows; cmp_row++)
			{
				if (uniquePairs(samples[row].sample, 
								samples[cmp_row].sample) == 4)
				{
					matches++;

					/*
					 *  Keep a record of the matching rows in the structure
					 *  passed into the thread.
					 *  Because each thread has its own matching array
					 *  structure, we don't need to mutex this.
					 */
					if (rm->matches == NULL)
						rm->matches = malloc(sizeof(RowMatch) * match_chunk);
					else if (rm->total % match_chunk == 0)
						rm->matches = realloc((void *) rm->matches, 
								(rm->total + match_chunk) * sizeof(RowMatch));

					rm->matches[rm->total].m1 = row;
					rm->matches[rm->total].m2 = cmp_row;
					rm->total++;
				}
			}
		}
	}

	pthread_mutex_lock(&total_matches_mutex);
	total_matches += matches;
	pthread_mutex_unlock(&total_matches_mutex);

	return NULL;
}


int main(int argc, char* argv[])
{
    char const* const fileName = argv[1]; /* should check that argc > 1 */
	FILE* file = fopen(fileName, "r"); /* should check the result */
	char line[256];
	char filename[128];
	char numbers[24];

	samples = calloc(200000, sizeof(Sample));

	total_rows = 0;
	while (fgets(line, sizeof(line), file)) 
	{
		/* note that fgets don't strip the terminating \n, checking its
		   presence would allow to handle lines longer that sizeof(line) */

		sscanf(line, 
			"%s " 
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c "
			"%c %c %c %c", 
			samples[total_rows].filename, 
			&numbers[0],
			&numbers[1],
			&numbers[2],
			&numbers[3],
			&numbers[4],
			&numbers[5],
			&numbers[6],
			&numbers[7],
			&numbers[8],
			&numbers[9],
			&numbers[10],
			&numbers[11],
			&numbers[12],
			&numbers[13],
			&numbers[14],
			&numbers[15],
			&numbers[16],
			&numbers[17],
			&numbers[18],
			&numbers[19],
			&numbers[20],
			&numbers[21],
			&numbers[22],
			&numbers[23]);

		/* convert numbers into an integer */
		unsigned long n = 0;

		samples[total_rows].sample = 0;

		for (short i = 0; i < 24; i++)
		{
			samples[total_rows].sample <<= 1;
			samples[total_rows].sample |= numbers[i] == '0' ? 0UL : 1UL;
		}

		/*
		printf("[%s] - %s - 0x%x\n", line, 
				samples[ns].filename, samples[ns].sample); 
		*/
		total_rows++;
	}

	/* may check feof here to make a difference between 
	   eof and io failure -- network
	   timeout for instance */

	fclose(file);

	/*
	 * Init current_row before starting up threads to process
	 */
	current_row = 0;

	pthread_t threads[MAX_THREADS];
	RowMatches rowMatches[MAX_THREADS];

	for (unsigned short i = 0; i < MAX_THREADS; i++)
	{
		int err = pthread_create(&(threads[i]), NULL, 
									&processRows, &(rowMatches[i]));
		if (err != 0)
			printf("\ncan't create thread :[%s]", strerror(err));
		else
			printf("\n Thread created successfully\n");
	}


 	unsigned long total_matches = 0;
	for (unsigned short i = 0; i < MAX_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
		total_matches += rowMatches[i].total;
	}

	fprintf(stderr, "\nFound %ld sample matches\n", total_matches);
	fprintf(stderr, "\nFound %ld row matches\n", total_matches);

	/*
	 *  Write out all the matches to stdout.
	 */
	char fileNames[128];
	char buffer[4096];
	buffer[0] = '\0';
	for (unsigned short i = 0; i < MAX_THREADS; i++)
	{
		for (unsigned j = 0; j < rowMatches[i].total; j++)
		{
			unsigned long r1 = rowMatches[i].matches[j].m1;
			unsigned long r2 = rowMatches[i].matches[j].m2;
			sprintf(fileNames, "%s %s\n", 
						samples[r1].filename, samples[r2].filename);
			strcat(buffer, fileNames);

			if (j % 75 == 0)
			{
				printf(buffer);
				buffer[0] = '\0';
			}

			/*
			printf("%s %s\n", 
						samples[r1].filename, samples[r2].filename);
			 */
		}
	}

	printf(buffer);
    return 0;

}
