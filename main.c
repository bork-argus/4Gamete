#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>

//----
// Set MAX_ALLELE_ELEMENTS to the maximum number of 64-bit storage elements you want to support for Alleles.
// Two elements supports up to 128 alleles/locus.
// Three elements would support up to 192 alleles/locus.
//----

// Maximum number of allele elements we store per locus
// Set this to the maximum number of 64-bit storage elements you want to support for Alleles.
#define MAX_ALLELE_ELEMENTS		(2)

//----
// Various macros to define the allele representation. Do not change any of these.
//----

// Number of bits per byte from limits.h
#define NUM_BITS_PER_BYTE		(CHAR_BIT)

// Maximum number of bits per allele element representation
#define NUM_BITS_PER_ELEMENT	(sizeof(unsigned long long) * NUM_BITS_PER_BYTE)

// Maximum number of alleles per locus we can handle
#define MAX_ALLELES 			(MAX_ALLELE_ELEMENTS * NUM_BITS_PER_ELEMENT)

//----
// Defines the temporary storage for each line read from the input file.
// If any will exceed this value, it will need to be raised.
// Keep in mind this is allocated on the stack, so, if it is raised
// beyond, say 1024, it might be wiser to dynamically allocated the space
// in the heap instead.
//----
#define MAX_LINE_LENGTH			(384)

//----
// Command line options
//----
static unsigned char l_no_output = 0;			// -n
static unsigned short l_num_threads = 1;		// -j <num>

// defines the maximum number of matching indices we can store in a LocusMatches structure.
#define LOCUS_INDEX_CHUNK		(1024)

// This structure is used to record matching indices into the sample array.
// Up to LOCUS_INDEX_CHUNK indexes can be stored in the structure and additional
// structures, in the case of overflow, are accessed through next_chunk.
typedef struct {
	unsigned long			match_index[LOCUS_INDEX_CHUNK];
	unsigned short			num_matches;
	struct LocusMatches		*next_chunk;
} LocusMatches;

// This structure defines a "locus", including its name and alleles stored as a
// bit map. It also contains a pointer to any matching samples (in field matches).
typedef struct {
	char 				locusName[128];
	unsigned long long 	allele[MAX_ALLELE_ELEMENTS];	// alleles array. Each element represents 64 alleles.
	LocusMatches		*matches;
} Locus;

/// @brief Compares numBits alleles associated with two loci, returning the number of unique matches found (1 - 4)
/// @param locus1 First locus to compare
/// @param locus2 Second locus to compare
/// @param numBits Number of alleles to compare in the loci
/// @return The number of unique matches in the alleles -- can only be up to 4.
unsigned short uniquePairsViaLocus(const Locus *locus1, const Locus *locus2, const unsigned short numBits)
{
	//----
	// Holds a record of now many unique pairs of bits we'd seen in the low four bits.
	//----
	unsigned char uniques = 0;

	register char index = 0;
	register unsigned long long bitmask = 1ULL;
	register unsigned long long p1 = locus1->allele[index], p2 = locus2->allele[index];
	for (unsigned short i = 0; i < numBits; i++)
	{
		register unsigned long long p1_bit = p1 & bitmask;
		register unsigned long long p2_bit = p2 & bitmask;

		if (p1_bit && p2_bit)
			uniques |= 0x01;
		else if (p1_bit)
			uniques |= 0x02;
		else if (p2_bit)
			uniques |= 0x04;
		else
			uniques |= 0x08;

		// no point in continuing if we already have 4 unique values.
		if (uniques == 0x0F)
			break;

		// check the next bit
		bitmask <<= 1;

		// if we've managed to rotate the bitmask all the way off the
		// high end, we've exhausted how may bits are available in the
		// current element. Advance to the next index and keep comparing.
		if (bitmask == 0ULL)
		{
			index++;
			bitmask = 1ULL;
			p1 = locus1->allele[index];
			p2 = locus1->allele[index];
		}
	}

	//---
	// Convert the record of number of unique pairs in the
	// bit pattern of uniques into numUniques. Basically, we
	// are counting the number of set bits in the lowest
	// four bits of uniques.
	//----
	unsigned short numUniques = 0;
	/*
	for (register unsigned char i = 0; i < 4; i++)
	{
		numUniques += uniques & 0x01;
		uniques >>= 1;
	}
	*/

	//----
	// This switch statement is a touch faster than looping and adding
	// to calculate the number of uniques.
	//----
	switch (uniques & 0x0F)
	{
		// This case should never occur
		case 0x0:
			break;

		// Cases where exactly 1 bit is set
		case 0x01:
		case 0x02:
		case 0x04:
		case 0x08:
			numUniques = 1;
			break;

		// All bits set	
		case 0x0F:
			numUniques = 4;
			break;
		
		// Exactly 2 bits set
		case 0x03:
		case 0x05:
		case 0x06:
		case 0x09:
		case 0x0A:
		case 0x0C:
			numUniques = 2;
			break;

		// Exactly 3 bits set
		case 0x07:
		case 0x0D:
		case 0x0E:
		case 0x0B:
			numUniques = 3;
			break;
		
		// This should be utterly impossible to get to
		default:
			// error?
			break;
	}

	return numUniques;
}

/// @brief Returns the number of alleles in a locus line
/// @param sampleLocus Locus line from input file
/// @return 
unsigned short numAlleles(const char *locusLine)
{
	// walk the sampleLocus until we find white space, after that,
	// count the number of alleles until we hit the end of the line
	register unsigned short numAlleles = 0;
	register unsigned char latch = 0x00;
	for (register const char *p = locusLine; *p; p++)
	{
		// find first separator between the locus name and the alleles
		if (*p == ' ' || *p == '\t')
			latch = 0x01;

		// count how many alleles we see, but only after the locus name
		if (latch && (*p == '0' || *p == '1'))
			numAlleles++;
	}

	return numAlleles;
}

/// @brief Parses a input file line into a Locus structure
/// @param locusLine 
/// @return 
Locus parseLocus(const char *locusLine)
{
	Locus locus;

	locus.matches = NULL;

	register const char *p = locusLine;
	register char *ln = locus.locusName;

	// get locus name
	for (; *p != ' ' && *p != '\t'; p++)
	{
		*ln++ = *p;
	}
	*ln = '\0';  // terminate the locus name.

	// initialize allelles
	for (register unsigned char i = 0; i < MAX_ALLELE_ELEMENTS; i++)
		locus.allele[i] = 0;

	unsigned short numAlleles;
	for (; *p; p++)
	{
		if (*p == '0' || *p == '1')
		{
			unsigned char index = numAlleles/NUM_BITS_PER_ELEMENT;
			locus.allele[index] <<= 1;		// shift to left
			locus.allele[index] |= (*p == '0' ? 0 : 1); // OR in the new allele
			numAlleles++;
		}
	}

	return locus;
}

/// @brief Adds a matching index to locus.
/// @param locus Locus structure to add match_index to
/// @param match_index Macthing index
void AddMatch(Locus *locus, unsigned long match_index)
{
	LocusMatches *p = locus->matches;

	//----
	// if this is the first time through, allocate the first match structure.
	//----
	if (p == NULL)
	{
		p = (LocusMatches *) calloc(1, sizeof(LocusMatches));
		if (p == NULL)
		{
			fprintf(stderr, "Failed to allocate memory to store locus matches. Aborting.\n\n");
			exit(-1);
		}

		locus->matches = p;
	}

	//----
	// walk the matches list until we find either a structure with room,
	// or a structure with no room and a NULL pointer.
	// if there is room, then just add the value into the existing match structure,
	// if there isn't room, then allocate a new one and add it in.
	//----
	while (p)
	{
		//----
		// Not full of match indices? Add her in and move on.
		//----
		if (p->num_matches < LOCUS_INDEX_CHUNK)
		{
			p->match_index[p->num_matches] = match_index;
			p->num_matches++;
			break;
		}
		else
		{
			//----
			// walk to the next one, if we are at the end, allocate a new one and 
			// add it into the list and then keep going.
			//----
			if (p->next_chunk == NULL)
			{
				p->next_chunk = (struct LocusMatches *) calloc(1, sizeof(LocusMatches));
				if (p->next_chunk == NULL)
				{
					fprintf(stderr, "Failed to allocate memory to store next locus matches. Aborting.\n\n");
					exit(-1);
				}
			}
			
			p = (LocusMatches *) p->next_chunk;
		}
	}
}

/// @brief Holds information about an input file
typedef struct
{
	Locus			*loci;
	unsigned long	num_loci;
	unsigned short	num_alleles;
} InputFileInfo;

#define IS_INPUT_VALID(input_file_info)		(input_file_info.loci != NULL)

//----
// Local variables for thread processing
// This will be get populated by the main thread and then
// passes to each child thread.
//----
typedef struct {
	InputFileInfo	file1;
	InputFileInfo	file2;
	unsigned long	current_row;
	pthread_mutex_t current_row_mutex;
} ThreadArgs;

/// @brief Thread to process loci rows
/// @param arg 
/// @return 
void *processLoci(void *arg)
{
	ThreadArgs *args = (ThreadArgs *) arg;

	const unsigned long process_rows = 100;

	Locus *samples = args->file1.loci;
	unsigned long total_samples = args->file1.num_loci;
	unsigned short nAlleles = args->file1.num_alleles;
	pthread_mutex_t *mutex = &(args->current_row_mutex);

	unsigned char file2_valid = IS_INPUT_VALID(args->file2);
	Locus *samples2 = args->file2.loci;
	unsigned long total_samples2 = args->file2.num_loci;

	//----
	// There is different end-game processing if we are handling
	// one file vs. two.
	// If we are processing only File1, then we can't compare the
	// last row (index = total_samples - 1) to anything -- there
	// are no rows after it to compare against, so ... we stop
	// if the starting row is index (total_samples - 1) or greater.
	// If we are comparing against a second file, then we need to
	// process all rows in file1 (i.e. until starting row advances
	// to total_samples).
	//----
	unsigned char last_index_from_end = 1;
	if (file2_valid)
	{
		//----
		// We need to process all rows in File1 so we process until
		// start_row is actually invalid (i.e. larger or equal to total_samples)
		//----
		last_index_from_end = 0;
	}

	const unsigned long final_processing_row = total_samples - last_index_from_end;

	while (1)
	{
		// only one thread can access current_row at a time
		// Everything else during processing ought to be thread safe.
		pthread_mutex_lock(mutex);
		unsigned long start_row = args->current_row;
		args->current_row += process_rows;
		pthread_mutex_unlock(mutex);
				
		// if the thread detects we are done ... exit
		if (start_row >= final_processing_row)
			break;

		//----
		// calculate the last row to process. If we fall off the end, adjust
		// back to the index beyond the last row we should touch.
		//----
		unsigned long end_row = start_row + process_rows;
		if (end_row > final_processing_row)
			end_row = final_processing_row;

		if (file2_valid)
		{
			//----
			// Compare each row in our current block of rows to every sample
			// appearing in the second file.
			//----
			for (register unsigned long row = start_row; row < end_row; row++)
			{
				for (register unsigned long cmp_row = 0; cmp_row < total_samples2; cmp_row++)
				{
					if (uniquePairsViaLocus(&(samples[row]), &(samples2[cmp_row]), nAlleles) == 4)
					{
						AddMatch(&(samples[row]), cmp_row);
					}				
				}
			}
		}
		else
		{
			//----
			// Compare each row in our current block of rows to every sample
			// appearing in the array after it. Yes, this is an exponential
			// comparison and the sheer number of comparisons will make this a
			// lengthy process.
			//----
			for (register unsigned long row = start_row; row < end_row; row++)
			{
				for (register unsigned long cmp_row = row + 1; cmp_row < total_samples; cmp_row++)
				{
					if (uniquePairsViaLocus(&(samples[row]), &(samples[cmp_row]), nAlleles) == 4)
					{
						AddMatch(&(samples[row]), cmp_row);
					}				
				}
			}
		}

		//----
		// print progress every process_rows lines
		//----
		fprintf(stderr, ".");
	}

	return NULL;
}

/// @brief Reads the input file and returns the Locus samples
/// @param filename 
/// @return Structure with loci and stats
InputFileInfo ReadInputFile(const char *fileName)
{
	// File information structure.
	InputFileInfo		fileInfo;

	// Initialize the file information
	fileInfo.loci = NULL;
	fileInfo.num_loci = 0;
	fileInfo.num_alleles = 0;

	//----
	// If filename is unspecified, bail.
	//----
	if (fileName == NULL || *fileName == '\0')
	{
		return fileInfo;
	}

	//----
	// Open the input file
	//----
	FILE* file = fopen(fileName, "r");
	if (file == NULL)
	{
		fprintf(stderr, "Could not open input file: %s\n", fileName);
		exit(-1);
	}

	//----
	// This will break if any lines in the input file are greater than
	// MAX_LINE_LENGTH. If that happens, raise MAX_LINE_LENGTH and
	// recompile.
	//----
	char line[MAX_LINE_LENGTH];

	//----
	// read in the file, discarding the contents as we go.
	// We are only interested in the number of records in the
	// file at this point.
	//----
	unsigned long numLoci = 0;
	while (fgets(line, sizeof(line), file))
		numLoci++;

	//----
	// Allocate space enough for all the records in the file.
	//----
	fileInfo.loci = (Locus *) calloc(numLoci, sizeof(Locus));
	if (fileInfo.loci == NULL)
	{
		fprintf(stderr, "Failed to allocate memory to store %s. Aborting.", fileName);
		exit(-1);
	}

	// reset file pointer to beginning to file
	rewind(file);

	//----
	// Read in the first line to determine the number of
	// allenes for this file.
	//----
	fgets(line, sizeof(line), file);
	fileInfo.num_alleles = numAlleles(line);
	if (fileInfo.num_alleles > MAX_ALLELES)
	{
		fprintf(stderr, "Input file [%s] specifies too many alleles (%d) which is larger than "
							"the maximum: [%ld]\n\n"
							"Either increase MAX_ALLELE_ELEMENTS or reduce alleles in the input file.", 
							fileName, fileInfo.num_alleles, MAX_ALLELES);
		exit(-1);
	}

	rewind(file);

	//----
	// Read in the file (again) converting the
	// lines as we go into Locus structures and
	// storing them in our samples array.
	//----
	unsigned long ns = 0;
	Locus *samples = fileInfo.loci;
	while (fgets(line, sizeof(line), file)) 
	{
		/* note that fgets don't strip the terminating \n, checking its
		   presence would allow to handle lines longer that sizeof(line) */

		// this should bit copy the parsed Locus structure into our array of
		// loci.
		samples[ns] = parseLocus(line);

		ns++;
	}

	fileInfo.num_loci = ns;

	/* may check feof here to make a difference between 
	   eof and io failure -- network
	   timeout for instance */

	fclose(file);

	return fileInfo;
}

/// @brief Entry point of application.
/// @param argc 
/// @param argv 
/// @return 
int main(int argc, char* argv[])
{
	//----
	// Handle command line options
	// -n = no output
	// -j <num> = number of threads to launch -- set to number of cores in processor for optimal results.
	//----
    int c;
    extern char *optarg;
    extern int optind, opterr, optopt;

    opterr = 0; // disable error messages

    while ((c = getopt(argc, argv, "nj:")) != -1) {
        switch (c) {
			case 'n':
				l_no_output = 1;
				break;
			case 'j':
				l_num_threads = atoi(optarg);
				break;
			case '?':
				printf("Unknown option %c\n", optopt);
				break;
        }
    }

	if (argc <= optind)
	{
		fprintf(stderr, "Usage: %s <input allele file data>\n", argv[0]);
		exit(-1);
	}

	//----
	// Read the input files specified on the command line
	//----
	InputFileInfo file1 = ReadInputFile(argv[optind]);
	InputFileInfo file2 = ReadInputFile(argc >= optind + 1 ? argv[optind + 1] : NULL);

	//----
	// Check if the two files match up with number of alleles
	// If not, bail.
	//----
	if (IS_INPUT_VALID(file2) && (file1.num_alleles != file2.num_alleles))
	{
		fprintf(stderr, "File %s number of alleles (%d) mismatches with file %s (%d)\n\n",
				argv[optind], file1.num_alleles, argv[optind + 1], file2.num_alleles);
		exit(-1);
	}
		
	fprintf(stderr, "Processing %ld samples\n", file1.num_loci);

	struct timeval t0, t1;
	gettimeofday(&t0, NULL);

	//----
	// peel off threads to handle the comparisons
	// ideally, -j on the command line set the number of threads to be
	// the number of processor cores available. We'll use each core to
	// simultaneously process various parts of the samples in parallel.
	// See thread function processLoci().
	//----
	pthread_t *thread = (pthread_t *) calloc(l_num_threads, sizeof(pthread_t));

	//----
	// This communicates various processing values to the threads.
	// Of course, each thread must know specifics about the data.
	// We pass samples, the number of samples, the number of alleles
	// per locus, and the current row for the threads to work on.
	// Note that the current_row will be advanced by each thread
	// as they do their work. Because each thread will advance this
	// value, it is very important to protect current_row from
	// race conditions. Thus the mutex which is used to ensure that
	// only one thread at a time accesses current_row (and more importantly
	// can modify it).
	//----
	ThreadArgs thread_args;
	thread_args.file1 = file1;
	thread_args.file2 = file2;
	thread_args.current_row = 0;
	pthread_mutex_init(&(thread_args.current_row_mutex), NULL);

	//----
	// Spin up the number of threads requested. As soon as the
	// thread is created, it will begin to do work.
	// Store the thread identifier in the thread array.
	// We'll use these identifiers to wait until all threads
	// are done before continuing here.
	//----
	for (unsigned short i = 0; i < l_num_threads; i++)
	{
		int err = pthread_create(&(thread[i]), NULL, 
								&processLoci, &thread_args);
		if (err != 0)
			fprintf(stderr, "\n Can't create thread :[%s]", strerror(err));
		else
			fprintf(stderr, "\n Thread created successfully.\n");
	}

	//----
	// wait for all the threads to finish
	// This loop cannot complete until all threads complete.
	//----
	for (unsigned short i = 0; i < l_num_threads; i++)
	{
		pthread_join(thread[i], NULL);
	}

	free(thread);

	gettimeofday(&t1, NULL);
	double elapsed = (double)(t1.tv_usec - t0.tv_usec) / 1000000 + (double)(t1.tv_sec - t0.tv_sec);
	fprintf(stderr, "\n\nElapsed time to perform comparisons: %f seconds\n", elapsed);

	//----
	// Walk the locus list and print matches (if output is allowed)
	// The processing threads identify matching pairs of locus and store
	// this information in structures attached to each locus in our samples.
	// Each block of matches can handle 1024 matches, but sometimes that might
	// not be enough, so, if it overflows, the system allocates another
	// block of 1024 slots and stores a pointer to it in next_chunk.
	// This block of code walks through this data structure and decodes the
	// matches, counting them and potentially printing to stdout as it
	// goes.
	//-----
	gettimeofday(&t0, NULL);

	Locus *samples = file1.loci;

	unsigned long uniques = 0;
	for (register unsigned long locus_index = 0; locus_index < file1.num_loci; locus_index++)
	{
		LocusMatches *p = samples[locus_index].matches;
		while (p)
		{
			for (register unsigned short i = 0; i < p->num_matches; i++)
			{
				uniques++;
				if (!l_no_output)
				{
					unsigned long match_index = p->match_index[i];
					printf("%s %s\n", 
							samples[locus_index].locusName, 
							IS_INPUT_VALID(file2) ? file2.loci[match_index].locusName : samples[match_index].locusName);
				}
			}

			//----
			// We no longer need the memory here, so free it up, but before we do, grab the next pointer.
			//----
			LocusMatches *next = (LocusMatches *) p->next_chunk;
			free(p);

			p = next;
		}
	}

	//----
	// Print out various statistics about the processing at the end.
	//----
	gettimeofday(&t1, NULL);
	elapsed = (double)(t1.tv_usec - t0.tv_usec) / 1000000 + (double)(t1.tv_sec - t0.tv_sec);
	fprintf(stderr, "\nElapsed time to write results: %f seconds\n", elapsed);

	fprintf(stderr, "\nFound %ld sample matches\n", uniques);

	fprintf(stderr, "\nStats:\n");

	if (!IS_INPUT_VALID(file2))
	{
		double matching_percentage = (((double) uniques * 2) / (((double) file1.num_loci) * ((double) file1.num_loci - 1))) * 100;
		fprintf(stderr, "Matching percentage: %.3f%%\n", matching_percentage);
	}

	fprintf(stderr, "\n");

	//----
	// Free the sample storage
	// Technically, we are killing the process anyway, the memory would get freed
	// regardless, but let's be cleaner about this.
	//----
	free(file1.loci);
	if (IS_INPUT_VALID(file2))
		free(file2.loci);

    return 0;
}
