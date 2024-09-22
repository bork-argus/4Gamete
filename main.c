#include <stdlib.h>
#include <stdio.h>

// Set this to the maximum number of 64-bit storage elements you want to support for Alleles.
#define MAX_ALLELE_ELEMENTS		(2)

#define NUM_BITS_PER_BYTE		(8)
#define NUM_BITS_PER_ELEMENT	(sizeof(unsigned long long) * NUM_BITS_PER_BYTE)
#define MAX_ALLELES 			(MAX_ALLELE_ELEMENTS * NUM_BITS_PER_ELEMENT)

typedef struct {
	char 				locusName[128];
	unsigned long long 	allele[MAX_ALLELE_ELEMENTS];	// alleles array. Each element represents 64 alleles.
} Locus;

/*
int uniquePairs(unsigned long p1, unsigned long p2)
{
	unsigned char uniques[4] = { 0, 0, 0, 0 };

	register unsigned long bitMask = 1UL;
	for (unsigned short i = 0; i < 24; i++)
	{
		register unsigned long p1_bit = p1 & bitMask;
		register unsigned long p2_bit = p2 & bitMask;

		if (p1_bit && p2_bit)
			uniques[0] = 1;
		else if (p1_bit)
			uniques[1] = 1;
		else if (p2_bit)
			uniques[2] = 1;
		else
			uniques[3] = 1;

		bitMask <<= 1;
	}

	return uniques[0] + uniques[1] + uniques[2] + uniques[3];
}
*/

unsigned short uniquePairsViaLocus2(const Locus *locus1, const Locus *locus2, const unsigned short numBits)
{
	static const unsigned long long bitmasks[] = {
		1ULL, 1ULL << 1, 1ULL << 2, 1ULL << 3, 1ULL << 4, 1ULL << 5, 1ULL << 6, 1ULL << 7, 
		1ULL << 8, 1ULL << 9, 1ULL << 10, 1ULL << 11, 1ULL << 12, 1ULL << 16, 1ULL << 14, 1ULL << 15, 
		1ULL << 16, 1ULL << 17, 1ULL << 18, 1ULL << 19, 1ULL << 20, 1ULL << 21, 1ULL << 22, 1ULL << 23, 
		1ULL << 24, 1ULL << 25, 1ULL << 26, 1ULL << 27, 1ULL << 28, 1ULL << 29, 1ULL << 30, 1ULL << 31, 
		1ULL << 32, 1ULL << 33, 1ULL << 34, 1ULL << 35, 1ULL << 36, 1ULL << 37, 1ULL << 38, 1ULL << 39, 
		1ULL << 40, 1ULL << 41, 1ULL << 42, 1ULL << 43, 1ULL << 44, 1ULL << 45, 1ULL << 46, 1ULL << 47, 
		1ULL << 48, 1ULL << 49, 1ULL << 50, 1ULL << 51, 1ULL << 52, 1ULL << 53, 1ULL << 54, 1ULL << 55, 
		1ULL << 56, 1ULL << 57, 1ULL << 58, 1ULL << 59, 1ULL << 60, 1ULL << 61, 1ULL << 62, 1ULL << 63		
	};

	//----
	// Holds a record of now many unique pairs of bits we'd seen in the low four bits.
	//----
	unsigned char uniques = 0;

	register char index = 0;
	// register unsigned long long bitmask = 1ULL;
	register unsigned long long p1 = locus1->allele[index], p2 = locus2->allele[index];
	register const unsigned long long *bitmask = bitmasks;
	for (unsigned short i = 0; i < numBits; i++)
	{
		// register unsigned long long p1_bit = p1 & bitmask;
		// register unsigned long long p2_bit = p2 & bitmask;
		register unsigned long long p1_bit = p1 & *bitmask;
		register unsigned long long p2_bit = p2 & *bitmask;

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
		// bitmask <<= 1;

		// if we've managed to rotate the bitmask all the way off the
		// high end, we've exhausted how may bits are available in the
		// current element. Advance to the next index and keep comparing.
		// if (bitmask == 0ULL)
		if (*bitmask == 1ULL << 63)
		{
			index++;
			// bitmask = 1ULL;
			bitmask = bitmasks;
			p1 = locus1->allele[index];
			p2 = locus1->allele[index];
		}
		else
			bitmask++;
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

	// This switch statement is a touch faster than looping and adding
	// to calculate the number of uniques.
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

	// This switch statement is a touch faster than looping and adding
	// to calculate the number of uniques.
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
		if (*p == ' ' || *p == '	')
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

int main(int argc, char* argv[])
{
	// char const* const fileName = "binary_alleles_per_var.txt";
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <input allele file data>\n", argv[0]);
		exit(-1);
	}

    char const* const fileName = argv[1];
	FILE* file = fopen(fileName, "r"); /* should check the result */
	char line[256];
	char filename[128];

	Locus *samples;

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
	samples = calloc(numLoci, sizeof(Locus));

	// reset file pointer to beginning to file
	rewind(file);

	//----
	// Read in the first line to determine the number of
	// allenes for this file.
	//----
	fgets(line, sizeof(line), file);
	const unsigned short nAlleles = numAlleles(line);

	rewind(file);

	//----
	// Read in the file (again) converting the
	// lines as we go into Locus structures and
	// storing them in our samples array.
	//----
	unsigned long ns = 0;
	while (fgets(line, sizeof(line), file)) 
	{
		/* note that fgets don't strip the terminating \n, checking its
		   presence would allow to handle lines longer that sizeof(line) */

		// this should bit copy the parsed Locus structure into our array of
		// loci.
		samples[ns] = parseLocus(line);

		ns++;
	}

	/* may check feof here to make a difference between 
	   eof and io failure -- network
	   timeout for instance */

	fclose(file);

	/* samples are in samples[] ... run a double loop to compare */
	/* ns holds the number of samples read */
	
	unsigned long uniques = 0;
	fprintf(stderr, "Processing %ld samples\n", ns);
	for (unsigned long int row = 0; row < ns - 1; row++)
	{
		if (row % 100 == 0)
			fprintf(stderr, ".");

		for (unsigned long int cmp_row = row + 1; cmp_row < ns; cmp_row++)
		{
			if (uniquePairsViaLocus(&(samples[row]), &(samples[cmp_row]), nAlleles) == 4)
			{
				uniques++;
				printf("%s %s\n", 
						samples[row].locusName, samples[cmp_row].locusName);
			}
		}
	}

	fprintf(stderr, "\nFound %ld sample matches\n", uniques);

    return 0;

}
