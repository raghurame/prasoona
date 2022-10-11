#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>

int findNentries (FILE *input)
{
	rewind (input);
	int nEntries;

	char lineString[2000];

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineString, 2000, input);
	}

	sscanf (lineString, "%d", &nEntries);

	return nEntries;
}

int findNtimeframes (FILE *input, int nEntries)
{
	rewind (input);
	int nTimeframes = 0, nLines = 0;
	char lineString[2000];

	while (fgets (lineString, 2000, input) != NULL)
	{
		nLines++;
	}

	nTimeframes = nLines / (nEntries + 9);

	return nTimeframes;
}

float computeBondAverage (FILE *input, int nEntries, int nTimeframes, float bondAverage)
{
	rewind (input);
	char lineString[2000];
	float bondAverage_local;
	bondAverage = 0;
	int denominator = 0;

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int i = 0; i < 9; ++i)
		{
			fgets (lineString, 2000, input);
		}

		for (int i = 0; i < nEntries; ++i)
		{
			fgets (lineString, 2000, input);
			sscanf (lineString, "%*d %*d %*d %*f %f", &bondAverage_local);
			denominator++;
			bondAverage += bondAverage_local;
		}
	}

	bondAverage /= denominator;

	return bondAverage;
}

int main(int argc, char const *argv[])
{
	FILE *input;
	input = fopen (argv[1], "r");

	int nEntries = findNentries (input), nTimeframes = findNtimeframes (input, nEntries);

	// printf("%d %d\n", nEntries, nTimeframes);
	float bondAverage;
	bondAverage = computeBondAverage (input, nEntries, nTimeframes, bondAverage);

	printf("Bond average: %f\n", bondAverage);

	return 0;
}