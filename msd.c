#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct atomicCoordinates
{
	float x, y, z;
} ATOM_COORDS;

typedef struct chainEnds
{
	float x1, y1, z1, x2, y2, z2;
} CHAIN_ENDS;

typedef struct centerOfMass
{
	float x, y, z;
} CENTER_OF_MASS;

typedef struct coordinates
{
	float x, y, z;
} COORDINATES;

int calculateNatoms (FILE *input)
{
	rewind (input);
	int nAtoms;
	char lineString[2000];

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, input); }

	sscanf (lineString, "%d", &nAtoms);

	return nAtoms;
}

int calculateNtimeframes (FILE *input, int nAtoms)
{
	rewind (input);
	int nTimeframes, nLines = 0;
	char lineString[2000];

	while (fgets (lineString, 2000, input) != NULL)
	{
		nLines++;
	}

	nTimeframes = nLines / (nAtoms + 9);

	return nTimeframes;
}

// CHAIN_ENDS *fetchChainEnds (CHAIN_ENDS *ends, FILE *input, int nAtoms, int degreeOfPolymerization, int nChains)
// {
// 	char lineString[2000];
// 	int currentLine = 0, currentTimestep = 0;

// 	while (fgets (lineString, 2000, input) != NULL)
// 	{
// 		currentLine++;

// 		if (currentLine > 9 && currentLine <= 9 + nAtoms)
// 		{
// 			sscanf ()
// 		}
// 	}

// 	return ends;
// }

CENTER_OF_MASS **computeCenterOfMass (CENTER_OF_MASS **chainCOMs, FILE *input, int nAtoms, int degreeOfPolymerization, int nChains, ATOM_COORDS **allAtomicCoordinates)
{
	rewind (input);

	char lineString[2000];
	int currentLine = 0, currentTimestep = 0, currentArrayCounter = 0, coordEntryNumber = 0;

	COORDINATES *chainCoordinates;
	chainCoordinates = (COORDINATES *) malloc (nAtoms * sizeof (COORDINATES));

	while (fgets (lineString, 2000, input) != NULL)
	{
		currentLine++;

		if (currentLine > 9 && currentLine <= 9 + nAtoms)
		{
			sscanf (lineString, "%*d %*d %*f %*f %*f %f %f %f\n", &chainCoordinates[currentLine - 10].x, &chainCoordinates[currentLine - 10].y, &chainCoordinates[currentLine - 10].z);

			(*allAtomicCoordinates)[coordEntryNumber].x = chainCoordinates[currentLine - 10].x;
			(*allAtomicCoordinates)[coordEntryNumber].y = chainCoordinates[currentLine - 10].y;
			(*allAtomicCoordinates)[coordEntryNumber].z = chainCoordinates[currentLine - 10].z;

			coordEntryNumber++;
		}

		if (currentLine == (nAtoms + 9))
		{
			currentArrayCounter = 0;

			for (int j = 0; j < nChains; ++j)
			{
				chainCOMs[currentTimestep][j].x = 0;
				chainCOMs[currentTimestep][j].y = 0;
				chainCOMs[currentTimestep][j].z = 0;
			}

			for (int j = 0; j < nChains; ++j)
			{
				for (int k = 0; k < degreeOfPolymerization; ++k)
				{
					chainCOMs[currentTimestep][j].x += chainCoordinates[currentArrayCounter].x;
					chainCOMs[currentTimestep][j].y += chainCoordinates[currentArrayCounter].y;
					chainCOMs[currentTimestep][j].z += chainCoordinates[currentArrayCounter].z;

					currentArrayCounter++;
				}

				chainCOMs[currentTimestep][j].x /= degreeOfPolymerization;
				chainCOMs[currentTimestep][j].y /= degreeOfPolymerization;
				chainCOMs[currentTimestep][j].z /= degreeOfPolymerization;
			}

			currentLine = 0;
			currentTimestep++;
		}
	}

	return chainCOMs;
}

float calculateDisplacement (float x1, float y1, float z1, float x2, float y2, float z2)
{
	float displacement, dx, dy, dz;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;
	displacement = dx * dx + dy * dy + dz * dz; 
	return displacement;
}

float **computeMeanSquareDisplacement (float **meanSquareDisplacement, CENTER_OF_MASS **chainCOMs, int nTimeframes, int nChains)
{
	float displacement;
	int delT;
	int **msdDenominator;
	int progress = 0, progressCheckFreq = (int) ceil (nChains * nTimeframes * 0.01);
	msdDenominator = (int **) malloc (nTimeframes * sizeof (int *));

	for (int i = 0; i < nTimeframes; ++i) {
		msdDenominator[i] = (int *) malloc (nChains * sizeof (int)); }

	for (int i = 0; i < nTimeframes; ++i) {
		for (int j = 0; j < nChains; ++j) {
			msdDenominator[i][j] = 0; } }

	for (int i = 0; i < nTimeframes; ++i) {
		for (int j = 0; j < nChains; ++j) {
			meanSquareDisplacement[i][j] = 0; } }

	printf("Computing mean square displacement...\n");

	for (int currentChain = 0; currentChain < nChains; ++currentChain)
	{
		for (int initTime = 0; initTime < nTimeframes; ++initTime)
		{
			for (int finalTime = initTime; finalTime < nTimeframes; ++finalTime)
			{
				displacement = calculateDisplacement (chainCOMs[initTime][currentChain].x, chainCOMs[initTime][currentChain].y, chainCOMs[initTime][currentChain].z, chainCOMs[finalTime][currentChain].x, chainCOMs[finalTime][currentChain].y, chainCOMs[finalTime][currentChain].z);
				delT = finalTime - initTime;
				meanSquareDisplacement[delT][currentChain] += pow (displacement, 2);
				msdDenominator[delT][currentChain]++;
			}

			progress++;
			if ((progress % progressCheckFreq) == 0) {
				printf("%.0f%% completed               \r", (float) (progress * 100) /(nChains * nTimeframes));
				fflush (stdout); }
		}
	}
	printf("100%% completed...                    \n");

	for (int i = 0; i < nChains; ++i)
	{
		for (int j = 0; j < nTimeframes; ++j)
		{
			meanSquareDisplacement[j][i] /= msdDenominator[j][i];
		}
	}

	return meanSquareDisplacement;
}

float *avgMSD (float **meanSquareDisplacement, float *meanSquareDisplacement_avg, int nTimeframes, int nChains)
{
	// Initialize the values to zero
	for (int i = 0; i < nTimeframes; ++i)
	{
		meanSquareDisplacement_avg[i] = 0;
	}

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			meanSquareDisplacement_avg[i] += meanSquareDisplacement[i][j];
		}

		meanSquareDisplacement_avg[i] /= nChains;
	}
	return meanSquareDisplacement_avg;
}

void printMSD (float *meanSquareDisplacement_avg, int nTimeframes, const char *inputFilename)
{
	char *outputFilename;
	outputFilename = (char *) malloc (100 * sizeof (char));
	snprintf (outputFilename, 100, "%s_2.msd", inputFilename);

	FILE *output;
	output = fopen (outputFilename, "w");

	for (int i = 0; i < nTimeframes; ++i)
	{
		fprintf(output, "%d %f\n", i, meanSquareDisplacement_avg[i]);
	}

	fclose (output);
	free (outputFilename);
}

float *computeMeanSquareDisplacement2 (float *meanSquareDisplacement2, ATOM_COORDS *allAtomicCoordinates, int nAtoms, int nTimeframes)
{
	float scalarDistance;
	ATOM_COORDS *distance;
	distance = (ATOM_COORDS *) malloc (nAtoms * sizeof (ATOM_COORDS));

	// cycle through 'tau'
	for (int i = 1; i < floor (nTimeframes / 2); ++i)
	{
		// timestep
		for (int j = 0; j < nTimeframes; ++j)
		{
			// nAtoms
			for (int k = 0; k < nAtoms; ++k)
			{
				distance[k].x = allAtomicCoordinates
				distance[k].y = 
				distance[k].z = 
			}
		}
	}

	return meanSquareDisplacement2;
}

int main(int argc, char const *argv[])
{
	if (argc == 1) {
		(void)printf("\nARGUMENTS:\n~~~~~~~~~~\n\n   {~} argv[0] = program\n   {~} argv[1] = input dump filename\n   {~} argv[2] = Degree of polymerization\n   {~} argv[3] = Number of chains present in the dump file\n\nHOW THIS WORKS?\n~~~~~~~~~~~~~~~\n\n   The program will calculate the center of mass of each polymers, then compute an average mean square displacement as a function of time. All molecules of interest must be present initially, followed by other molecules (solvent/nanoparticles). Output is saved as \"argv[1].msd\". Column 1 in the output contains the si.no. for the time (multiply this number by your dump frequency and integration timestep) and column 2 contains the msd values. \n\nNOTE:\n~~~~~\n\n   Download a fresh copy of this code from https://github.com/raghurame/sameer before usage.\n\n");
		exit (1); }

	FILE *input;
	input = fopen (argv[1], "r");

	int degreeOfPolymerization = atoi (argv[2]), nChains = atoi (argv[3]);

	int nAtoms = calculateNatoms (input), nTimeframes = calculateNtimeframes (input, nAtoms);

	// Calculating MSD of center of mass
	CENTER_OF_MASS **chainCOMs;
	chainCOMs = (CENTER_OF_MASS **) malloc (nTimeframes * sizeof (CENTER_OF_MASS *));

	ATOM_COORDS *allAtomicCoordinates;
	allAtomicCoordinates = (ATOM_COORDS *) malloc (nAtoms * nTimeframes * sizeof (ATOM_COORDS));

	for (int i = 0; i < nTimeframes; ++i) {
		chainCOMs[i] = (CENTER_OF_MASS *) malloc (nChains * sizeof (CENTER_OF_MASS)); }

	chainCOMs = computeCenterOfMass (chainCOMs, input, nAtoms, degreeOfPolymerization, nChains, &allAtomicCoordinates);

	// for (int i = 0; i < (nAtoms * nTimeframes); ++i)
	// {
	// 	printf("%f %f %f\n", allAtomicCoordinates[i].x, allAtomicCoordinates[i].y, allAtomicCoordinates[i].z);
	// 	sleep (1);
	// }

	float **meanSquareDisplacement, *meanSquareDisplacement2;
	meanSquareDisplacement = (float **) malloc (nTimeframes * sizeof (float *));
	meanSquareDisplacement2 = (float *) malloc (nTimeframes * sizeof (float));

	for (int i = 0; i < nTimeframes; ++i) {
		meanSquareDisplacement[i] = (float *) malloc (nChains * sizeof (float)); }

	meanSquareDisplacement = computeMeanSquareDisplacement (meanSquareDisplacement, chainCOMs, nTimeframes, nChains);

	meanSquareDisplacement2 = computeMeanSquareDisplacement2 (meanSquareDisplacement2, allAtomicCoordinates, nAtoms, nTimeframes);

	float *meanSquareDisplacement_avg;
	meanSquareDisplacement_avg = (float *) malloc (nTimeframes * sizeof (float));

	meanSquareDisplacement_avg = avgMSD (meanSquareDisplacement, meanSquareDisplacement_avg, nTimeframes, nChains);

	printMSD (meanSquareDisplacement_avg, nTimeframes, argv[1]);

	free (chainCOMs);
	free (meanSquareDisplacement);
	free (meanSquareDisplacement_avg);
	fclose (input);
	return 0;
}