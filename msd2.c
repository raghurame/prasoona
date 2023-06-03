#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#define NPOLYMERBEADS 1
#define NCOUNTERIONS 0
#define OUTPUTFILENAME "msd.output"
#define INVERSE_TIMESTEP_LIMIT 1

typedef struct trajdata
{
	float x, y, z;
	int type, ix, iy, iz;
} TRAJ_DATA;

int countNAtoms (FILE *inputDumpFile, char lineContent[])
{
	int nAtoms;

	for (int i = 0; i < 4; ++i)
	{
		fgets (lineContent, 2000, inputDumpFile);
	}

	sscanf (lineContent, "%d\n", &nAtoms);
	rewind (inputDumpFile);
	printf("%d atoms are present...\n", nAtoms);

	return nAtoms;
}

int countTimeframes (const char *filename, int nAtoms)
{
	int totalTimesteps;
	FILE *countingLines;
	char *pipeString, lineContent[2000];
	pipeString = (char *) malloc (100 * sizeof (char));
	snprintf (pipeString, 100, "wc -l %s", filename);
	countingLines = popen (pipeString, "r");

	fgets (lineContent, 2000, countingLines);
	sscanf (lineContent, "%d", &totalTimesteps);

	totalTimesteps /= (nAtoms + 9);
	printf("Identified %d timesteps in traj file...\n", totalTimesteps);

	return totalTimesteps;
}

TRAJ_DATA **readTrajData (TRAJ_DATA **atoms, FILE *inputDumpFile, int totalTimesteps, int nAtoms)
{
	printf("Saving trajectory data...\n");
	char lineContent[2000];
	rewind (inputDumpFile);

	for (int i = 0; i < totalTimesteps; ++i)
	{
		for (int j = 0; j < 9; ++j)
		{
			fgets (lineContent, 2000, inputDumpFile);
		}
		for (int j = 0; j < nAtoms; ++j)
		{
			fgets (lineContent, 2000, inputDumpFile);
			sscanf (lineContent, "%*d %d %f %f %f %d %d %d\n", 
				&atoms[i][j].type, 
				&atoms[i][j].x, 
				&atoms[i][j].y, 
				&atoms[i][j].z, 
				&atoms[i][j].ix, 
				&atoms[i][j].iy, 
				&atoms[i][j].iz);
		}
	}

	return atoms;
}

TRAJ_DATA **initializeAtoms (TRAJ_DATA **atoms, int totalTimesteps, int nAtoms)
{
	printf("Initializing mem for atom trajectory...\n");
	atoms = (TRAJ_DATA **) malloc (totalTimesteps * sizeof (TRAJ_DATA *));

	for (int i = 0; i < totalTimesteps; ++i)
	{
		atoms[i] = (TRAJ_DATA *) malloc (nAtoms * sizeof (TRAJ_DATA));
	}

	return atoms;
}

float **initializeDisplacementSquare (float **displacementSquare, int totalTimesteps, int nAtoms)
{
	printf("Initializing mem for storing displacement square...\n");
	displacementSquare = (float **) malloc (totalTimesteps * sizeof (float *));

	for (int i = 0; i < totalTimesteps; ++i)
	{
		displacementSquare[i] = (float *) calloc (nAtoms, sizeof (float));
	}

	return displacementSquare;
}


float **computeDisplacementSquare (float **displacementSquare, TRAJ_DATA **atoms, int totalTimesteps, int nAtoms)
{
	int currentProgress = 0, **denominator;
	denominator = (int **) malloc (totalTimesteps * sizeof (int *));

	for (int i = 0; i < totalTimesteps; ++i) {
		denominator[i] = (int *) calloc (nAtoms, sizeof (int)); }

	#pragma omp parallel for
	for (int i = 0; i < ceil ((float)totalTimesteps / (float)INVERSE_TIMESTEP_LIMIT); ++i)
	{
		currentProgress++;
		if ((currentProgress%100) == 0)
		{
			fprintf(stdout, "Calculating displacement square: %.2f%% complete...           \r", ((float)currentProgress * (float)INVERSE_TIMESTEP_LIMIT * (float)100 / (float)totalTimesteps), currentProgress);
			fflush (stdout);
		}

		for (int j = 0; j < (totalTimesteps - i); ++j)
		{
			for (int k = 0; k < nAtoms; ++k)
			{
				displacementSquare[i][k] += 
					((atoms[j + i][k].x - atoms[j][k].x) * (atoms[j + i][k].x - atoms[j][k].x)) + 
					((atoms[j + i][k].y - atoms[j][k].y) * (atoms[j + i][k].y - atoms[j][k].y)) + 
					((atoms[j + i][k].z - atoms[j][k].z) * (atoms[j + i][k].z - atoms[j][k].z));
				denominator[i][k]++;
			}
		}
	}

	for (int i = 0; i < totalTimesteps; ++i)
	{
		for (int j = 0; j < nAtoms; ++j) {
			displacementSquare[i][j] /= denominator[i][j]; }
	}

	return displacementSquare;
}


int main(int argc, char const *argv[])
{
	if (argc != 4)
	{
		printf("Incorrect args passed. The program requires input dump file as argv[1], dt as argv[2], and output frequency as argv[3].\n");
		exit (1);
	}

	FILE *inputDumpFile, *outputMSDFile;
	inputDumpFile = fopen (argv[1], "r");
	outputMSDFile = fopen (OUTPUTFILENAME, "w");

	char lineContent[2000];

	int nAtoms = countNAtoms (inputDumpFile, lineContent), totalTimesteps = countTimeframes (argv[1], nAtoms);

	TRAJ_DATA **atoms;
	atoms = initializeAtoms (atoms, totalTimesteps, nAtoms);
	atoms = readTrajData (atoms, inputDumpFile, totalTimesteps, nAtoms);

	float **displacementSquare;
	displacementSquare = initializeDisplacementSquare (displacementSquare, totalTimesteps, nAtoms);
	displacementSquare = computeDisplacementSquare (displacementSquare, atoms, totalTimesteps, nAtoms);

	float *meanSquareDisplacementPolymer, *meanSquareDisplacementCounterions;
	meanSquareDisplacementPolymer = (float *) calloc (totalTimesteps, sizeof (float));
	meanSquareDisplacementCounterions = (float *) calloc (totalTimesteps, sizeof (float));

	for (int i = 0; i < totalTimesteps; ++i)
	{
		for (int j = 0; j < NPOLYMERBEADS; ++j)
		{
			meanSquareDisplacementPolymer[i] += displacementSquare[i][j];
		}

		meanSquareDisplacementPolymer[i] /= 92;
	}

	for (int i = 0; i < totalTimesteps; ++i)
	{
		for (int j = 93; j < NCOUNTERIONS; ++j)
		{
			meanSquareDisplacementCounterions[i] += displacementSquare[i][j];
		}

		meanSquareDisplacementCounterions[i] /= 55;
	}

	for (int i = 0; i < ceil ((float)totalTimesteps / (float)INVERSE_TIMESTEP_LIMIT); ++i)
	{
		fprintf(outputMSDFile, "%f %f %f\n", (i + 1) * atof (argv[2]) * atof (argv[3]), meanSquareDisplacementPolymer[i], meanSquareDisplacementCounterions[i]);
	}

	fclose (inputDumpFile);
	fclose (outputMSDFile);
	return 0;
}