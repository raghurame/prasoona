#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

int countTimeframes (FILE *inputfile, const char *filename, int nAtoms)
{
	rewind (inputfile);
	FILE *wordCountPipe;
	char wordCountPipeString[200], lineString[2000];
	int nLines = 0;
	snprintf (wordCountPipeString, 200, "wc -l %s", filename);
	wordCountPipe = popen (wordCountPipeString, "r");
	fgets (lineString, 2000, wordCountPipe);
	sscanf (lineString, "%d\n", &nLines);
	int nTimeframes = nLines / (nAtoms + 9);
	rewind (inputfile);
	return nTimeframes;
}

int countNAtoms (FILE *inputfile)
{
	rewind (inputfile);
	int nAtoms = 0;
	char lineString[2000];
	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, inputfile); }
	sscanf (lineString, "%d\n", &nAtoms);
	rewind (inputfile);
	return nAtoms;
}

typedef struct atomposition
{
	float x, y, z;
} ATOM_POSITION;

float calculateDistance (ATOM_POSITION atom1, ATOM_POSITION atom2)
{
	float distance;

	distance = sqrt (pow ((atom2.x - atom1.x), 2) + pow ((atom2.y - atom1.y), 2) + pow ((atom2.z - atom1.z), 2));

	return distance;
}

float *computeEndToEndDistance (FILE *inputfile, int nTimeframes, int nAtoms, int chainEnd1, int chainEnd2, float *endToEndDistance)
{
	rewind (inputfile);
	char lineString[2000];
	ATOM_POSITION atom1, atom2;
	int sino;

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < 9; ++j) {
			fgets (lineString, 2000, inputfile); }
		for (int j = 0; j < nAtoms; ++j)
		{
			fgets (lineString, 2000, inputfile);
			sscanf (lineString, "%d\n", &sino);
			if (sino == chainEnd1) {
				sscanf (lineString, "%*d %*d %f %f %f\n", &atom1.x1, &atom1.y1, &atom1.z1); }
			else if (sino == chainEnd2) {
				sscanf (lineString, "%*d %*d %f %f %f\n", &atom2.x2, &atom2.y2, &atom2.z2); }
		}
		endToEndDistance[i] = calculateDistance (atom1, atom2);
	}

	rewind (inputfile);
	return endToEndDistance;
}

int main(int argc, char const *argv[])
{
	FILE *inputfile;
	inputfile = fopen (argv[1], "r");

	int nAtoms = countNAtoms (inputfile), nTimeframes = countTimeframes (inputfile, argv[1], nAtoms), chainEnd1 = atoi (argv[2]), chainEnd2 = atoi (argv[3]);
	float *endToEndDistance;
	endToEndDistance = (float *) malloc (nTimeframes * sizeof (float));

	endToEndDistance = computeEndToEndDistance (inputfile, nTimeframes, nAtoms, chainEnd1, chainEnd2, endToEndDistance);

	free (endToEndDistance);
	fclose (inputfile);
	return 0;
}