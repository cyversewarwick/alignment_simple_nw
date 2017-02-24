// Implementation of the Needleman-Wunsch-algorithm and
// search routine for conserved elements
// adapted to g++

#include <stdio.h>
#include <time.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct conservedregion
{
  unsigned int startfirstspecies;
  unsigned int startsecondspecies;
  struct conservedregion *nextconservedregion; //Annika: pointer to next conserved region
};

//Annika: struct for scorearray:
struct scores
{
  int numberregions;
  struct conservedregion *firstconservedregion;
};

//Annika: catch the NWscore of an alignment
struct NWblock
{
  double NWmin;
  double NWmax;
};

unsigned int get_next_value(FILE *file) {
  char nextcharacter[2];
  unsigned int result;

  nextcharacter[1] = 0;
  nextcharacter[0] = getc(file);
  result = atoi(nextcharacter);
  nextcharacter[0] = getc(file);
  while (nextcharacter[0] != '\t') {
    result = 10*result;
    result = result+atoi(nextcharacter);
    nextcharacter[0] = getc(file);
  }
  return result;
}


int main(int argc, char* argv[]) {
 typedef long int dpfieldtype;

/*  const char *firstsequencefile = "/home/sascha/TempFilesPerl/firstsequencetoCcmdA.txt";
  const char *secondsequencefile = "/home/sascha/TempFilesPerl/secondsequencetoCcmdA.txt";
  const char *resultfile = "/home/sascha/TempFilesPerl/resultfromCcmdA.txt";
  const char *firstprofilefilename = "/home/sascha/TempFilesPerl/firstprofileA.txt";
  const char *secondprofilefilename = "/home/sascha/TempFilesPerl/secondprofileA.txt";
*/

  const int matchscore = 1;    // may not be changed
  const int mismatchscore = 0; // may not be changed
  
  // the next three values must fulfill: (effectivegapscore/effectivematchscore) = gapscore
  const int doublegapscore = -1;     // must be negative
  const dpfieldtype effectivematchscore = 2;
  const dpfieldtype effectivegapscore = -1;

  const unsigned int maxconsregionsoffset = 3000;    // number of conserved regions returned is at most
  const double maxconsregionsperbase = 0.4; 	     // maxconsregionsoffset + maxconsregionsperbase * #bases
  double *firstprofile;
  double *secondprofile;
  FILE *firstinputfile, *secondinputfile;
  FILE *outputfile, *outputfilefirstprofile, *outputfilesecondprofile;
  dpfieldtype *dpfield;
  char *a, *b;
  double cutoffthreshold;
  dpfieldtype charscore;
  double nextscoreboundary, largerstepwidth, lowerstepwidth;
  int offset1, offset2, i, j;
  int numberblocks1, numberblocks2, offsetcount1, offsetcount2;
  //Annika: integers for corridor and column in blockfield:
  int startRow, endpoint, distance, newcolumn;
							  dpfieldtype ScoreLeftPlusGap;
							  dpfieldtype ScoreAboveLeft;
							  dpfieldtype ScoreAboveLeftPlusMatch;
							  dpfieldtype * EntryBeingUpdated;
							  dpfieldtype ScoreAbovePlusGap;
  struct NWblock **blockfield;
  unsigned int firststepwidth, secondstepwidth, windowlength;
  unsigned int nonacgts, limitnonacgts, maxpairs;
  //Annika: for result data:
  int stackcounter, positioninscores;
  int lastpoint;
  struct scores *scorearray; // [(int)windowlengthPlus1]; 
  struct conservedregion *stackaddress;
  struct conservedregion* *pointers;

   if (argc < 6 ) {
	printf("Too few arguments\n");
	exit(13);
 }

 
 //~ for (i=0; i<argc; i++) {
	//~ printf("%s\n", argv[i]);
 //~ }
 
  const char *firstsequencefile = argv[1];
  const char *secondsequencefile = argv[2];
  const char *resultfile = argv[3];
  const char *firstprofilefilename = argv[4];
  const char *secondprofilefilename = argv[5];

 
  // read input from files
  firstinputfile = fopen(firstsequencefile,"r");
  if (!firstinputfile) {
    exit(13);
  }
  secondinputfile = fopen(secondsequencefile,"r");
  if (!secondinputfile) {
    exit(13);
  }
  firststepwidth = get_next_value(firstinputfile);
  secondstepwidth = get_next_value(firstinputfile);
  windowlength = get_next_value(firstinputfile);
  int windowlengthPlus1 = windowlength+1;

  cutoffthreshold = (double) get_next_value(firstinputfile);
  a = (char *) malloc(sizeof(char)*1000000);
  fgets(a,1000000,firstinputfile); // rest of first line
  fgets(a,1000000,firstinputfile);
  b = (char *) malloc(sizeof(char)*1000000);
  fgets(b,1000000,secondinputfile);
  a[strlen(a)-1] = 0; // remove newline-character
  b[strlen(b)-1] = 0; // remove newline-character
  fclose(firstinputfile);
  fclose(secondinputfile);  

  // do computation
  limitnonacgts = (unsigned int) (((double) windowlength)/((double) 4.0));
  numberblocks1 = (int) ((strlen(a))/firststepwidth);
  numberblocks1 -= (int) ((windowlength/firststepwidth)-1);
  while ((numberblocks1>0)&&
	 (((numberblocks1-1)*firststepwidth)+windowlength>strlen(a))) {
    numberblocks1--;
  };
  numberblocks2 = (int) ((strlen(b))/secondstepwidth);
  numberblocks2 -= (int) ((windowlength/secondstepwidth)-1);
  while ((numberblocks2>0)&&
	 (((numberblocks2-1)*secondstepwidth)+windowlength>strlen(b))) {
    numberblocks2--;
  };
  if (numberblocks1<0) {
    numberblocks1 = 0;
  }  
  if (numberblocks2<0) {
    numberblocks2 = 0;
  }
  firstprofile = (double *) malloc(sizeof(double)*numberblocks1);
  secondprofile = (double *) malloc(sizeof(double)*numberblocks2);
  for (i=0;i<numberblocks1;i++) {
    firstprofile[i] = 0.0;
  }
  for (i=0;i<numberblocks2;i++) {
    secondprofile[i] = 0.0;
  }
  //compute maxpairs
  maxpairs = maxconsregionsoffset + ((unsigned int) (maxconsregionsperbase * ((double) ((strlen(a)+strlen(b))/2))));

  //Annika: allocate memory for stack
  stackaddress=(struct conservedregion*)malloc(sizeof(struct conservedregion)*(int)maxpairs);
    
  //Annika: allocate memory for pointers:
  pointers = (struct conservedregion**)malloc(sizeof(struct conservedregion*)*(int)maxpairs);
  //Annika: initialize the pointers:
  for (i=0;i<maxpairs;i++) {
    pointers[i]=stackaddress+i;
  }
  //Annika: set stackcounter:
  stackcounter=0;

  //Annika: allocate memory for scorearray:
  scorearray = (struct scores*) malloc(sizeof(struct scores)*((int)windowlengthPlus1));
  //initialize scorearray.numberregions=0:
  for (i=0;i<=windowlength;i++) {
	  scorearray[i].numberregions=0;
  }    
  dpfield = (dpfieldtype *) malloc( windowlengthPlus1 * windowlengthPlus1 *sizeof(dpfieldtype));
  lastpoint =  windowlengthPlus1 *windowlength+windowlength;
  //Annika: allocate memory for blockfield:
  blockfield = (struct NWblock **) malloc(sizeof(struct NWblock*)*(numberblocks2));
  for (i=0;i<(int)numberblocks2;i++) {
	  blockfield[i] = (struct NWblock *) malloc(sizeof(struct NWblock)*(2)); 
  }

//  FILE * f = fopen("data.csv","w");

  for (offsetcount1 = 0;offsetcount1<numberblocks1;offsetcount1++) {
	  newcolumn = offsetcount1%2;//Annika: newcolumn for blockfield; oldcolumn=(newcolumn+1)%2;
	  offset1 = offsetcount1*firststepwidth;
	  nonacgts = 0;
//	  fprintf(f,"%100s",a +offset1);
	  for (i=0;i< (int) windowlength;i++) {
		  switch(a[offset1+i])
		  {
		  case 'A':
		  case 'C':
		  case 'G':
		  case 'T':
			  break;
		  default:
			  a[offset1+i] = '@';
			  nonacgts++;
		  }
	  }
//	  fprintf(f,"\n");
	  if (nonacgts<limitnonacgts) {
		  for (offsetcount2 = 0;offsetcount2<numberblocks2;offsetcount2++) {
			  offset2 = offsetcount2*secondstepwidth;
			  nonacgts = 0;      
//			  fprintf(f,"%100s",b +offset2);
			  for (i=0;i< (int) windowlength;i++) {
				  switch(b[offset2+i])
				  {
				  case 'A':
				  case 'C':
				  case 'G':
				  case 'T':
					  break;
				  default:
					  b[offset2+i] = '#';
					  nonacgts++;
				  }
			  }
//			  fprintf(f,"\n");
//			  fclose(f);
			  struct NWblock * boff2NewCol = &blockfield[offsetcount2][newcolumn];

			  if (nonacgts<limitnonacgts) {
				  //Annika: decision to compute the Needleman-Wunsch block, = NWblock or not; 
				  //fill NWmax, NWmin
				  if (offsetcount1>0) {
					  if (offsetcount2>0) {

						  //	
						  boff2NewCol->NWmax = (blockfield[offsetcount2][(newcolumn+1)%2].NWmax+((double)firststepwidth))-((double)firststepwidth)*doublegapscore;
						  boff2NewCol->NWmin = (blockfield[offsetcount2][(newcolumn+1)%2].NWmin-((double)firststepwidth))+((double)firststepwidth)*doublegapscore;

						  if (((blockfield[offsetcount2-1][newcolumn].NWmax+((double)secondstepwidth))-((double)secondstepwidth)*doublegapscore)<boff2NewCol->NWmax) {			  
							  boff2NewCol->NWmax = (blockfield[offsetcount2-1][newcolumn].NWmax+((double)secondstepwidth))-((double)secondstepwidth)*doublegapscore;
						  }

						  if (((blockfield[offsetcount2-1][newcolumn].NWmin-((double)secondstepwidth))+((double)secondstepwidth)*doublegapscore)>boff2NewCol->NWmin) {			   
							  boff2NewCol->NWmin = (blockfield[offsetcount2-1][newcolumn].NWmin-((double)secondstepwidth))+((double)secondstepwidth)*doublegapscore;
						  }
						  largerstepwidth = ((double)firststepwidth);
						  lowerstepwidth = ((double)secondstepwidth);
						  if (lowerstepwidth>largerstepwidth) {
							  largerstepwidth = ((double)secondstepwidth);
							  lowerstepwidth = ((double)firststepwidth);
						  }
						  nextscoreboundary = (blockfield[offsetcount2-1][(newcolumn+1)%2].NWmax+largerstepwidth)-(largerstepwidth-lowerstepwidth)*doublegapscore;
						  if (nextscoreboundary<boff2NewCol->NWmax) {
							  boff2NewCol->NWmax = nextscoreboundary;
						  }
						  nextscoreboundary = (blockfield[offsetcount2-1][(newcolumn+1)%2].NWmin-largerstepwidth)+(largerstepwidth-lowerstepwidth)*doublegapscore;
						  if (nextscoreboundary>boff2NewCol->NWmin) {
							  boff2NewCol->NWmin = nextscoreboundary;
						  }
						  if (offsetcount2+1<numberblocks2) {
							  nextscoreboundary = (blockfield[offsetcount2+1][(newcolumn+1)%2].NWmax+largerstepwidth+lowerstepwidth)+(-(largerstepwidth+lowerstepwidth)*doublegapscore);
							  if (nextscoreboundary<boff2NewCol->NWmax) {
								  boff2NewCol->NWmax = nextscoreboundary;
							  }
							  nextscoreboundary = (blockfield[offsetcount2+1][(newcolumn+1)%2].NWmin-(largerstepwidth+lowerstepwidth))+((largerstepwidth+lowerstepwidth)*doublegapscore);
							  if (nextscoreboundary>boff2NewCol->NWmin) {
								  boff2NewCol->NWmin = nextscoreboundary;
							  }
						  }
					  }
					  else {
						  boff2NewCol->NWmax = (blockfield[offsetcount2][(newcolumn+1)%2].NWmax+((double)firststepwidth))-((double)firststepwidth)*doublegapscore;
						  boff2NewCol->NWmin = (blockfield[offsetcount2][(newcolumn+1)%2].NWmin-((double)firststepwidth))+((double)firststepwidth)*doublegapscore;
					  }
				  }
				  else { //offsetcount1==0
					  if (offsetcount2>0) {
						  boff2NewCol->NWmax = (blockfield[offsetcount2-1][newcolumn].NWmax+((double)secondstepwidth))-((double)secondstepwidth)*doublegapscore;
						  boff2NewCol->NWmin = (blockfield[offsetcount2-1][newcolumn].NWmin-((double)secondstepwidth))+((double)secondstepwidth)*doublegapscore;
					  }
					  else {
						  boff2NewCol->NWmax = ((double)windowlength);
						  boff2NewCol->NWmin = 0;
					  }
				  }
				  if (boff2NewCol->NWmax>((double)windowlength)) {
					  boff2NewCol->NWmax = ((double)windowlength);
				  }
				  if (boff2NewCol->NWmin<0) {
					  boff2NewCol->NWmin = 0.0;
				  }		  
				  if (( boff2NewCol->NWmax>cutoffthreshold)||( boff2NewCol->NWmax>firstprofile[offsetcount1])||( boff2NewCol->NWmax>secondprofile[offsetcount2])) 
				  {
					  //Annika: calculating the distance for the corridor:
					  distance = ((int)(((double)windowlength)-boff2NewCol->NWmin)/(1.0-doublegapscore));
					  distance++;
					  if (distance>windowlength) {
						  distance = windowlength;
					  }
					  if (distance<1) {
						  distance = 1;
					  }		      
					  for (j=-1;j<=distance-1;j++) {
						  dpfield[j+1] = effectivegapscore*((dpfieldtype) (j+1));
					  }
//					  		    			  MyTimer t;

					  int row;
					  for (row=0;row<((int)windowlength);row++) {
						  startRow = row-distance;
						  int nextRow = row+1;
						  if (startRow<-1) {
							  startRow = -1;
						  }
						  int startPointPlusOne = startRow+1;
						
						  if ((row+distance) < (int)windowlength)
							  endpoint = row+distance;
						  else
							  endpoint = (int)windowlength;

						  if (startRow==-1) {
							  dpfield[ nextRow * windowlengthPlus1 ] = effectivegapscore*((dpfieldtype)  nextRow );
						  }
						  else {
							  if (a[offset1+row] == b[offset2+startRow]) {
								  charscore = effectivematchscore;
							  }
							  else {
								  charscore = (dpfieldtype) mismatchscore;
							  }
							  dpfieldtype a = dpfield[row* windowlengthPlus1 +startRow]+charscore;
							  dpfieldtype b = dpfield[row* windowlengthPlus1 + startPointPlusOne ]+effectivegapscore;
							  if (a > b)
								  dpfield[ nextRow * windowlengthPlus1 + startPointPlusOne ] = a;
							  else
								  dpfield[ nextRow * windowlengthPlus1 + startPointPlusOne ] = b;

//							  dpfield[ nextRow * windowlengthPlus1 + startPointPlusOne ] = dpfield[row* windowlengthPlus1 +startRow]+charscore;
//							  if (dpfield[row* windowlengthPlus1 + startPointPlusOne ]+effectivegapscore>dpfield[ nextRow * windowlengthPlus1 + startPointPlusOne ]) {
//								  dpfield[ nextRow * windowlengthPlus1 + startPointPlusOne ] = dpfield[row* windowlengthPlus1 + startPointPlusOne ]+effectivegapscore;
//							  }
						  }

						  for (j=startPointPlusOne;j<endpoint;j++) {
							  ScoreAboveLeft = dpfield[row* windowlengthPlus1 +j];
							  ScoreAbovePlusGap = dpfield[row* windowlengthPlus1 +j +1] + effectivegapscore;
							  ScoreLeftPlusGap = dpfield[ nextRow * windowlengthPlus1  +j ] + effectivegapscore;
							  ScoreAboveLeftPlusMatch = ScoreAboveLeft + effectivematchscore;
							  EntryBeingUpdated = &dpfield[ nextRow * windowlengthPlus1  +j +1];;
							  if (a[offset1+row]==b[offset2+j]) {

								  if (ScoreAbovePlusGap>ScoreAboveLeftPlusMatch) {
									  if (ScoreLeftPlusGap>ScoreAbovePlusGap) {
										  *EntryBeingUpdated = ScoreLeftPlusGap;
									  }
									  else {
										  *EntryBeingUpdated = ScoreAbovePlusGap;
									  }
								  }
								  else {
									  if (ScoreLeftPlusGap > ScoreAboveLeftPlusMatch)
									  {
										  *EntryBeingUpdated = ScoreLeftPlusGap;
									  }
									  else {
										  *EntryBeingUpdated = ScoreAboveLeftPlusMatch;
									  }
								  }
							  }
							  else {
								  if (ScoreAbovePlusGap>ScoreAboveLeft) {
									  if (ScoreLeftPlusGap>ScoreAbovePlusGap) {
										  *EntryBeingUpdated = ScoreLeftPlusGap;
									  }
									  else {
										  *EntryBeingUpdated = ScoreAbovePlusGap;
									  }
								  }
								  else {
									  if (ScoreLeftPlusGap>ScoreAboveLeft) {
										  *EntryBeingUpdated = ScoreLeftPlusGap;
									  }
									  else {
										  *EntryBeingUpdated = ScoreAboveLeft;
									  }
								  }
							  }
						  }
						  if (endpoint<((int)windowlength)) {
							  if (a[offset1+row]==b[offset2+endpoint]) {
								  charscore = effectivematchscore;
							  }
							  else {
								  charscore = (dpfieldtype) mismatchscore;
							  }
							  int a = dpfield[row* windowlengthPlus1 +endpoint]+charscore;
							  int b = dpfield[ nextRow * windowlengthPlus1 +endpoint]+effectivegapscore;
							  if (a > b)
								  dpfield[ nextRow * windowlengthPlus1 +(endpoint+1)] = a;
							  else
								  dpfield[ nextRow * windowlengthPlus1 +(endpoint+1)] = b;


//							  dpfield[ nextRow * windowlengthPlus1 +(endpoint+1)] = dpfield[row* windowlengthPlus1 +endpoint]+charscore;
//							  if (dpfield[ nextRow * windowlengthPlus1 +endpoint]+effectivegapscore>dpfield[ nextRow * windowlengthPlus1 +(endpoint+1)]) {
//								  dpfield[ nextRow * windowlengthPlus1 +(endpoint+1)] = dpfield[ nextRow * windowlengthPlus1 +endpoint]+effectivegapscore;
//							  }
						  }
					  }
					  dpfield[lastpoint] = dpfield[lastpoint]/effectivematchscore;

					  //Annika: fill blockfield:
					  boff2NewCol->NWmax = (double) dpfield[lastpoint];
					  boff2NewCol->NWmin = (double) dpfield[lastpoint];

					  if (dpfield[lastpoint]>firstprofile[offsetcount1]) {
						  firstprofile[offsetcount1] = dpfield[lastpoint];
					  }
					  if (dpfield[lastpoint]>secondprofile[offsetcount2]) {
						  secondprofile[offsetcount2] = dpfield[lastpoint];
					  }
					  if ((double) dpfield[lastpoint] > cutoffthreshold) {
						  //Annika: save results:
						  positioninscores = ((int)ceil((double) dpfield[lastpoint]));
						  if (scorearray[positioninscores].numberregions>0) {
							  pointers[stackcounter]->nextconservedregion=scorearray[positioninscores].firstconservedregion;
						  }
						  else {
							  pointers[stackcounter]->nextconservedregion=NULL;
						  }
						  scorearray[positioninscores].numberregions++;
						  scorearray[positioninscores].firstconservedregion=pointers[stackcounter];
						  scorearray[positioninscores].firstconservedregion->startfirstspecies=(unsigned int)offset1;
						  scorearray[positioninscores].firstconservedregion->startsecondspecies=(unsigned int)offset2;

						  if (stackcounter<(maxpairs-1)) {
							  stackcounter++;
						  }
						  else { //last pointer used, therefore increase cutoffthreshold
							  cutoffthreshold = cutoffthreshold+1.0;
							  while (scorearray[(int)cutoffthreshold].numberregions==0) {
								  cutoffthreshold = cutoffthreshold + 1.0;
							  }
							  //Annika: put freed pointers to conservedregions into "pointers":
							  stackcounter++;
							  for (row=1;row<=scorearray[(int)cutoffthreshold].numberregions;row++) {
								  pointers[(int)maxpairs-row]=scorearray[(int)cutoffthreshold].firstconservedregion;
								  scorearray[(int)cutoffthreshold].firstconservedregion=scorearray[(int)cutoffthreshold].firstconservedregion->nextconservedregion;
								  stackcounter--;
							  }
						  }//increase cutoffthreshold
					  }
				  }//if compute NWblock		 
			  }	      
			  else {
				  boff2NewCol->NWmin = 0.0;
				  boff2NewCol->NWmax = ((double)windowlength);
			  }
		  }
	  }
	  else {
		  for (i=0;i<numberblocks2;i++) {
			  blockfield[i][newcolumn].NWmin = 0.0;
			  blockfield[i][newcolumn].NWmax = ((double)windowlength);
		  }
	  }
//	  TimeReport();
  }
  free(dpfield);

  // write result to file
  outputfile = fopen(resultfile,"w");
  outputfilefirstprofile = fopen(firstprofilefilename,"w");
  outputfilesecondprofile = fopen(secondprofilefilename,"w");
  if ((!outputfile)||(!outputfilefirstprofile)||(!outputfilesecondprofile)) {
	  exit(13);
  }
  fprintf(outputfile,"%d\n",(int) cutoffthreshold);

  //Annika: output results from scorearray:
  for (i=1+(int)cutoffthreshold;i<=(int)windowlength;i++) {
	  if (scorearray[i].numberregions!=0) {
		  for (j=1;j<=(scorearray[i].numberregions);j++) {
			  fprintf(outputfile,"%d\t",scorearray[i].firstconservedregion->startfirstspecies);
			  fprintf(outputfile,"%d\t",scorearray[i].firstconservedregion->startsecondspecies);
			  fprintf(outputfile,"%d\t",i);
			  fprintf(outputfile,"\n");
			  scorearray[i].firstconservedregion=scorearray[i].firstconservedregion->nextconservedregion;
		  }//for j
	  }
  }//for i

  for (i=0;i<numberblocks1;i++) {
//	  fprintf(outputfilefirstprofile,"%f\n",firstprofile[i]);
	  fprintf(outputfilefirstprofile,"%f\n",firstprofile[i]);
  }
  for (i=0;i<numberblocks2;i++) {
	  fprintf(outputfilesecondprofile,"%f\n",secondprofile[i]);
  }
  fclose(outputfile);
  fclose(outputfilefirstprofile);
  fclose(outputfilesecondprofile);

  //Annika: free blockfield:
  for (i=0;i<(int)numberblocks2;i++) {
	  free (blockfield[i]);
  }
  free (blockfield);

  //Annika: free conservedregions in scorearray and stack:
  free (stackaddress);

  //Annika: free scorearray:
  free (scorearray);

  //Annika: free pointers:
  free (pointers);

  free(firstprofile);
  free(secondprofile);
  free(a);
  free(b);

  return 0;
} // main //