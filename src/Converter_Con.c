#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define MATCH 2
#define MISMATCH -2
#define GAP -1
#define VERBOSE 1
#define SIMILARITY 0.4
#define MAXLEN 1000

#define PPVTHRESH 0

#define max(A,B) (((A)>(B))?(A):(B))
#define score(A,B) (((A)==(B))?(MATCH):(MISMATCH))

int A[2*MAXLEN][2*MAXLEN];
int Contact[MAXLEN][MAXLEN];
int n,m;

void fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

double global_align(char Seq1[MAXLEN],char Seq2[MAXLEN])
{
		int i,j,max=0;

		for(i=0;i<2*MAXLEN;i++)	
			for(j=0;j<2*MAXLEN;j++)	
				A[i][j]=0;

		m=strlen(Seq1);
		n=strlen(Seq2);
		for(i=1;i<=m;i++)
			for(j=1;j<=n;j++)
				A[i][j] = max( A[i-1][j]+GAP , max(A[i][j-1]+ GAP , max( A[i-1][j-1]+score(Seq1[i-1],Seq2[j-1]),0)));

		for(i=1;i<strlen(Seq1);i++)
			for(j=1;j<strlen(Seq2);j++)
			{
				if( A[i][j] > max )
				{
					m =i;
					n =j;
					max = A[m][n];
				}
			}
		return A[m][n];
}

void appendchar(char c,char String[MAXLEN])
{
		char aux[MAXLEN];
		strcpy(aux,String);
		sprintf(String,"%c%s",c,aux);
}

int main(int argc, char *argv[])
{
	int sequence_found=0;
	int i,j, Score, ScoreDiag, ScoreUp, ScoreLeft,size;
	int fasta2pdb[MAXLEN], k;
	int total_contacts=0,true_positives,false_positives,res1,res2,aux1,aux2;
	char AlignmentA[2*MAXLEN],AlignmentB[2*MAXLEN];
	char Fasta_Seq[MAXLEN],AUX[81], holder[2*MAXLEN];
	char ALN_Seq[MAXLEN], ALN_Seq_Cand[MAXLEN];
	double dir_inf, aux, candScore = 0, highScore = 0;
	FILE *input_fasta,*input_aln, *input_contacts, *input_map;

	true_positives=0;
	false_positives=0;

	if(argc<5)
	{
		printf("Usage: %s fasta_file aln_file psicov_file contact_map_file\n",argv[0]);
		return 0;
	}

	input_fasta    = fopen(argv[1],"r");
	if (! input_fasta)
            fail("Unable to open fasta file!");

	input_aln      = fopen(argv[2],"r");
	if (! input_aln)
            fail("Unable to open aln file!");

	input_contacts = fopen(argv[3],"r");
	if (! input_contacts)
            fail("Unable to open psicov file!");

	input_map      = fopen(argv[4],"r");
	if (! input_map)
            fail("Unable to open contact_map file!");
	

	/* THROW AWAY THE FIRST LINE (HEADER) OF FASTA FILE   */
	fscanf(input_fasta,"%s",Fasta_Seq);
	/* READ INPUT FASTA SEQUENCE 	  					   */
	fscanf(input_fasta,"%s",Fasta_Seq);
	while(fscanf(input_fasta,"%s",AUX)!=EOF) 	
		strcat(Fasta_Seq,AUX);

	highScore = 10;
	/* FIND A SIMILAR SEQUENCE IN THE ALIGNMENT	 */
	for(i=1;fscanf(input_aln,"%s",ALN_Seq_Cand)!=EOF;i++)
	{
		/* If the score of the global alignment is significantly high */	
		candScore = global_align(Fasta_Seq,ALN_Seq_Cand);
	//	printf("%s\n%s\n%.2lf\n",Fasta_Seq,ALN_Seq_Cand,candScore);
		if( candScore > highScore )		
		{
			sequence_found=1;
			highScore = candScore;
			strcpy(ALN_Seq, ALN_Seq_Cand);
			if (VERBOSE)
				fprintf(stderr,"Similar sequence found in line %d\n\n",i);
            break;
		}
	}
/*
	if(!sequence_found)
	{
		fprintf(stderr,"Could not find a sequence similar to %s in the alignment file %s\n",argv[1],argv[2]);
		fclose(input_map);
		fclose(input_fasta);
		fclose(input_aln);
		fclose(input_contacts);
		return 0;
	}*/
	
    // Have to redo the winning alignment
	global_align(Fasta_Seq,ALN_Seq); 
	/* Print the two sequences */
	if (VERBOSE)
		fprintf(stderr,"%s\n%s\n\n",Fasta_Seq,ALN_Seq);

	/***** PERFORM THE TRACEBACK STEP OF THE GLOBAL ALIGNMENT *****/
	// The unconventional char '&' is used to distinguish gaps that are part of pfam from gaps that are introduced in the above alignment
	strcpy(AlignmentA,"\0");	
	strcpy(AlignmentB,"\0");
	i = m;
	j = n;

	while (/*i > 0 && j > 0*/ A[i][j] != 0)
	{
		Score     = A[i][j];
    		ScoreDiag = A[i - 1][j - 1];
	    	ScoreUp   = A[i][j - 1];
    		ScoreLeft = A[i - 1][j];

		if (Score == ScoreUp + GAP)
    		{
	    		appendchar('&', AlignmentA);
	 	    	appendchar(ALN_Seq[j-1], AlignmentB);
	 	 	j = j - 1;
	    	}
		else if (Score == ScoreLeft + GAP)
    		{
			appendchar(Fasta_Seq[i-1], AlignmentA);
   			appendchar('&', AlignmentB);
    	    		i = i - 1;
	    	}
    		else if ( Score ==  ScoreDiag  +  score(Fasta_Seq[i-1], ALN_Seq[j-1]) )
	    	{
		    appendchar(Fasta_Seq[i-1], AlignmentA);
		    appendchar(ALN_Seq[j-1], AlignmentB);
		    i = i - 1;
		    j = j - 1;
	    	}
	}
	while (i > 0)
	{
		appendchar(Fasta_Seq[i-1], AlignmentA);
   		appendchar('&', AlignmentB);
	    	i = i - 1;
	}
	while (j > 0)
	{
	    appendchar('&', AlignmentA);
	 	appendchar(ALN_Seq[j-1], AlignmentB);
    		j = j - 1;
	}

	
	/* Print the aligned sequences! */
	if (VERBOSE)
		fprintf(stderr,"### SEQUENCE ALIGNMENT: ###\n\n%s\n%s\n\n",AlignmentA,AlignmentB);
	
	/***** DETERMINE THE OFFSET FOR EACH POSITION *****/
    /* Replace offset with a fasta2pdb vector such that fasta[0] contains the pdb index of the first fasta residue
     * Each entry is either a number >0 or 0 (ie no match). */ 

	for(i=1,j=1,k=0; k<strlen(AlignmentA); k++)
	{
		// If struc is gapped, but aln isn't
		if(AlignmentA[k]!='&' && AlignmentB[k]=='&')
		{
			j++;
		}
		// If pfam is gapped, but struc isn't
		if(AlignmentA[k]=='&' && AlignmentB[k]!='&')
		{
			fasta2pdb[i]=0;
			i++;
		}
		if(AlignmentA[k]!='&' && AlignmentB[k]!='&')
		{
			fasta2pdb[i]=j;
			i++;
			j++;
		}
	}
	for(j=1;j<i;j++)
		fprintf(stderr,"%d ",fasta2pdb[j]);
	fprintf(stderr,"\n");	

	/* Print the Offset array */
    if (VERBOSE)
    {
        fprintf(stderr,"\n\n%s\n",ALN_Seq);  
    	for(i=1;strlen(ALN_Seq);i++)
	    {
		    fprintf(stderr,"%d %d\n",i,fasta2pdb[i]);
    	}
	    printf("\n");
    }

	/***** LOAD THE TRUE CONTACT MAP *****/
    fprintf(stderr, "Now reading in contacts\n");
    fscanf(input_map, "%d", &size);
    fprintf(stderr, "Size is %d\n", size);
    while (fscanf(input_map,"%d %d",&i,&j)!=EOF)
    {
        //fprintf(stderr, "%d %d\n", i, j);
        Contact[i][j] = 1;
        Contact[j][i] = 1;
	    total_contacts++;
    }
	fprintf(stderr,"Contact map was read correctly.\n");

    /* Make sure contacts are not counted twice! */
	total_contacts/=2; 

	/***** CORRECT THE POSITIONS OF THE PREDICTED CONTACTS *****/
	while(fscanf(input_contacts,"%d %d %d %lf %lf",&i,&j,&aux2,&aux,&aux) != EOF)
	{
        /* Weird indexing for predicted contacts, skip pair */
		if(i >= strlen(AlignmentA) || j >= strlen(AlignmentA) || i>MAXLEN-1 || j > MAXLEN-1 )
			continue;

        /* If contact score is too low, stop reading in predicted contacts */        
        if(aux<0.5) break;

		//fprintf(stderr,"Predicted: i=%d j=%d\n",i,j);
		res1 = fasta2pdb[i];
		res2 = fasta2pdb[j];

        if(VERBOSE)
        {
            fprintf(stderr,"Adjusted: res1=%d res2=%d\n",res1,res2);
        	fprintf(stderr,"Predicted: i=%d j=%d\n",i,j);	
        }

        /* Adjusted contacts not present in the PDB! */
		if(res1 <= 0 || res2 <= 0 || res1> MAXLEN-1 || res2 > MAXLEN-1) 
                printf("%d %d %lf -1\n",res1,res2,aux);
        else
        { 
    		printf("%d %d %lf %d\n",res1,res2,aux,Contact[res1-1][res2-1]);
    		if(Contact[res1-1][res2-1] && aux >= 0.5 )	true_positives++;
            else false_positives++;
        }
	}

	/* Print total accuracy */
	if (VERBOSE)
	{
		fprintf(stderr,"True Positives: %.2lf %d %d\n",100.0*true_positives/(double)(true_positives+false_positives),true_positives,true_positives+false_positives);
		fprintf(stderr,"False Positives: %.2lf %d %d\n",100.0*false_positives/(double)(true_positives+false_positives),false_positives,true_positives+false_positives);
		fprintf(stderr,"Total protein contacts: %d\n",total_contacts);
	}

	fclose(input_map);
	fclose(input_fasta);
	fclose(input_aln);
	fclose(input_contacts);
	return 0;
}
