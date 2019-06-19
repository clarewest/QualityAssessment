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
#define min(A,B) (((A)<(B))?(A):(B))

#define score(A,B) (((A)==(B))?(MATCH):(MISMATCH))

int A[2000][2000];
int Contact[2000][2000];
int n,m;

void fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

double global_align(char Seq1[2000],char Seq2[2000])
{
		int i,j,max=0;

		for(i=0;i<2000;i++)	
			for(j=0;j<2000;j++)	
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

void appendchar(char c,char String[2000])
{
		char aux[2000];
		strcpy(aux,String);
		sprintf(String,"%c%s",c,aux);
}

void openfile(FILE **input,char *argv[], char suffix[100])
{
    char AUX[200],AUX2[200];
    strcpy(AUX,argv[1]);
    strcat(AUX,suffix);
    *input  = fopen(AUX,"r");
    sprintf(AUX2,"Unable to open %s file!",AUX);
    if (!*input)
        fail(AUX2);
}

int main(int argc, char *argv[])
{
	int sequence_found=0,status;
	int i,j, Score, ScoreDiag, ScoreUp, ScoreLeft,size;
	int pfam2pdb[2000], k,l;
	int total_contacts=0,true_positives,false_positives,res1,res2,aux1,aux2, true_positives_short,true_positives_long,false_positives_short,false_positives_long;
    int true_positives_h,false_positives_h,true_positives_short_h,true_positives_long_h,false_positives_short_h,false_positives_long_h;

    int phi,phi_l,psi,psi_l;
	char AlignmentA[2000],AlignmentB[2000],AlignmentC[2000],AlignmentD[2000];
    char c;
	char Fasta_Seq[2000],AUX[181],Proxy_Seq[2000],Proxy_SS[2000], holder[MAXLEN], SS[2000];
	char ALN_Seq[2000], ALN_Seq_Cand[2000];
	double dir_inf, aux, candScore = 0, highScore = 0, TA[2][2000],Pred_TA[2][2000],MAE=0.0;
    double diff_phi,diff_phi_l,diff_psi,diff_psi_l;
	FILE *input_fasta,*input_proxy,*input_dssp,*input_psipred,*input_angles,*input_pred_angles,*input_aln, *input_contacts, *input_map;
    
    diff_phi = 0.0;
    diff_phi_l = 0.0;
    diff_psi = 0.0;
    diff_psi_l = 0.0;
    phi=0;
    phi_l=0;
    psi=0;
    psi_l=0;
	true_positives=0;
	false_positives=0;

	if(argc!=2)
	{
		printf("Usage: %s PDB_ID\n",argv[0]);
		return 0;
	}

    openfile(&input_fasta,argv,".fasta.txt");
    /* THROW AWAY THE FIRST LINE OF FASTA FILE   */
    for(c=fgetc(input_fasta); c!='\n'; c=fgetc(input_fasta));
    /* READ INPUT FASTA SEQUENCE                           */
    fscanf(input_fasta,"%s",Fasta_Seq);
    while(fscanf(input_fasta,"%s",AUX)!=EOF)
        strcat(Fasta_Seq,AUX);

    openfile(&input_proxy,argv,".proxy_fasta");
    /* THROW AWAY THE FIRST LINE OF PROXY_FASTA FILE   */
    for(c=fgetc(input_proxy); c!='\n'; c=fgetc(input_proxy) );
    /* READ INPUT PROXY_FASTA SEQUENCE                 */
    fscanf(input_proxy,"%s",Proxy_Seq);
    
    /* READ TRUE SS - SAME LENGTH AS PROXY SEQ */
    openfile(&input_dssp,argv,".dssp_ss");
    fscanf(input_dssp,"%s",Proxy_SS);

    /* READ PREDICTED SS - SAME LENGTH AS FASTA SEQ */
    openfile(&input_psipred,argv,".psipred_ss");
    fscanf(input_psipred,"%s",SS);

    openfile(&input_angles,argv,".angles");
    /* THROW AWAY THE FIRST LINE OF ANGLE FILE   */
    for(c=fgetc(input_angles); c!='\n'; c=fgetc(input_angles) );
    for(i=0; i< strlen(Proxy_Seq);i++)
        fscanf(input_angles,"%d %c %lf %lf",&j,&c,&TA[0][i],&TA[1][i]);

    openfile(&input_pred_angles,argv,".pred_angles");
    for(i=0; i< strlen(Fasta_Seq);i++)
        fscanf(input_pred_angles,"%lf %lf",&Pred_TA[0][i],&Pred_TA[1][i]);

    openfile(&input_aln,argv,".aln");    
    openfile(&input_contacts,argv,".metapsicov.stage1");
	openfile(&input_map,argv,".proxy_map");
	
	/* ALING THE TWO SEQUENCES */
	global_align(Fasta_Seq,Proxy_Seq);

	/***** PERFORM THE TRACEBACK STEP OF THE GLOBAL ALIGNMENT *****/
	strcpy(AlignmentA,"\0");	
	strcpy(AlignmentB,"\0");
    strcpy(AlignmentC,"\0");
    strcpy(AlignmentD,"\0");

	i = m;
	j = n;

	while (A[i][j] != 0 && i>0 && j>0)
	{
    		Score     = A[i][j];
    		ScoreDiag = A[i - 1][j - 1];
	    	ScoreUp   = A[i][j - 1];
    		ScoreLeft = A[i - 1][j];
            if (Score == ScoreUp + GAP)
    		{
	    		appendchar('&', AlignmentA);
                appendchar('&', AlignmentC);
	 	    	appendchar(Proxy_Seq[j-1], AlignmentB);
                appendchar(Proxy_SS[j-1], AlignmentD);
	 	 	    j = j - 1;
	    	}
		    else
            {
                if (Score == ScoreLeft + GAP)
    		    {
    			    appendchar(Fasta_Seq[i-1], AlignmentA);
                    appendchar(SS[i-1], AlignmentC);
   	        		appendchar('&', AlignmentB);
                    appendchar('&', AlignmentD);
    	   	    	i = i - 1;
	    	    }
            	else
                {
                    if ( Score ==  ScoreDiag  +  score(Fasta_Seq[i-1], Proxy_Seq[j-1]) )
        	    	{
    	        	    appendchar(Fasta_Seq[i-1], AlignmentA);
                        appendchar(SS[i-1], AlignmentC);
        	    	    appendchar(Proxy_Seq[j-1], AlignmentB);
                        appendchar(Proxy_SS[j-1], AlignmentD);
                        if(TA[0][j-1]>-370 && Pred_TA[0][i-1]>-370)
                        {
                            diff_phi += min(fabs(TA[0][j-1]-Pred_TA[0][i-1]), fabs(fabs(TA[0][j-1])+fabs(Pred_TA[0][i-1]) - 360));
                            phi++;
                            if(Proxy_SS[j-1]!='H'&&Proxy_SS[j-1]!='E')
                            {   
                                diff_phi_l += min(fabs(TA[0][j-1]-Pred_TA[0][i-1]), fabs(fabs(TA[0][j-1])+fabs(Pred_TA[0][i-1]) - 360));
                                phi_l++;
                            }
                        }

                        if(TA[1][j-1]>-370 && Pred_TA[1][i-1]>-370)
                        {   
                            diff_psi += min(fabs(TA[1][j-1]-Pred_TA[1][i-1]), fabs(fabs(TA[1][j-1])+fabs(Pred_TA[1][i-1]) - 360));
                            psi++;
                            if(Proxy_SS[j-1]!='H'&&Proxy_SS[j-1]!='E')
                            {   
                                diff_psi_l += min(fabs(TA[1][j-1]-Pred_TA[1][i-1]), fabs(fabs(TA[1][j-1])+fabs(Pred_TA[1][i-1]) - 360));
                                psi_l++;
                            }   
                        }


    		            i = i - 1;
        	    	    j = j - 1;
                    }
                }
	    	}
	}
	while (i > 0)
	{
		appendchar(Fasta_Seq[i-1], AlignmentA);
        appendchar(SS[i-1], AlignmentC);
   		appendchar('&', AlignmentB);
        appendchar('&', AlignmentD);
       	i = i - 1;
	}
	while (j > 0)
	{
	    appendchar('&', AlignmentA);
        appendchar('&', AlignmentC);
	 	appendchar(Proxy_Seq[j-1], AlignmentB);
        appendchar(Proxy_SS[j-1], AlignmentD);
   		j = j - 1;
	}

	
	/* Print the aligned sequences! */
//	if (VERBOSE)
//		fprintf(stderr,"### SEQUENCE ALIGNMENT: ###\n\n%s\n%s\n%s\n%s\n\n",AlignmentA,AlignmentB,AlignmentC,AlignmentD);
    

    for(i=0,j=0,k=0,l=0; i<strlen(AlignmentA);i++)
    {
        if(AlignmentC[i]!='&'&& AlignmentD[i]!='&')
        {
            j++;
            if (AlignmentC[i]==AlignmentD[i])
            {
                k++;
                l++;
            }
            else
            {
                if (AlignmentC[i]!='E'&&AlignmentC[i]!='H'&&AlignmentD[i]!='E'&&AlignmentD[i]!='H')
                    l++;
            }
        }
    }
    printf("SS Precision = %d% (%d/%d)\n",k*100/j,k,j);
    printf("SS Precision (3-state) = %d% (%d/%d)\n",l*100/j,l,j);

    printf("MAE Phi = %.2f \t MAE Psi = %.2f\n",diff_phi/phi,diff_psi/psi);
    printf("MAE Phi (loops) = %.2f \t MAE Psi (loops) = %.2f\n",diff_phi_l/phi_l,diff_psi_l/psi_l);

    highScore = 10;
    /* FIND A SIMILAR SEQUENCE IN THE ALIGNMENT  */
    for(i=1;fscanf(input_aln,"%s",ALN_Seq_Cand)!=EOF;i++)
    {
        /* If the score of the global alignment is significantly high */
        candScore = global_align(Fasta_Seq,ALN_Seq_Cand);
        if( candScore > highScore )
        {
            sequence_found=1;
            highScore = candScore;
            strcpy(ALN_Seq, ALN_Seq_Cand);
            break;
        }
    }
    global_align(Fasta_Seq,ALN_Seq);

    /***** PERFORM THE TRACEBACK STEP OF THE GLOBAL ALIGNMENT *****/
    strcpy(AlignmentA,"\0");
    strcpy(AlignmentB,"\0");

    i = m;
    j = n;

    while (A[i][j] != 0)
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

	for(i=1,j=1,k=0; k<strlen(AlignmentA); k++)
	{
		// If struc is gapped, but pfam isn't
		if(AlignmentA[k]!='&' && AlignmentB[k]=='&')
		{
			j++;
		}
		// If pfam is gapped, but struc isn't
		if(AlignmentA[k]=='&' && AlignmentB[k]!='&')
		{
			pfam2pdb[i]=0;
			i++;
		}
		if(AlignmentA[k]!='&' && AlignmentB[k]!='&')
		{
			pfam2pdb[i]=j;
			i++;
			j++;
		}
	}

    fscanf(input_map, "%d", &size);
    while (fscanf(input_map,"%d %d",&i,&j)!=EOF)
    {
           Contact[i][j] = 1;
           Contact[j][i] = 1;
	       total_contacts++;
    }

	total_contacts/=2; 
	/***** CORRECT THE POSITIONS OF THE PREDICTED CONTACTS *****/
	/* **** PSICOV **** */
	while(fscanf(input_contacts,"%d %d %d %lf %lf",&i,&j,&aux2,&aux,&aux) != EOF)
	{
		if(i >= strlen(AlignmentA) || j >= strlen(AlignmentA) || i>1999 || j > 1999 )
			continue;
		res1 = pfam2pdb[i];
		res2 = pfam2pdb[j];
		
		if(aux<0.5) break;


		if(res1 <= 0 || res2 <= 0 || res1> 1999 || res2 > 1999) continue;

		//printf("%d %d %lf %d\n",res1,res2,aux,Contact[res1-1][res2-1]);

		if(Contact[res1-1][res2-1] && aux > 0.5 )	true_positives++;
		if(!Contact[res1-1][res2-1] && aux > 0.5 )	false_positives++;
        if(Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) < 23) true_positives_short++;
        if(!Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) < 23) false_positives_short++;
        if(Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) >= 23) true_positives_long++;
        if(!Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) >= 23) false_positives_long++;
        if(res1-1 < n/2 && res2-1 < n/2)
        {
                    if(Contact[res1-1][res2-1] && aux > 0.5 )   true_positives_h++;
                    if(!Contact[res1-1][res2-1] && aux > 0.5 )  false_positives_h++;
                    if(Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) < 23) true_positives_short_h++;
                    if(!Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) < 23) false_positives_short_h++;
                    if(Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) >= 23) true_positives_long_h++;
                    if(!Contact[res1-1][res2-1] && aux > 0.5 && abs(res1-res2) >= 23) false_positives_long_h++;

        }
	}

	/* Print total accuracy */
	if (VERBOSE)
	{
		printf("Contacts PPV: %.2lf %d %d\n",100.0*true_positives/(double)(true_positives+false_positives),true_positives,true_positives+false_positives);
        printf("Short-range Contacts PPV: %.2lf %d %d\n",100.0*true_positives_short/(double)(true_positives_short+false_positives_short),true_positives_short,true_positives_short+false_positives_short);
        printf("Long-range Contacts PPV: %.2lf %d %d\n",100.0*true_positives_long/(double)(true_positives_long+false_positives_long),true_positives_long,true_positives_long+false_positives_long);



        printf("Contacts at 50% PPV: %.2lf %d %d\n",100.0*true_positives_h/(double)(true_positives_h+false_positives_h),true_positives_h,true_positives_h+false_positives_h);
        printf("Short-range Contacts at 50% PPV: %.2lf %d %d\n",100.0*true_positives_short_h/(double)(true_positives_short_h+false_positives_short_h),true_positives_short_h,true_positives_short_h+false_positives_short_h);
        printf("Long-range Contacts at 50% PPV: %.2lf %d %d\n",100.0*true_positives_long_h/(double)(true_positives_long_h+false_positives_long_h),true_positives_long_h,true_positives_long_h+false_positives_long_h);



		fprintf(stderr,"Total protein contacts: %d\n",total_contacts);
	}

	fclose(input_map);
	fclose(input_fasta);
	fclose(input_aln);
	fclose(input_contacts);
	return 0;
}
