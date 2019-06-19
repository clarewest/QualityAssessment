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
	int i,j;
    int phi,phi_l,psi,psi_l;
    char c;
	char Fasta_Seq[2000],AUX[181],Proxy_Seq[2000],Proxy_SS[2000], SS[2000];
	double TA[2][2000],Pred_TA[2][2000];
    double diff_phi,diff_phi_l,diff_psi,diff_psi_l;
	FILE *input_fasta,*input_proxy,*input_dssp,*input_psipred,*input_angles,*input_pred_angles;
    
    diff_phi = 0.0;
    diff_phi_l = 0.0;
    diff_psi = 0.0;
    diff_psi_l = 0.0;
    phi=0;
    phi_l=0;
    psi=0;
    psi_l=0;

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


    for(i=0;i<strlen(Fasta_Seq);i++)
    {
        if(TA[0][i-1]>-370 && Pred_TA[0][i-1]>-370)
        {
                diff_phi += min(fabs(TA[0][i-1]-Pred_TA[0][i-1]), fabs(fabs(TA[0][i-1])+fabs(Pred_TA[0][i-1]) - 360));
                phi++;
                if(Proxy_SS[i-1]!='H'&&Proxy_SS[i-1]!='E')
                {   
                    diff_phi_l += min(fabs(TA[0][i-1]-Pred_TA[0][i-1]), fabs(fabs(TA[0][i-1])+fabs(Pred_TA[0][i-1]) - 360));
                    phi_l++;
                }
        }
        if(TA[1][i-1]>-370 && Pred_TA[1][i-1]>-370)
        {   
                diff_psi += min(fabs(TA[1][i-1]-Pred_TA[1][i-1]), fabs(fabs(TA[1][i-1])+fabs(Pred_TA[1][i-1]) - 360));
                psi++;
                if(Proxy_SS[i-1]!='H'&&Proxy_SS[i-1]!='E')
                {   
                    diff_psi_l += min(fabs(TA[1][i-1]-Pred_TA[1][i-1]), fabs(fabs(TA[1][j-1])+fabs(Pred_TA[1][i-1]) - 360));
                    psi_l++;
                }   
        }
    }

    printf("%.2f\t%.2f\t",diff_phi/phi,diff_psi/psi);
    printf("%.2f\t%.2f\n",diff_phi_l/phi_l,diff_psi_l/psi_l);

	fclose(input_fasta);
    fclose(input_proxy);
    fclose(input_dssp);
    fclose(input_psipred);
    fclose(input_angles);
    fclose(input_pred_angles);
	return 0;
}
