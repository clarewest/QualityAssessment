#include<stdio.h>
#include<stdlib.h>

int main(int argc, char* argv[])
{
    int i,j,aux1,aux2;
    float score,minscore;
    FILE *input;

    input = fopen(argv[1],"r");
    minscore = atof(argv[3]);
    printf("LEN\t%d\n",atoi(argv[2]));
    while(fscanf(input,"%d %d %d %d %f",&i,&j,&aux1,&aux2,&score)!=EOF)
    {
        if(score>=minscore)
        {
        printf("CON\t%d\t%d\t%.3f\n",i-1,j-1,score);
        }
    }
    return 0;

}
