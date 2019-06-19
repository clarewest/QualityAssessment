#include<stdio.h>
#include<stdlib.h>

int main(int argc, char* argv[])
{
    int i,j,aux;
    float score;
    FILE *input;

    input = fopen(argv[1],"r");
    printf("LEN\t%d\n",atoi(argv[2]));
    while(fscanf(input,"%d %d %f %d",&i,&j,&score,&aux)!=EOF)
    {
        printf("CON\t%d\t%d\t%.3f\n",i-1,j-1,score);
    }
    return 0;

}
