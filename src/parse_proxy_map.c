#include<stdio.h>
#include<stdlib.h>

int main(int argc, char* argv[])
{
    int length,i,j;
    FILE *input;

    input = fopen(argv[1],"r");
    fscanf(input,"%d",&length);
    printf("LEN\t%d\n",length);
    while(fscanf(input,"%d %d",&i,&j)!=EOF)
    {
        printf("CON\t%d\t%d\t1\n",i,j);
    }
    return 0;

}
