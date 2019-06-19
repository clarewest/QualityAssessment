#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

int main (int argc, char* argv[] ){
  if (argc<=3){
    std::cout << "Usage: "  << std::endl;
    exit(1); 
  }



  std::string filename = argv[1];
  std::string target;
  int start = atoi(argv[2]);
  int end = atoi(argv[3]);
  
  
  std::string decoy;
  float currentscore;
  float cumulativescore = 0;
  float meanscore; 
  int i = 0;

  std::ifstream is(filename);
  if (is.fail()){
    std::cout << "File does not exist!!!!!!!!!!!!" << std::endl;
    exit(1); 
  }
  std::string line; 

  while (getline(is,line)){
  
    std::string read;
    if(line[0] == 'F'){
      int index = 0;
      while(index < line.size() && line[index] != ' '){
        read += line[index];
        index++;
      }
      if(read == "File:"){
        index++;
        decoy.clear();
        while(index < line.size()){
          decoy += line[index];
          index++;
        }
      }
    }

    if(line[0] == 'A'){
      cumulativescore = 0.0;
      i = 0;
      char chain;
      int resnum; 
      char resname[64], assess[64], dummy[64];
      sscanf(line.c_str(),"%c %s %d %s %f %s\n", &chain, &resname, &resnum, &assess, &currentscore, &dummy);
      if(resnum>=start && resnum<=end){
        cumulativescore += currentscore;
        i++;
      }
      if(resnum==end){
        meanscore = cumulativescore / ((float) i);
        std::cout << decoy << " " << meanscore << std::endl; 
      }

    }

  }

  is.close();

}



/*
char decoy
float currentscore
float cumulativescore = 0
float meanscore 
int i = 0

open filename
while $1 != "File: " 
continue

decoy = line[1]

while $1 != "A "
continue

while $3 < start
continue

while $3 <= end
currentscore = $5
cumulativescore=+currentscore
i=+1

meanscore = cumulativescore/i

print target, decoy, meanscore
*/
