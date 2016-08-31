#include "fasta.h"
#include <fstream>
#include <iostream>

// this is primarily used to load full genomes into memory; hence the
// use of the genome identifier 
std::map<std::string, std::string> loadFasta(const char* fasta_file){
  const char *null_id = "null_sequence"; // this is a dirty hack.
  std::map<std::string, std::string> genome;
  std::ifstream is(fasta_file);
  if(!is){
    std::cerr << "loadFasta: Unable to open " << fasta_file << std::endl;
    genome.erase(null_id);
    return(genome);
  }
  std::string line;
  // this will allow me to avoid an if in the read part
  genome.insert(std::make_pair(null_id, ""));
  std::map<std::string, std::string>::iterator it = genome.find("null_sequence");
  while(getline( is, line)){
    if(line[0] == '>'){
      size_t p = line.find_first_of(" \t");
      std::string id;
      // I shouldn't need this check, as putting npos into 
      // substr should work, but I have a feeling that it screwed up before sometime;
      if(p != std::string::npos){
	id = line.substr(1, p-1);
      }else{
	id = line.substr(1);
      }
      genome.insert( std::make_pair(id, "") );
      it = genome.find( id );
      continue;
    }
    it->second += line;
  }
  genome.erase(null_id);
  return(genome);
}
