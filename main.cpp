#include "fasta.h"
#include "sam_util.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

// specify location of the bam file, genome file, chromosome, and minQ
// to identify variants.. 
int main(int argc, char **argv){
  if(argc < 4){
    std::cerr << "Usage: " << argv[0] << " bam_file genome_file minQ" << std::endl;
    exit(1);
  }
  
  bamFile bf( argv[1] );
  
  // htsFile *bf = hts_open( argv[1], "rb" );
  // hts_idx_t *bf_index = hts_idx_load( argv[1], HTS_FMT_BAI );
  
  // if( !hts_check_EOF( bf ) ){
  //   std::cerr << "Failed to find end of " << argv[1] << " corrupted file?" << std::endl;
  //   exit(2);
  // }
  
  // if( bf_index == NULL ){
  //   std::cerr << "Failed to load index for " << argv[1] << std::endl;
  //   exit(3);
  // }
  
  // then we wish to load the genome
  std::map<std::string, std::string> genome = loadFasta( argv[2] );
  if(!genome.size()){
    std::cerr << "Failed to load the genome sequence" << std::endl;
  }else{
    std::cout << "Loaded genome with " << genome.size() << " sequences" << std::endl;
  }
  
  unsigned int minQ = atoi( argv[3] );

  // bam_hdr_t *header = bam_hdr_read( bf->fp.bgzf );
  // sam_hdr_map bam_header( header );
  // bam1_t *bam_entry = bam_init1();
  
  // const char *chrom = argv[3];
  // if(!genome.count( chrom )){
  //   std::cerr << "No sequence data for chromosome : " << chrom << std::endl;
  //   exit(4);
  // }
  
  // if( bam_header.name2id( chrom ) < 0 ){
  //   std::cerr << "No bam data for chromosome : " << chrom << std::endl;
  //   exit(5);
  // }

  // sam_hdr_map bam_header = bf.headerMap();
  // hts_itr_t *itr = hts_itr_querys( bf.bamIndex(), "LG01", &sam_hdr_name2id, &bam_header, &hts_itr_query, &read_record );
  // bam1_t *bam_entry = bam_init1();
  // unsigned int count = 0;
  // htsFile *f = bf.file();
  // while( hts_itr_next( f->fp.bgzf, itr, bam_entry, f) > 0 && ++count < 10 )
  //   print_bam1(bam_entry);
  
  // exit(0);

  // hts_itr_t *itr = hts_itr_querys( bf_index, chrom, &sam_hdr_name2id, &bam_header, &hts_itr_query, &read_record );
  // chromVariants cVar;
  // unsigned int count = 0;
  // while( hts_itr_next( bf->fp.bgzf, itr, bam_entry, bf ) >= 0 ){
  //   //    print_bam1( bam_entry );
  //   checkMutations( bam_entry, genome[ chrom ], cVar, minQ );
  //   if(!( ++count % 1000000 )){
  //     std::cout << count << std::endl;
  //     print_bam1( bam_entry );
  //   }
  // }
  std::map<std::string, chromVariants> genomeVariants;
  unsigned long *match_dist = new unsigned long[256];
  unsigned long *mismatch_dist = new unsigned long[256];
  memset( match_dist, 0, sizeof(unsigned long) * 256 );
  memset( mismatch_dist, 0, sizeof(unsigned long) * 256 );
  //std::map<std::string, std::string>::iterator it=genome.find("LG01");
  for(std::map<std::string, std::string>::iterator it=genome.begin(); it != genome.end(); it++){
    chromVariants cVar;
    genomeVariants.insert(make_pair( it->first, cVar));
    std::map<std::string, chromVariants>::iterator cit = genomeVariants.find( it->first );
    bf.scanForVariants( genome, it->first, cit->second, minQ, match_dist, mismatch_dist );
    std::cout << "chromosome: " << it->first <<  "  " << cit->second.size() << " variants" << std::endl;
  }
  for(unsigned int i=0; i < 256; ++i)
    std::cout << i << " :\t" << match_dist[i] << "\t" << mismatch_dist[i] << std::endl;
}
