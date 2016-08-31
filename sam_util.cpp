#include "sam_util.h"
#include <iostream>

int sam_hdr_name2id( void *hdr, const char* name ){
  sam_hdr_map *hm = (sam_hdr_map*)hdr;
  return( hm->name2id( name ) );
}

int read_record( BGZF *fp, void *data, void *r, int *tid,
		 int *beg, int *end){
  bam1_t *bt = (bam1_t*)r;
  // I don't know what is the meaning of the br code, but should be ok.
  int br = bam_read1((BGZF*)fp, bt); 
  *tid = bt->core.tid;
  *beg = bt->core.pos;
  // this is strictly speaking incorrect. To give the correct answer we need to 
  // parse the cigar string. There may be a function / macro that does that.
  // However, it seems unreasonable to parse the SAM line here without using it.
  *end = bt->core.pos + bt->core.l_qseq;  
  return(br);
}

void insertSeqVariant( chromVariants& c_var, uint32_t pos, seq_variant *var){
  chromVariants::iterator it = c_var.find( pos );
  if(it == c_var.end()){
    std::deque<seq_variant*> vars;
    vars.push_back( var );
    c_var.insert(std::make_pair( pos, vars ));
  }else{
    it->second.push_back( var );
  }
}
// this function incremements counters for the distribution of qualities
// It may be possible to parallelise this function, but that depends on whether bam_seqi is re-entrant
// and whether the increment operator can be considered atomic (which I think it can't).
// safest way is probably to do individual chromosomes at the same time.
inline void scanCigarQualities(std::string& ref_seq, uint32_t r_pos, uint32_t q_pos, uint8_t *q_seq, uint8_t *q_qual,
			       uint32_t op_l, unsigned long *match_dist, unsigned long *mismatch_dist){
  for(uint32_t i=0; i < op_l; ++i){
    if( ref_seq[ r_pos + i ] == iupac[ bam_seqi( q_seq, q_pos + i ) ] ){
      match_dist[ q_qual[ q_pos + i ] ]++;
    }else{
      mismatch_dist[ q_qual[ q_pos + i] ]++;
    }
  }
}

inline void scanCigarMatches(std::string& ref_seq, uint32_t& r_pos, uint32_t& q_pos, uint8_t *q_seq, uint8_t *q_qual,
			     uint32_t op_l, chromVariants& c_var, unsigned char minQ, uint32_t& var_count,
			     unsigned long *match_dist, unsigned long *mismatch_dist){
  if(match_dist && mismatch_dist)
    scanCigarQualities( ref_seq, r_pos, q_pos, q_seq, q_qual, op_l, match_dist, mismatch_dist );

  for(uint32_t i = 0; i < op_l; ++i){
    uint32_t v_length = 1;
    if( ref_seq[ r_pos ] != iupac[ bam_seqi( q_seq, q_pos ) ] && q_qual[ q_pos ] >= minQ ){
      while( r_pos + v_length < ref_seq.size() && q_pos + v_length < op_l && ref_seq[ r_pos + v_length ] != iupac[ bam_seqi(q_seq, q_pos + v_length) ] 
	     && q_qual[ q_pos + v_length ] >= minQ){
	++v_length;
      }
      char *var_seq = new char[ v_length ];
      char *var_q = new char[ v_length ];
      for(uint32_t j=0; j < v_length; ++j){
	var_seq[ j ] = iupac[ bam_seqi( q_seq, q_pos + j ) ];
	var_q[ j ] = q_qual[ q_pos + j ];
      }
      seq_variant *var = new seq_variant( seq_variant::mismatch, v_length, var_seq, var_q );
      insertSeqVariant( c_var, r_pos + i, var );
      var_count++;
    }
    q_pos += v_length;
    r_pos += v_length;
    i += (v_length - 1);  // as it will anyway get incremented by the loop .. ??
  }
}


// this is ugly, but I'm not sure how to avoid it
// to be honest, branch prediction should work reasonably well here, but still 
void checkMutations( bam1_t *bt, std::string& ref_seq, chromVariants& c_var, unsigned char minQ,
		     unsigned long *match_dist, unsigned long *mismatch_dist){
  // bam_cigar_op gives a value of 0 to 9 with the values mapping to
  // operations : "MIDNSHP=XB"
  uint32_t *cigar = bam_get_cigar( bt );
  uint32_t r_pos = bt->core.pos;
  uint32_t q_pos = 0;
  uint8_t *q_seq = bam_get_seq( bt );
  uint8_t *q_qual = bam_get_qual( bt );
  uint32_t var_count = 0;
  for( unsigned int i=0; i < bt->core.n_cigar; ++i ){
    uint32_t op = bam_cigar_op( *(cigar + i) );
    uint32_t op_l = bam_cigar_oplen( *(cigar + i) );
    char *var_seq = 0;
    char *var_q = 0;
    seq_variant *var = 0;
    //    uint32_t v_length = 0;
    switch (op) {
    case 0: // M  match, can be mismatch go through all and identify matches
      scanCigarMatches( ref_seq, r_pos, q_pos, q_seq, q_qual, op_l,
			c_var, minQ, var_count, match_dist, mismatch_dist );
      break;
    case 1: // I Insertion
      // and insert an insertion event here.
	// these are more difficult to apply a quality call to. Because we may
	// end up splitting or shortening the insert into several places..
	var_seq = new char[ op_l ];
	var_q = new char[ op_l ];
	for(unsigned int i=0; i < op_l; ++i){
	  var_seq[i] = iupac[ bam_seqi(q_seq, q_pos + i) ];
	  var_q[i] = q_qual[ q_pos + i ];
	}
	var = new seq_variant( seq_variant::insertion, op_l, var_seq, var_q );
	insertSeqVariant(c_var, r_pos, var);
	var_count++;
	q_pos += op_l;
	break;
      case 2: // D Deletion
	var = new seq_variant( seq_variant::deletion, op_l, 0, 0 );
	insertSeqVariant(c_var, r_pos, var);
	var_count++;
	r_pos += op_l;
	// and insert deletion even here
	break;
      case 3: // N Region skipped from reference (similar to D, but long)
	var = new seq_variant( seq_variant::deletion, op_l, 0, 0 );
	insertSeqVariant(c_var, r_pos, var);
	var_count++;
	r_pos += op_l;
	// and insert a deletion event here
	break;
      case 4: // S soft clipped (sequence present in query sequence but not aligned)
	q_pos += op_l;
	r_pos += op_l;  // depends on what is defined as the pos.
	// there is no variant, so leave it here
	break;
      case 5: // H hard clipped (sequence not present in the query sequence)
	r_pos += op_l;
	// also no variant
	break;
      case 6: // P padding (silent deletion from padded reference ?)
	r_pos += op_l; // though I'm not really sure here
	break;
      case 7: // = a real match
	if( match_dist && mismatch_dist )
	  scanCigarQualities( ref_seq, r_pos, q_pos, q_seq, q_qual, op_l, match_dist, mismatch_dist );
	q_pos += op_l;
	r_pos += op_l;
	break;
      case 8: // X a real mismatch
	// insert a variant of length op_l, then increment counters.
	scanCigarMatches( ref_seq, r_pos, q_pos, q_seq, q_qual, op_l,
			  c_var, minQ, var_count, match_dist, mismatch_dist );
	break;
      case 9: // B backspace, like N, but negative.. ?? rarely used
	r_pos -= op_l;  // but what does this mean
	break;
      default:
	std::cerr << "Unknown cigar operation : " << op << std::endl;
      }
  }
  // std::cout << "\nInserted total of " << var_count << " variants" << std::endl;
}


void print_bam1(bam1_t *bt){
  std::cout << "************************************" << std::endl;;
  std::cout << "tid: " << bt->core.tid << "  pos: " << bt->core.pos << "  bin: " << bt->core.bin
	    << "  qual: " << bt->core.qual << "  cigar length: " << bt->core.n_cigar << std::endl;
  std::cout << "mate position: " << bt->core.mtid << ", " << bt->core.mpos << "  isize? " << bt->core.isize << std::endl;
  std::cout << "query name: " << bt->data << std::endl;
  uint32_t *cigar = bam_get_cigar( bt );
  for(unsigned int i=0; i < bt->core.n_cigar; ++i){
    std::cout << bam_cigar_opchr( *(cigar + i) ) << bam_cigar_oplen( *(cigar + i) ) << std::endl;
  }
  uint8_t *seq = bam_get_seq( bt );
  uint8_t *qual = bam_get_qual( bt );
  for(int32_t i=0; i < bt->core.l_qseq; ++i){
    std::cout << (char)(33 + qual[i]);
  }
  std::cout << std::endl;
  for(int32_t i=0; i < bt->core.l_qseq; ++i){
    std::cout << (int)qual[i] << ",";
  }
  std::cout << std::endl;
  for(int32_t i=0; i < bt->core.l_qseq; ++i){
    std::cout << iupac[ bam_seqi( seq, i ) ];
  }
  std::cout << std::endl;
  std::cout << "************************************\n" << std::endl;;
}

bamFile::bamFile(){
  bam_file = 0;
  bf_index = 0;
  header = 0;
  itr = 0;
}

bamFile::bamFile( const char *file_name ) : file_name(file_name)
{
  // better to have a small initialise to null function..
  bam_file = 0;
  bf_index = 0;
  header = 0;
  itr = 0;
  
  // open file and set up headers and indices
  bam_file = hts_open( file_name, "rb" );
  bf_index = hts_idx_load( file_name, HTS_FMT_BAI );
  
  if( !hts_check_EOF( bam_file ) ){
    std::cerr << "bamFile::bamFile Failed to find end of " << file_name << " corrupted file?" << std::endl;
    return;
  }
  
  if( bf_index == NULL ){
    std::cerr << "bamFile::bamFile Failed to load index for " << file_name << std::endl;
    return;
  }

  header = bam_hdr_read( bam_file->fp.bgzf );
  header_map = sam_hdr_map( header );
  // we have no need to set up the iterator here. We can leave it as 0.
}

bamFile::~bamFile()
{
  // note that calling hts_close on a null pointer is invalid. However safe for the
  // index and the iterator. To be safe we'll destroy in reverse order of creation:
  hts_itr_destroy( itr );
  hts_idx_destroy( bf_index );
  bam_hdr_destroy( header );
  if(bam_file)
    hts_close( bam_file );
}

htsFile* bamFile::file(){
  return( bam_file );
}

const hts_idx_t* bamFile::bamIndex(){
  return( bf_index );
}

const bam_hdr_t* bamFile::bamHeader(){
  return( header );
}

sam_hdr_map bamFile::headerMap(){
  return( header_map );
}

void bamFile::scanForVariants( std::map<std::string, std::string>& genome, std::string chromosomeRange,
			       chromVariants& chrom_variants, unsigned char minQ,
			       unsigned long *match_dist, unsigned long *mismatch_dist){
  // first determine the chromosome identity
  size_t col_pos = chromosomeRange.find_first_of(":");
  std::string chrom = chromosomeRange.substr(0, col_pos);
  
  // and check that we have actually got some sequence for it
  if(!genome.count(chrom)){
    std::cerr << "bamFile::scanForVariants : no sequence defined for " << chrom << std::endl;
    return;
  }
  // is this safe on a null index ?
  if(!bf_index){
    std::cerr << "bamFile::scanForVariants : no index defined for " << file_name << std::endl;
    return;
  }
  itr = hts_itr_querys( bf_index, chromosomeRange.c_str(), &sam_hdr_name2id, &header_map, &hts_itr_query, &read_record );
  bam1_t *bam_entry = bam_init1();
  unsigned int count = 0;
  std::string& seq = genome[chrom]; 
  while( hts_itr_next( bam_file->fp.bgzf, itr, bam_entry, bam_file ) >= 0 ){
    checkMutations( bam_entry, seq, chrom_variants, minQ, match_dist, mismatch_dist );
    if(!( ++count % 100000 ))
      std::cerr << "checked " << count << " alignments" << std::endl;
  }
  std::cerr << "checked a total of " << count << " alignments" << std::endl;
  bam_destroy1( bam_entry );
  hts_itr_destroy( itr );
  itr = 0;
}
