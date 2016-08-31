#ifndef SAM_UTIL_H
#define SAM_UTIL_H

#include <sam.h>
#include <map>
#include <deque>
#include <string>

// convert 4 bit encoding to IUPAC encoding (via 8 bit code.. well,)
const char* const iupac = "XACMGRSVTWYHKDBN";


// hold a sam/bam header as a pair of its tid and its size
// mapped by it's name to enable quick lookup of a tid for
// a given chromosome or contig name.
struct sam_hdr_map {
  std::map<std::string, std::pair<size_t, uint32_t> > hdr;
  std::string text;

  sam_hdr_map(){};
  sam_hdr_map(bam_hdr_t *bh){
    text = std::string( bh->text, bh->l_text );
    for(int i=0; i < bh->n_targets; ++i){
      hdr.insert(std::make_pair( std::string(bh->target_name[i]), std::make_pair( i, bh->target_len[i] )));
    }
  }
  int name2id( const char* name ){
    if( !hdr.count( std::string(name) ) )
      return(-1);
    return( hdr[ std::string(name) ].first );
  }
};

// holds simple variants.
// Warning; this struct does not attempt to manage its memory resources.
// In particular the user must avoid to copy objects around as
// the destructor will destroy the v_seq and v_qual objects.
// This is most eaasily handled by only using pointers to it
// If we implement reference counting then we can safely use it;
struct seq_variant {  
  enum varType {
    mismatch,
    insertion,
    deletion,
    null
  };
  varType v_type;
  unsigned int length;
  // space for v_seq and v_qual must be allocated during
  // construction. It is not null-terminated
  // and if length is 0, or we have a deletion
  // type it is simply set to 0. (note that these actually can have a non-zero length)
  char* v_seq;
  char* v_qual;
  // Then we have 4 + 4 + 8 = 16 bytes per
  // entry, plus whatever the length of the
  // vseq is.
  // Note that v_seq will have to be taken care of when appropriate...
  seq_variant(){
    v_type=null;
    length = 0;
    v_seq = 0;
    v_qual = 0;
  }
  seq_variant( varType v_type, unsigned int length, char* v_seq, char* v_qual ) :
    v_type(v_type), length(length), v_seq(v_seq), v_qual(v_qual) {
  }
  ~seq_variant(){
    delete []v_seq;
    delete []v_qual;
    // do not set these to 0; we want the code to segment fault if
    // the user performs a double delete. That will make it easier to
    // debug code.
  }
};


// a typedefs to hold variants by chromosome
// we hold these as pointers as the seq_variant object will 
// delete the sequence and quality data fields when destroyed.
typedef std::map<unsigned int, std::deque<seq_variant*> > chromVariants;

// This is an unsafe function that is required to
// be passed to hts_itr_querys.
// It's probably better to not use hts_itr_querys, but instead to make
// a wrapper function that takes a pointer to an sam_hdr_map, a string and 
// calls hts_itr_query instead.
// However, hts_itr_querys makes use of a nice region identifier so that
// one can either pass "chr" or "chr:beg-end" using all sorts of text styles.
// IMPORTANT:
// hdr must be a pointer to a sam_hdr_map as given above
int sam_hdr_name2id( void *hdr, const char* name );
/* { */
/*   sam_hdr_map *hm = (sam_hdr_map*)hdr; */
/*   return( hm->name2id( name ) ); */
/* } */

// this function should be passed to hts_itr_query and hts_itr_querys
// it will be called when hts_itr_next is called on the the iterator.
int read_record( BGZF *fp, void *data, void *r, int *tid,
		 int *beg, int *end);

void insertSeqVariant( chromVariants& c_var, uint32_t pos, seq_variant *var);

inline void scanCigarQualities(std::string& ref_seq, uint32_t r_pos, uint32_t q_pos, uint8_t *q_seq, uint8_t *q_qual,
			       uint32_t op_l, unsigned long *match_dist, unsigned long *mismatch_dist);

inline void scanCigarMatches(std::string& ref_seq, uint32_t& r_pos, uint32_t& q_pos, uint8_t *q_seq, uint8_t *q_qual,
			     uint32_t op_l, chromVariants& c_var, unsigned char minQ, uint32_t& var_count,
			     unsigned long *match_dist, unsigned long *mismatch_dist);

void checkMutations( bam1_t *bt, std::string& ref_seq, chromVariants& c_var, unsigned char minQ=0, 
		     unsigned long *match_dist=0, unsigned long *mismatch_dist=0);

void print_bam1(bam1_t *bt);

// some classes for convenience

// is it really useful to have a bamFile class  or should it be a hts_file?
class bamFile {
  const char *file_name;
  htsFile *bam_file;
  hts_idx_t *bf_index;
  bam_hdr_t *header;
  sam_hdr_map header_map;
  // we can also keep an iterator around; might be useful in order to maintain
  // state of the object.
  hts_itr_t *itr;

 public:
  bamFile();
  bamFile(const char *file_name);
  ~bamFile();

  // some getter functions
  htsFile* file();
  const hts_idx_t* bamIndex();
  const bam_hdr_t* bamHeader();
  sam_hdr_map headerMap();
  
  // and function to look for variants...
  // this uses hts_itr_querys to find the region. Does not guarantee that the mutations fall within
  // the specified region.
  void scanForVariants( std::map<std::string, std::string>& genome, std::string chromosomeRange,
			chromVariants& chrom_variants, unsigned char minQ,
			unsigned long *match_dist=0, unsigned long *mismatch_dist=0);
};

#endif
