#include <string>

#include "protein_pep_hash.pb.h"

struct enzyme_cut_params {
   int         missed_cleavage;
   string      precut_amino;
   string      postcut_amino;
   string      prenocut_amino;
   string      postnocut_amino;
   int         semi_tryptic;
};

struct protein_hash_db_ {
   peptide_hash_database::phd_file phd_file_entry;
   vector<string*>* phd_get_peptides_ofmass(int mass);
};

typedef protein_hash_db_* protein_hash_db_t;

const protein_hash_db_t phd_retrieve_hash_db (const char *protein_file, 
                              enzyme_cut_params params, const char *phd_file);

