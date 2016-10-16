#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

#include "protein_pep_hash.pb.h"

using namespace std;

struct enzyme_cut_params {
   int         missed_cleavage;
   string      precut_amino;
   string      postcut_amino;
   string      prenocut_amino;
   string      postnocut_amino;
   int         semi_tryptic;
};

struct protein_pep_hash_db_hdr {
   int   version;
};

struct protein_data {
   string prd_id;
   string prd_pep_sequence;
};

typedef vector<protein_data *> data_t;

#define MAX_EDGES 30
#define MIN_FRAGMENT_MASS 600
#define MAX_FRAGMENT_MASS 600000

int fragment_mass[MAX_FRAGMENT_MASS];

static int ppsg_num_peptides = 0;

int pp_amino_acid_mass[MAX_EDGES] = {  71, //A 
               99, //B not there
               103, //C
               115, //D
               129, //E
               147, //F
               57, //G
               137, //H
               113, //I
               99, //J
               128, //K
               113, //L
               131, //M
               114, //N
               132, //O
               97, //P
               128, //Q
               156, //R
               87, //S
               101, //T
               99, //U
               99, //V
               186, //W
               99, //X
               163, //Y
               99 //Z 
            };

void phd_split_string(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
   /* Store the original string in the array, so we can loop the rest
    * of the algorithm. */
   tokens.push_back(str);

   // Store the split index in a 'size_t' (unsigned integer) type.
   size_t splitAt;
   // Store the size of what we're splicing out. One character is good enough.
   size_t splitLen = 1;
   // Create a string for temporarily storing the fragment we're processing.
   std::string frag;
   // Loop infinitely - break is internal.
   while(true)
   {
      /* Store the last string in the vector, which is the only logical
       * candidate for processing. */
      frag = tokens.back();
      /* The index where the split is. */
      splitAt = frag.find_first_of(splitBy);
      // If we didn't find a new split point...
      if(splitAt == string::npos)
      {
         // Break the loop and (implicitly) return.
         break;
      }
      /* Put everything from the left side of the split where the string
       * being processed used to be. */
      tokens.back() = frag.substr(0, splitAt + 1);
      /* Push everything from the right side of the split to the next empty
       * index in the vector. */
      tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
   }
}

void phd_read_protein_database(string file, peptide_hash_database::phd_file& data)
{
   // For every record we can read from the file, append it to our resulting data
   peptide_hash_database::phd_protein *record, *old_record;
   int pcount = 0;

   std::ifstream fstream(file.c_str());
   std::string line, prd_id, pr_seq;   

   string split_char = " ";
   vector<string> split_tokens;

   while (std::getline(fstream, line)) {
      if (line.at(0) == '>') {
         // terminate old sequence if it was started 
         //record = data.add_phdpro();
         pcount++;
         phd_split_string(line, split_char, split_tokens);
         record->set_phdpro_name(split_tokens[0]);
         record->set_phdpro_id(pcount);
      } else {
         // append or start a new one
         (record->mutable_phdpro_pepseq())->append(line);
      }
      split_tokens.clear();
   }

   return;
}

struct range {
   int start, length, missed, left, right;
};

void phd_basic_cut_pre_post(enzyme_cut_params params, const string protein_seq,
                              vector<range *> &splits)
{
   const string spres = params.precut_amino;
   const string sposts = params.postcut_amino;
   size_t pos = 0, length = protein_seq.length(), npres, nposts, cpos, start = pos;

   cout << "Protein length: " << length << " and pre split string " << spres << " and post split string is " << sposts << " and position is " << pos << endl;
   cout << "Protein sequence is: " << protein_seq << endl;

   while ((pos != std::string::npos) && (start < length)) {
      range *r = new range; r->start = pos;

      cpos = protein_seq.find_first_of(spres + sposts, start);

      npres = protein_seq.find_first_of(spres, start);
      nposts = protein_seq.find_first_of(sposts, start);

      if (cpos == std::string::npos) {
         r->length = length - r->start;
         splits.push_back(r);
         break;
      } else if (cpos == npres) {
         pos = cpos;
         r->length = pos - r->start;
         start = pos + 1;
      } else if (cpos == nposts) {
         // Check for K followed by Proline
         pos = cpos + 1;
         r->length = pos - r->start;
         start = pos;
      }
      cout << "Splits are " << npres << " " << nposts << " " << " final position is " << pos << endl;

      splits.push_back(r);
   }

   // Now we get the basic peptides
   cout << "We got the basic peptides based on the pre and post conditions" << endl;

   for (range *r: splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << endl;
   }

}

void phd_handle_post_merge(enzyme_cut_params params, const string protein_seq,
                           vector<range *> &splits, vector<range *> &nocut_splits)
{
   // We have to take care of pre nocut amino acids and postcut amino acids and merge
   const string nocut_spres = params.prenocut_amino;
   const string nocut_spost = params.postnocut_amino;

   vector<range *>::iterator split_iterator;

   cout << "Post no cut " << nocut_spost << endl;

   for (split_iterator = splits.begin(); split_iterator != splits.end(); split_iterator++) {
      range *r = *split_iterator;

      range *r_new = new range;
      r_new->start = r->start;
      r_new->length = r->length;
      r_new->missed = r_new->left = r_new->right = 0;

      while (split_iterator + 1 != splits.end()) {
         range *r_next = *(split_iterator + 1);
         //cout << "Split is " << r_next->start << " length is " << r_next->length << " the character is " << (protein_seq.c_str())[r_next->start] << endl;
         if (nocut_spost.find_first_of((protein_seq.c_str())[r_next->start]) != std::string::npos) {
            //cout << "Found a split" << endl;
            r_new->length += r_next->length;
            split_iterator++;
         } else break;
      }
      
      nocut_splits.push_back(r_new);
   }
   
   cout << "No cut splits" << endl;

   for (range *r: nocut_splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << endl;
   }

}

void phd_handle_missed_cleavage(enzyme_cut_params params, const string protein_seq,
                                 vector<range *> &nocut_splits)
{
   // We have to handle missing cleavage
   int missed_cleavage = params.missed_cleavage;
   cout << "Missed cleavage " << missed_cleavage << endl;

   if (missed_cleavage) {
      int num_cuts = nocut_splits.size();
      for (int i = 0; i < num_cuts - 1; i++) {
         range *r = nocut_splits.at(i), *r_next = nocut_splits.at(i+1);
         cout << "Peeking at the end " << (protein_seq.c_str())[r->start + r->length -1] << endl;
         if ((protein_seq.c_str())[r->start + r->length -1] == 'K') {
            range *r_new = new range;
            r_new->start = r->start;
            r_new->length = r->length + r_next->length;
            r_new->missed = r->length;
            r_new->left = r_new->right = 0;
            nocut_splits.push_back(r_new);
         }
      }
   }

   for (range *r: nocut_splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << 
               " and missed cleavage " << r->missed << " Left is " << r->left << 
               " Right is " << r->right << endl;
   }

}

void phd_add_missed_semi_tryptic(const string protein_seq, range* r,
                                 vector<range *> &nocut_splits)
{
   // Left tryptic
   cout << "Left trypic peptides for missed cleavages " << r->start << " " << r->length << " " << r->missed << endl;
   for (int pos = r->start + r->missed + 1; pos < (r->start + r->length); pos++) {
      range *r_new = new range;
      r_new->start = r->start;
      r_new->length= pos - r->start;
      r_new->missed = 1;
      r_new->right = 0;
      r_new->left = 1;
      nocut_splits.push_back(r_new);
   }

   // Right tryptic
   for (int start = r->start + 1; start < (r->start + r->missed); start++) {
      range *r_new = new range;
      r_new->start = start;
      r_new->length= r->start + r->length - start;
      r_new->missed = 1;
      r_new->right = 1;
      r_new->left = 0;
      nocut_splits.push_back(r_new);
   } 
}

void phd_add_semi_tryptic(const string protein_seq, range* r,
                                 vector<range *> &nocut_splits)
{
   // Left tryptic
   int start = r->start;
   for (int len = 1; len < r->length; len++) {
      range *r_new = new range;
      r_new->start = r->start;
      r_new->length= len;
      r_new->missed = 0;
      r_new->right = 0;
      r_new->left = 1;
      nocut_splits.push_back(r_new);
   }

   // right tryptic
   for (start = r->start + 1; start < (r->start + r->length -1); start++) {
      range *r_new = new range;
      r_new->start = start;
      r_new->length= r->start + r->length - start;
      r_new->missed = 0;
      r_new->right = 1;
      r_new->left = 0;
      nocut_splits.push_back(r_new);
   } 
}

void phd_handle_semi_tryptic(enzyme_cut_params params, const string protein_seq,
                              vector<range *> &nocut_splits)
{
   // We have to add sem-tryptic peptides
   int semi_tryptic = params.semi_tryptic;
   cout << "Semi-tryptic " << semi_tryptic << endl;

   if (semi_tryptic) {
      int num_count = nocut_splits.size();
      for (int i = 0; i < num_count; i++) {
         if ((nocut_splits.at(i))->missed) {
            phd_add_missed_semi_tryptic(protein_seq, nocut_splits.at(i), nocut_splits);
         } else {
            phd_add_semi_tryptic(protein_seq, nocut_splits.at(i), nocut_splits);
         }
      }
   }

   for (range *r: nocut_splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << 
               " and missed cleavage " << r->missed << " Left is " << r->left << 
               " Right is " << r->right << " and the peptide is " << protein_seq.substr(r->start, r->length) << endl;
   }
}

void phd_split_protein_sequence_peptides(enzyme_cut_params params, 
                                          peptide_hash_database::phd_protein pro_seq)
{
   const string protein_seq = pro_seq.phdpro_pepseq();
   vector<range *> splits, nocut_splits;

   phd_basic_cut_pre_post(params, protein_seq, splits);

   phd_handle_post_merge(params, protein_seq, splits, nocut_splits);

   phd_handle_missed_cleavage(params, protein_seq, nocut_splits);

   phd_handle_semi_tryptic(params, protein_seq, nocut_splits);

   // clean up the memory for splits, nocut_splits
}

void phd_params_copy_from_header(peptide_hash_database::phd_header hdr, enzyme_cut_params &cut_params)
{
   cut_params.precut_amino = hdr.phdhdr_precut_amino();
   cut_params.postcut_amino = hdr.phdhdr_postcut_amino();
   cut_params.prenocut_amino = hdr.phdhdr_prenocut_amino();
   cut_params.postnocut_amino = hdr.phdhdr_postnocut_amino();
   cut_params.missed_cleavage = hdr.phdhdr_missed_cleavage();
   cut_params.semi_tryptic = hdr.phdhdr_semi_tryptic();
}

void phd_add_peptide_hash_database (peptide_hash_database::phd_file &pfile)
{
   const peptide_hash_database::phd_header hdr = pfile.phdhdr();

   enzyme_cut_params cut_params;
   phd_params_copy_from_header(hdr, cut_params);

   //for (int i = 0; i < pfile.phdpro_size(); i++) {
     // cout << "Splitting peptides for protein id : " << i << endl;
      //phd_split_protein_sequence_peptides(cut_params, pfile.phdpro(i));
   //}

   // If the parameter is semi-tryptic, add all left and right semi-tryptic peptides
}

int phd_save_hash(peptide_hash_database::phd_file &pfile)
{
   std::string filename = "phd.txt";

   /*
   peptide_hash_database::phd_header hdr;

   // Create header which is mandatory
   hdr = pfile.phdhdr();
   */
   cout << "Creating file" << endl;

   std::ofstream ofs;
   ofs.open (filename, ios::out | ios::trunc | ios::binary);
   if (!pfile.SerializeToOstream(&ofs)) {
      cout << "Cannot write to hash file" << endl;
      exit(1);
   }

   cout << "Wrote file" << endl;
   ofs.close();
}

int phd_load_hash()
{
}

void phd_params_copy_to_header(peptide_hash_database::phd_header *phdr, 
                                 string prot_file, int semi_tryptic,
                                 string precut_amino, string postcut_amino,
                                 string prenocut_amino, string postnocut_amino,
                                 int internal_lysine)
{
   phdr->set_phdhdr_version(1);
   phdr->set_phdhdr_protein_source_filename(prot_file);
   phdr->set_phdhdr_protein_source_file_digest("DEADBEEF");
   phdr->set_phdhdr_num_proteins(10);
   phdr->set_phdhdr_hash_file_name("DEADBEEF");
   phdr->set_phdhdr_hash_file_digest("DEADBEEF");

   phdr->set_phdhdr_precut_amino(precut_amino);
   phdr->set_phdhdr_postcut_amino(postcut_amino);
   phdr->set_phdhdr_prenocut_amino(prenocut_amino);
   phdr->set_phdhdr_postnocut_amino(postnocut_amino);
   phdr->set_phdhdr_missed_cleavage(internal_lysine);
   phdr->set_phdhdr_semi_tryptic(semi_tryptic);
}

void phd_read_cmdline(char *argv[], 
                        peptide_hash_database::phd_header *phdr)
{
   string prot_file(argv[1]);
   int semi_tryptic = atoi(argv[2]);
   string precut_amino(argv[3]);
   string postcut_amino(argv[4]);
   string prenocut_amino(argv[5]);
   string postnocut_amino(argv[6]);
   int internal_lysine = atoi(argv[7]);
   string hash_file(argv[8]);

   // Validate the parameters: no intersection of pre and post and break and no-break

   cout << endl 
         << "Parameters of the run: " << endl
         << "Protein database file name is: " << prot_file << endl
         << "Tryptic/Semi-tryptic run: " << semi_tryptic << endl
         << "Pre break amino acids: " << precut_amino << endl
         << "Post break amino acids: " << postcut_amino << endl
         << "Pre no break amino acids: " << prenocut_amino << endl
         << "Post no break amino acids: " << postnocut_amino << endl
         << "Internal lysine (missed cleavage): " << internal_lysine << endl 
         << "Hash file given is: " << hash_file << endl 
         <<endl;

   phd_params_copy_to_header(phdr, prot_file, semi_tryptic,
                              precut_amino, postcut_amino,
                              prenocut_amino, postnocut_amino,
                              internal_lysine);
}

void phd_print_params(peptide_hash_database::phd_header hdr)
{
   cout << endl 
         << "Parameters of the file: " << endl
         << "Protein database file name is: " << hdr.phdhdr_precut_amino()<< endl
         << "Tryptic/Semi-tryptic run: " << hdr.phdhdr_semi_tryptic() << endl
         << "Pre break amino acids: " << hdr.phdhdr_precut_amino() << endl
         << "Post break amino acids: " << hdr.phdhdr_postcut_amino() << endl
         << "Pre no break amino acids: " << hdr.phdhdr_prenocut_amino() << endl
         << "Post no break amino acids: " << hdr.phdhdr_postnocut_amino() << endl
         << "Internal lysine (missed cleavage): " << hdr.phdhdr_missed_cleavage() << endl 
         <<endl;

}

int phd_read_hash_file_and_compare(char *hash_file, peptide_hash_database::phd_file pfile)
{
   cout << "File name is " << hash_file << endl;
   fstream input(hash_file, ios::in | ios::binary);
   if (!input) {
      cout << hash_file << ": File not found.  Creating a new file." << endl;
      return 1;
   }

   if (pfile.ParseFromIstream(&input)) {
      cout << "Read the hash file: " << hash_file << endl;
      peptide_hash_database::phd_header phdr = pfile.phdhdr();
      phd_print_params(phdr);
      return 0;
   } else {
      cout << "File read has problems" << endl;
   }

   return 1;
}

int phd_save_load_hash(char* argv[], peptide_hash_database::phd_file &pfile)
{
   // Check for the existance of the file pdb_cut_ignore
   // If the file is present
   //    read the file and create in memory structures
   // else
   //    create the data structures and save

   phd_read_cmdline(argv, pfile.mutable_phdhdr());

   int not_matched = phd_read_hash_file_and_compare(argv[8], pfile);
   if (not_matched) {
      phd_save_hash(pfile);
   }

   /*
   if (phd_read_hash_file_and_compare(argv[8], pfile)) {
      // Here is the file containing the data. Read it into data.
      phd_read_protein_database(argv[1], pfile);

      for (int i = 0; i < MAX_FRAGMENT_MASS; i++) {
      }

      phd_add_peptide_hash_database(pfile);

      cout << "Read this many records: " << pfile.phdpro_size() << endl;

      phd_save_hash(pfile);
   } 
   */

}

int main(int argc, char *argv[]) 
{
   peptide_hash_database::phd_file pfile;

   /*
    Usage: <protein_database> <missed_cleavage> <cut_aminoacid> <ignore_prolin> <peptide_hash_file>
   */

   GOOGLE_PROTOBUF_VERIFY_VERSION;

   if (argc != 9) {
      cout << "./a.out <protein_database> <semi/tryptic> <precut_aa> " 
            << "<postcut_aa> <prenocut_aa> <postnocut_aa> <internal_lysine> <saved_hash_file>" << endl;
      exit(1);
   }
   
   // Here is the data we want.
   phd_save_load_hash(argv, pfile);

/*
   // For each protein, split the peptide sequence based on the cut amino acid
   // For each cut add to the graph 
   ppsg_add_protein_database_hash(data);

   for (int i = MIN_FRAGMENT_MASS; i < MAX_FRAGMENT_MASS; i++) {
      cout << "Number of peptide sequences at this mass are " << i << " " << fragment_mass[i] << endl;
   }
*/
   google::protobuf::ShutdownProtobufLibrary();
}
