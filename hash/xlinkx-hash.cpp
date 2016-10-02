#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

#include "protein_pep_hash.pb.h"

using namespace std;

struct protein_pep_hash_db_hdr {
   int   version;
};

struct protein_data {
   string prd_id;
   string prd_pep_sequence;
};

typedef vector<protein_data *> data_t;

struct protein_pep_sequence {
   protein_data*  pps_protein;
   int      pps_start, pps_length;
};

#define MAX_EDGES 30
#define MIN_FRAGMENT_MASS 600
#define MAX_FRAGMENT_MASS 600000

//vector<protein_pep_sequence*> fragment_mass[MAX_FRAGMENT_MASS];
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

void ppsg_split_string(std::string str, std::string splitBy, std::vector<std::string>& tokens)
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

void ppsg_split_peptides(protein_data* protein, std::string splitBy, std::vector<protein_pep_sequence*> &tokens)
{
   // Store the split index in a 'size_t' (unsigned integer) type.
   size_t splitAt;
   // Store the size of what we're splicing out. One character is good enough.
   size_t splitLen = 1;
   // Create a string for temporarily storing the fragment we're processing.
   std::string frag = protein->prd_pep_sequence;

   int cur_start = 0;

   // Loop infinitely - break is internal.
   while(true)
   {
      protein_pep_sequence* pps = new protein_pep_sequence;
      ppsg_num_peptides++;

      pps->pps_protein = protein;
      pps->pps_start = cur_start;
      tokens.push_back(pps);

      splitAt = frag.find_first_of(splitBy);
      // If we didn't find a new split point...
      if((splitAt == string::npos) || (splitAt == frag.size() - 1))
      {
         // Break the loop and (implicitly) return.
         pps->pps_length = frag.size();
         break;
      }
      pps->pps_length = splitAt+splitLen;
      /* Push everything from the right side of the split to the next empty
       * index in the vector. */
      int length = frag.length()-(splitAt+splitLen);
      frag = frag.substr(splitAt+splitLen, length);
      cur_start += pps->pps_length;
   }
}

void ppsg_read_protein_database(string file, data_t& data)
{
   // For every record we can read from the file, append it to our resulting data
   protein_data *record, *old_record;

   std::ifstream fstream(file.c_str());
   std::string line, prd_id, pr_seq;   

   string split_char = " ";
   vector<string> split_tokens;

   while (std::getline(fstream, line)) {
      if (line.at(0) == '>') {
         // terminate old sequence if it was started 
         record = new protein_data;
         data.push_back(record);
         ppsg_split_string(line, split_char, split_tokens);
         record->prd_id = split_tokens[0];
      } else {
         // append or start a new one
         (record->prd_pep_sequence).append(line);
      }
      split_tokens.clear();
   }

   return;
}

#define PRINT_PROTEIN_STATUS 10000


void ppsg_add_peptide_sequence_hash (protein_pep_sequence *peptide)
{
   int pos = peptide->pps_start, mass = 0;
   for (; pos < (peptide->pps_length + peptide->pps_start); pos++) {
      char amino_acid = ((peptide->pps_protein)->prd_pep_sequence).at(pos);
      mass += pp_amino_acid_mass[amino_acid - 'A'];
      fragment_mass[mass]++;
   }
}

void ppsg_add_protein_database_hash (data_t data)
{
   string split_char = "KR";
   vector<protein_pep_sequence*> split_peptides;

   std::ofstream ofs;
   ofs.open("peptide_sequence_dump.txt", std::ofstream::out | std::ofstream::app);

   int count = 0;
   for (protein_data *protein : data) {
      // cout << "Protein number is " << count << " and it is " << protein->prd_pep_sequence << endl;
      if (count % PRINT_PROTEIN_STATUS == 0)
         cout << "Reading protein " << count << endl;

      ppsg_split_peptides(protein, split_char, split_peptides);
      for (protein_pep_sequence* peptide : split_peptides) {
//       cout << "Peptide is " << peptide->pps_start << endl;
         ppsg_add_peptide_sequence_hash(peptide);
         ofs << "Protein ID: " << peptide->pps_protein->prd_id << " Start = " << 
            peptide->pps_start << " End = " << peptide->pps_start + peptide->pps_length - 1 << " " <<  
            (peptide->pps_protein->prd_pep_sequence).substr(peptide->pps_start, peptide->pps_length) << endl;
      }
      split_peptides.clear();
      count++;
   }
}

int main(int argc, char *argv[]) 
{
   /*
    Usage: <protein_database> <missed_cleavage> <cut_aminoacid> <ignore_prolin> <peptide_hash_file>
   */

   GOOGLE_PROTOBUF_VERIFY_VERSION;

   if (argc != 6) {
      cout << "./a.out <protein_database> <missed_cleavage> <cut_aminoacid> <ignore_prolin> <peptide_hash_file>" << endl;
      exit(1);
   }
   
   // Here is the data we want.
   data_t data;

   // Here is the file containing the data. Read it into data.
   ppsg_read_protein_database("protein_sequence_database.txt", data);

   // Otherwise, list some basic information about the file.
   cout << "Your input CSV file contains " << data.size() << " records.\n";

   for (protein_data *protein : data) {
      cout << "Protein is " << protein->prd_id << endl;
      cout << "Its amino acid sequence is " << protein->prd_pep_sequence << endl;
   }

/*
   // For each protein, split the peptide sequence based on the cut amino acid
   // For each cut add to the graph 
   ppsg_add_protein_database_hash(data);

   for (int i = MIN_FRAGMENT_MASS; i < MAX_FRAGMENT_MASS; i++) {
      cout << "Number of peptide sequences at this mass are " << i << " " << fragment_mass[i] << endl;
   }
*/
}
