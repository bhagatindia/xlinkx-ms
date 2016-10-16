#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

#include "protein_pep_hash.pb.h"

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
   peptide_hash_database::phd_parameters *pparams = phdr->mutable_phdhdr_params();

   phdr->set_phdhdr_protein_source_filename(prot_file);

   pparams->set_phdparam_precut_amino(precut_amino);
   pparams->set_phdparam_postcut_amino(postcut_amino);
   pparams->set_phdparam_prenocut_amino(prenocut_amino);
   pparams->set_phdparam_postnocut_amino(postnocut_amino);
   pparams->set_phdparam_missed_cleavage(internal_lysine);
   pparams->set_phdparam_semi_tryptic(semi_tryptic);
}

void phd_read_cmdline(char *argv[], 
                        peptide_hash_database::phd_header *phdr)
{
/*
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
*/
}

void phd_print_params(peptide_hash_database::phd_parameters params)
{
   cout << endl 
         << "Parameters of the file: " << endl
         << "Protein database file name is: " << params.phdparam_precut_amino()<< endl
         << "Tryptic/Semi-tryptic run: " << params.phdparam_semi_tryptic() << endl
         << "Pre break amino acids: " << params.phdparam_postcut_amino() << endl
         << "Post break amino acids: " << params.phdparam_postcut_amino() << endl
         << "Pre no break amino acids: " << params.phdparam_prenocut_amino() << endl
         << "Post no break amino acids: " << params.phdparam_postnocut_amino() << endl
         << "Internal lysine (missed cleavage): " << params.phdparam_missed_cleavage() << endl 
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

   if(!pfile.ParseFromIstream(&input)) {
      cout << "Read the hash file: " << hash_file << endl;
      peptide_hash_database::phd_header phdr = pfile.phdhdr();
      peptide_hash_database::phd_parameters params = phdr.phdhdr_params();
      phd_print_params(params);
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

   phd_read_hash_file_and_compare(argv[1], pfile);
   phd_save_hash(pfile);

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

   if (argc != 2) {
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
