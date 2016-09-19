#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

#define MAX_EDGES 26

// B is not present as an amino acid. So, reusing it to store protein name
#define PEPG_NAME_NODE_INDEX 1
// Z is not present as an amino acid. So, reusing it to store flags
#define PEPG_FLAG_NODE_INDEX 25

#define PEPG_FLAG_TERMINAL 0x01

#define IS_PEPG_NODE_TERMINAL(x) ((((protein_pep_sequence_graph *)x->pg_edge_neighbor)[PEPG_FLAG_NODE_INDEX]) & (PEPG_FLAG_TERMINAL))
#define SET_PEPG_NODE_TERMINAL(x) ((((protein_pep_sequence_graph *)x->pg_edge_neighbor)[PEPG_FLAG_NODE_INDEX]) |= (PEPG_FLAG_TERMINAL))

using namespace std;

struct protein_data {
	string prd_id;
	string prd_pep_sequence;
};

typedef vector<protein_data *> data_t;

struct protein_pep_sequence {
	protein_data*	pps_protein;
	int		pps_start, pps_length;
};

/*
 * Graph node to handle the adjacency list
 */
struct protein_pep_sequence_graph {
	protein_pep_sequence_graph* 			pg_edge_neighbor[MAX_EDGES];
	set<protein_pep_sequence*>			pg_pps_list;
};

int ppsg_add_peptide_sequence (protein_pep_sequence_graph *pg, protein_pep_sequence *peptide_seq, int pos);

static int ppsg_num_graph_nodes = 0, ppsg_num_peptides = 0;

int pp_amino_acid_mass[MAX_EDGES] = { 	71, //A 
					99999, //B not there
					103, //C
					115, //D
					129, //E
					147, //F
					57, //G
					137, //H
					113, //I
					99999, //J
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
					99999, //U
					99, //V
					186, //W
					99999, //X
					163, //Y
					99999 //Z 
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
		frag = frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen));
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

void ppsg_add_protein_database_graph (protein_pep_sequence_graph *pg, data_t data)
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
			// cout << "Token is " << token << endl;
			ppsg_add_peptide_sequence(pg, peptide, 0);
			ofs << "Protein ID: " << peptide->pps_protein->prd_id << " Start = " << 
				peptide->pps_start << " End = " << peptide->pps_start + peptide->pps_length - 1 << " " <<  
				(peptide->pps_protein->prd_pep_sequence).substr(peptide->pps_start, peptide->pps_length) << endl;
		}
		split_peptides.clear();
		count++;
	}
}

protein_pep_sequence_graph* ppsg_create_new_node()
{
	protein_pep_sequence_graph *node = new protein_pep_sequence_graph;

	for (int i = 0; i < MAX_EDGES; i++)
		node->pg_edge_neighbor[i] = NULL;

	ppsg_num_graph_nodes++;

	return node;
}

int ppsg_add_peptide_sequence (protein_pep_sequence_graph *pg, protein_pep_sequence *peptide_seq, int pos)
{
	if (pg == NULL) {
		// This should not happen as 
		cout << "There is a BUG. Dumb a$$. Panic" << endl;
		exit(1);
	}

	// get the edge with the mass or amino acid present at pos
	int index = (peptide_seq->pps_protein->prd_pep_sequence).at(peptide_seq->pps_start + pos) - 'A';

	// cout << "Index is " << index << endl;
	protein_pep_sequence_graph *node = pg->pg_edge_neighbor[index];
	if (node == NULL) {
		// if there is none, allocate a node and insert
		protein_pep_sequence_graph *new_node = ppsg_create_new_node();
		node = pg->pg_edge_neighbor[index] = new_node;
	}

	if (pos == peptide_seq->pps_length - 1) {
		// If this is the last one in the sequence, insert protein into the terminal list
		(node->pg_pps_list).insert(peptide_seq);
	} else {
		// recurse
		ppsg_add_peptide_sequence(node, peptide_seq, pos + 1);
	}
}

void ppsg_find_proteins_precursor_mass (protein_pep_sequence_graph *pg, int pmass, int cmass, string pep_seq)
{
	for (int i = 0; i < MAX_EDGES; i++) {
		protein_pep_sequence_graph* pg_tobe = pg->pg_edge_neighbor[i];
		// For all edges, if there is a valid edge
		if (pg_tobe != NULL) {
			// If the mass with that amino acid is less than or equal to, recurse and check terminal
			if ((cmass + pp_amino_acid_mass[i]) < pmass) {

				if ((pg_tobe->pg_pps_list).size() > 0) {
					cout << "The size of the number of peptide sequences in this terminal node are " << (pg_tobe->pg_pps_list).size() << endl;
					set<protein_pep_sequence*>::iterator it;

					char end_char = 'A' + i;
					for (it = (pg_tobe->pg_pps_list).begin(); it != (pg_tobe->pg_pps_list).end(); it++) {
						cout << (*it)->pps_protein->prd_id<< " and the sequence is " << pep_seq + end_char << 
							 " starting from " << (*it)->pps_start << endl;
					}
				}

				char end_char = 'A' + i;
				ppsg_find_proteins_precursor_mass(pg_tobe, pmass, cmass + pp_amino_acid_mass[i], pep_seq + end_char);

			}
		}
	}
}

int main() 
{
	cout << "Size of the node in PepGraph is " << sizeof(protein_pep_sequence_graph) << endl;
	cout << "Size of a peptide sequence is " << sizeof(protein_pep_sequence) << endl;
	// Read file

	// Here is the data we want.
	data_t data;

	// Here is the file containing the data. Read it into data.
	ppsg_read_protein_database("protein_sequence_database.txt", data);

	// Otherwise, list some basic information about the file.
	cout << "Your input CSV file contains " << data.size() << " records.\n";

	/*
	for (protein_data *protein : data) {
		cout << "Protein is " << protein->prd_id << endl;
		cout << "Its amino acid sequence is " << protein->prd_pep_sequence << endl;
	}
	*/

	// For each protein, split the peptide sequence based on the cut amino acid
	// For each cut add to the graph 
	protein_pep_sequence_graph *pg = ppsg_create_new_node();
	ppsg_add_protein_database_graph(pg, data);

	cout << "Number of graph nodes " << ppsg_num_graph_nodes << endl;

	cout << "Number of peptide sequences " << ppsg_num_peptides << endl;

	int precursor_mass;

	while (1) {
		cout << "Enter precursor mass to search: ";
		cin >> precursor_mass;

		ppsg_find_proteins_precursor_mass(pg, precursor_mass, 0, "");
	}
}
