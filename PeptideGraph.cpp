#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>

#define MAX_EDGES 26

// B is not present as an amino acid. So, reusing it to store protein name
#define PEPG_NAME_NODE_INDEX 1
// Z is not present as an amino acid. So, reusing it to store flags
#define PEPG_FLAG_NODE_INDEX 25

#define PEPG_FLAG_TERMINAL 0x01

#define IS_PEPG_NODE_TERMINAL(x) ((((pep_graph *)x->pg_edge_neighbor)[PEPG_FLAG_NODE_INDEX]) & (PEPG_FLAG_TERMINAL))
#define SET_PEPG_NODE_TERMINAL(x) ((((pep_graph *)x->pg_edge_neighbor)[PEPG_FLAG_NODE_INDEX]) |= (PEPG_FLAG_TERMINAL))

using namespace std;

struct protein_data {
	string pr_id;
	string pr_pep_sequence;
};

typedef vector<protein_data *> data_t;

/*
 * Graph node to handle the adjacency list
 */
struct pep_graph {
	pep_graph 		*pg_edge_neighbor[MAX_EDGES];
	vector<protein_data*>	pg_protein_list;
};

int ppsg_add_peptide_sequence (pep_graph *pg, protein_data *protein, std::string &peptide_seq, int pos);

static int ppsg_num_graph_nodes = 0;

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

void ppsg_split(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
	/* Store the original string in the array, so we can loop the rest
	 * of the algorithm. */
	tokens.push_back(str);

	// Store the split index in a 'size_t' (unsigned integer) type.
	size_t splitAt;
	// Store the size of what we're splicing out.
	size_t splitLen = splitBy.size();
	// Create a string for temporarily storing the fragment we're processing.
	std::string frag;
	// Loop infinitely - break is internal.
	while(true)
	{
		/* Store the last string in the vector, which is the only logical
		 * candidate for processing. */
		frag = tokens.back();
		/* The index where the split is. */
		splitAt = frag.find(splitBy);
		// If we didn't find a new split point...
		if(splitAt == string::npos)
		{
			// Break the loop and (implicitly) return.
			break;
		}
		/* Put everything from the left side of the split where the string
		 * being processed used to be. */
		tokens.back() = frag.substr(0, splitAt);
		/* Push everything from the right side of the split to the next empty
		 * index in the vector. */
		tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
	}
}

void ppsg_read_protein_database(string file, data_t& data)
{
	// For every record we can read from the file, append it to our resulting data
	protein_data *record, *old_record;

	std::ifstream fstream(file.c_str());
	std::string line, pr_id, pr_seq;	

	string split_char = " ";
	vector<string> split_tokens;

	while (std::getline(fstream, line)) {
		if (line.at(0) == '>') {
			// terminate old sequence if it was started 
			record = new protein_data;
			data.push_back(record);
			ppsg_split(line, split_char, split_tokens);
			record->pr_id = split_tokens[0];
		} else {
			// append or start a new one
			(record->pr_pep_sequence).append(line);
		}
	}

	return;
}

void ppsg_add_protein_database_graph (pep_graph *pg, data_t data)
{
	string split_char = "K";
	vector<string> split_tokens;

	int count = 0;
	for (protein_data *protein : data) {
		cout << "Protein number is " << count++ << " and it is " << protein->pr_pep_sequence << endl;
		ppsg_split(protein->pr_pep_sequence, split_char, split_tokens);
		for (string token : split_tokens) {
			// cout << "Token is " << token << endl;
			if (token.length() == 0) continue;
			ppsg_add_peptide_sequence(pg, protein, token, 0);
		}
		split_tokens.clear();
	}
}

int ppsg_add_protein_sequence (pep_graph pg, protein_data protein)
{
}

pep_graph* ppsg_create_new_node()
{
	pep_graph *node = new pep_graph;

	for (int i = 0; i < MAX_EDGES; i++)
		node->pg_edge_neighbor[i] = NULL;

	ppsg_num_graph_nodes++;

	return node;
}

int ppsg_add_peptide_sequence (pep_graph *pg, protein_data *protein, std::string &peptide_seq, int pos)
{
	if (pg == NULL) {
		// This should not happen as 
		cout << "There is a BUG. Dumb a$$. Panic" << endl;
		exit(1);
	}

	// get the edge with the mass or amino acid present at pos
	int index = peptide_seq.at(pos) - 'A';

	// cout << "Index is " << index << endl;
	pep_graph *node = pg->pg_edge_neighbor[index];
	if (node == NULL) {
		// if there is none, allocate a node and insert
		pep_graph *new_node = ppsg_create_new_node();
		node = pg->pg_edge_neighbor[index] = new_node;
	}

	if (pos == peptide_seq.length() - 1) {
		// If this is the last one in the sequence, 
		// mark and insert protein into the terminal list
		// SET_PEPG_NODE_TERMINAL(node);
		(node->pg_protein_list).push_back(protein);
	} else {
		// recurse
		ppsg_add_peptide_sequence(node, protein, peptide_seq, pos + 1);
	}
}

void ppsg_find_proteins_precursor_mass (pep_graph *pg, int pmass, int cmass, string pep_seq)
{
	for (int i = 0; i < MAX_EDGES; i++) {
		// For all edges, if there is a valid edge
		if ((pg->pg_edge_neighbor)[i] != NULL) {
			// If the mass with that amino acid is less than or equal to, recurse and check terminal
			if ((cmass + pp_amino_acid_mass[i]) < pmass) {
				char end_char = 'A' + pp_amino_acid_mass[i];
				ppsg_find_proteins_precursor_mass((pg->pg_edge_neighbor)[i], pmass, cmass + pp_amino_acid_mass[i], pep_seq + end_char);

				if ((pg->pg_protein_list).size() > 0) {
					for (protein_data *protein : pg->pg_protein_list) {
						cout << protein->pr_id << endl;
					}
				}
			}
		}
	}
}

int main() 
{
	cout << "Size of the node PepGraph is " << sizeof(pep_graph) << endl;
	// Read file

	// Here is the data we want.
	data_t data;

	// Here is the file containing the data. Read it into data.
	ppsg_read_protein_database("protein_sequence_database.txt", data);

	// Otherwise, list some basic information about the file.
	cout << "Your input CSV file contains " << data.size() << " records.\n";

	/*
	for (protein_data *protein : data) {
		cout << "Protein is " << protein->pr_id << endl;
		cout << "Its amino acid sequence is " << protein->pr_pep_sequence << endl;
	}
	*/

	// For each protein, split the peptide sequence based on the cut amino acid
	// For each cut add to the graph 
	pep_graph *pg = ppsg_create_new_node();
	ppsg_add_protein_database_graph(pg, data);

	cout << "Number of graph nodes " << ppsg_num_graph_nodes << endl;

	int precursor_mass;

	while (1) {
		cin >> precursor_mass;

		ppsg_find_proteins_precursor_mass(pg, precursor_mass, 0, "");
	}
}
