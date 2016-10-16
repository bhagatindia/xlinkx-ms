Protein hash database (PHD) is an utility to vend peptides of a particular mass corresponding 
to a protein database. In order to achieve this, it takes protein database, enzyme digestion 
criteria as input and produces a file which records these parameters along with the peptides
corresponding to a particular mass.

This document outlines the usage of protein hash database. In order to use this for retrieving
peptides of a particular mass in a protein database you have to follow the below steps:
1. Initialize the library, and retrieve or construct the protein hash database using 
   phd_retrieve. This will return a phd_file pointer.
2. Using the phd_file pointer, get the peptides corresponding to a particular mass using phd_get

You can pass appropriate flags to phd_retrieve for it to create a new hash database to to use an
exsting one. A new file is created when:
* a file with the specified hash database file name is not present
* enzyme cut paramaeters doesnt match
* the protein database file name doesnt match
* the protein database file digest doesnt match the digest the digest present in phd
