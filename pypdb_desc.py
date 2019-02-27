import pypdb

# https://github.com/williamgilpin/pypdb/blob/master/pypdb/pypdb.py

# Return a list of all PDB entries
a = pypdb.get_all()[:10]

# Return info about protein
print(pypdb.get_all_info(a[0]))

# Return ligands of given PDB ID
print(pypdb.get_ligands(a[0]))

# Return BLAST search results for a given PDB ID
print(pypdb.get_blast(a[0]))

# pypdb.get_pdb_file - 404 Error

