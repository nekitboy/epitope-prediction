from Bio.Blast import NCBIWWW
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PPBuilder, PDBList
from Bio.Blast import NCBIXML

pdb_parser = PDBParser(PERMISSIVE=True)
pp_builder = PPBuilder()
pdbl = PDBList()


def get_structure(file, p_id):
    """
    Get pdb structure
    :param file: path to the pdb file
    :param p_id: id for the structure
    :return pdb structure object:
    """
    return pdb_parser.get_structure(p_id, file)


def get_sequences(protein):
    """
    Get list of protein sequences objects from pdb structure (model, chain, residue, atom)
    :param protein:
    :return:
    """
    seqs = pp_builder.build_peptides(protein)

    return [seq.get_sequence() for seq in seqs]


def get_sequences_from_file(file, p_id):
    """
    Get list of protein sequences objects from pdb file
    :param file: path to the pdb file
    :param p_id: id for the structure
    :return: list of protein sequences:
    """
    return get_sequences(get_structure(file, p_id))


def blast(seqs):
    """
    Blast sequences
    :param seqs: list of seqs objects
    :return: list of blast records
    """
    records = []

    for seq in seqs:
        result = NCBIWWW.qblast("blastp", "pdb", seq)
        result_record = list(NCBIXML.parse(result))
        if len(result_record) > 1:
            raise Exception('Length of blast samples bigger than 1')
        record = result_record[0]
        records.append(record)
    return records


def download_pdb(pdb_code, d_path=None):
    """
    Download pdb file for pdb id code
    :param pdb_code: pdb id code
    :param d_path: path to download
    :return: path to downloaded pdb file
    """
    return pdbl.retrieve_pdb_file(pdb_code, pdir=d_path, file_format='pdb')


def get_chain(structure, chain_key):
    """
    Get chain from pdb structure corresponds to chain_key
    :param structure: pdb structure object
    :param chain_key: chain key id ('A', 'B')
    :return: chain object
    """
    models = structure.get_models()
    if len(models) > 1:
        raise Exception('Amount of models in pdb structure more than 1')
    return models[0][chain_key]