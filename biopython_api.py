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


def get_pdbs_from_blast_result(records, exclude_list):
    """
    Get pdb id's with chain id's with which blast input well relates with blast output
    :param records: list of blast records
    :param exclude_list: list of pdb id's which should be exclude from result list
    :return:
    """

    pdb_chain_list = []

    for record in records:
        pdb_chain_dict = {}
        for alignment in record.alignments:

            pdb_chain_tuples = [(i.split('|')[3], i.split('|')[4].split(' ')[0]) for i in alignment.title.split('>')]
            # Example: ('1AHW', 'A'), ('1AHW', 'D'), ('1FGN', 'L')

            # Check exclude list
            if sum([1 if (t[0] in exclude_list) else 0 for t in pdb_chain_tuples]) > 0:
                continue

            # Assumption that we get only one hsps
            if len(alignment.hsps) > 1:
                raise Exception("hsp objects more than 1 error in code. Change hsp handler")
            # HSP is object of hit
            hsp = alignment.hsps[0]
            for pdb_chain in pdb_chain_tuples:
                if pdb_chain[0] not in pdb_chain_dict:
                    pdb_chain_dict[pdb_chain[0]] = dict()
                pdb_chain_dict[pdb_chain[0]][pdb_chain[1]] = hsp
                """
                HSP Members:
                    - score           BLAST score of hit.  (float)
                    - bits            Number of bits for that score.  (float)
                    - expect          Expect value.  (float)
                    - num_alignments  Number of alignments for same subject.  (int)
                    - identities      Number of identities (int) if using the XML parser.
                      Tuple of number of identities/total aligned (int, int)
                      if using the (obsolete) plain text parser.
                    - positives       Number of positives (int) if using the XML parser.
                      Tuple of number of positives/total aligned (int, int)
                      if using the (obsolete) plain text parser.
                    - gaps            Number of gaps (int) if using the XML parser.
                      Tuple of number of gaps/total aligned (int, int) if
                      using the (obsolete) plain text parser.
                    - align_length    Length of the alignment. (int)
                    - strand          Tuple of (query, target) strand.
                    - frame           Tuple of 1 or 2 frame shifts, depending on the flavor.
            
                    - query           The query sequence.
                    - query_start     The start residue for the query sequence.  (1-based)
                    - query_end       The end residue for the query sequence.  (1-based)
                    - match           The match sequence.
                    - sbjct           The sbjct sequence.
                    - sbjct_start     The start residue for the sbjct sequence.  (1-based)
                    - sbjct_end       The end residue for the sbjct sequence.  (1-based)
                """

        pdb_chain_list.append(pdb_chain_dict)

    result_pdb_with_chains = intersect_chains_list(pdb_chain_list)

    return result_pdb_with_chains


def intersect_chains_list(pdb_chain_list):
    result_pdb_with_chains = {}
    for pdb_id in pdb_chain_list[0]:
        is_good = True
        chains = [pdb_chain_list[0][pdb_id]]

        for i in range(1, len(pdb_chain_list)):
            if pdb_id not in pdb_chain_list[i]:
                is_good = False
            else:
                chains.append(pdb_chain_list[i][pdb_id])

        if is_good:
            result_pdb_with_chains[pdb_id] = chains

    return result_pdb_with_chains
