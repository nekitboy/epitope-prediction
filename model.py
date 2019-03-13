import modeller
from modeller.automodel import *


def get_model(first_seq, second_seq):
    env = modeller.environ()
    aln = modeller.alignment(env)

    path = 'test_files/'
    first_file = path + first_seq + '.ali'
    second_file = path + second_seq + '.pdb'

    mdl = modeller.model(env, file=second_file)  # , model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, align_codes=second_seq, atom_files=second_file)
    aln.append(file=first_file, align_codes=first_seq)
    aln.align2d()
    aln.write(file='ans.ali', alignment_format='PIR')

    a = automodel(env, alnfile='ans.ali',
                  knowns=second_seq, sequence=first_seq,
                  assess_methods=(assess.DOPE,
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 6
    a.make()
    max_en = 0
    max_out = None
    for out in a.outputs:
        if out['molpdf'] > max_en:
            max_en = out['molpdf']
            max_out = out

    if max_out is not None:
        return max_out['name']
    else:
        return None


print(get_model("TvLDH", "1bdm"))
