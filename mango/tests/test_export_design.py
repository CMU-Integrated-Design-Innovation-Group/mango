from mango.utils.design_io import export_DNA_design
import dill

sample_result = './test_outputs/ramp_test.aj1'
with open(sample_result, 'rb') as f:
    data = dill.load(f)

# NOTE: If running these tests, you must specify this path:
DAED_PATH = '/Users/kodak/Desktop/PERDIX-Mac/DAEDALUS2'
TALOS_PATH = ''
SEQUENCE_FILE = '/Users/kodak/Desktop/PERDIX-Mac/M13.txt'
export_DNA_design(automated_scaffold_executable_path=DAED_PATH, design=data.design_space, export_path='./test_outputs',
                  savename_no_extension='shape_annealing_test', scaffold_sequence_filepath=SEQUENCE_FILE,
                  ply_or_obj='obj')