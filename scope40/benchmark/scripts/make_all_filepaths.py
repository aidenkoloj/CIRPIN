import os
import sys


directory = sys.argv[1]

cp_filepaths = sys.argv[2]

def make_all_v_all_filepath(directory):
    with open(os.path.join(directory, cp_filepaths), 'a') as fp:
        for filename in os.listdir(directory):
            if filename.endswith('.pdb'):
                print("Added", filename, "to path")
                fp.write(filename + ' ' + filename[:-4]+ '\n')

make_all_v_all_filepath(directory)

