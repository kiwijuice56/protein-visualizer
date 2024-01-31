import os
import subprocess

# Any folder with a collection of .pdb or .cif files
directory = "test_proteins/"

for file in os.listdir(directory):
    path = os.path.join(directory, file)
    if path.endswith(".cif"):
        subprocess.call(['python', 'main.py', path])
