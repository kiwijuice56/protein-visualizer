import os
import subprocess

# Any folder with a collection of .pdb files
directory = "proteins/"

for file in os.listdir(directory):
    path = os.path.join(directory, file)
    if path.endswith(".pdb"):
        subprocess.call(['python', 'main.py', path])
