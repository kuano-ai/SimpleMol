# Simple Mols Programming Exercise

In this folder are files containing python code to manipulate molecules. 


## Intro

Two of the most common ways of storing molecules as data are smiles strings and sdf files. A [smiles string](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) is a series of characters used to represent a molecule, follow the link for more information. A smiles file ("*.smi") contains a list of smiles strings. An [sdf file](https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics/02%3A_Representing_Small_Molecules_on_Computers/2.05%3A_Structural_Data_Files) is a data format to store molecules with coordinates, to record a 2D or 3D configuration of a molecule, as well as it's basic structural information.

## Files in this folder:

- The file `simple_mol.py` provides a class which holds molecules in basic python datatypes using numpy.
- The file `example_script.py` shows how to use the SimpleMol class and read functions to read in molecules and display their properties


## Exercise 1
- Write a python script which reads the `basic_mols.smi` file from the molecules folder, creates 2 SimpleMol class objects from it, calculate the number of atoms in the molecule and print the result.

## Ideas for other exercises
- Find an sdf file from someone, and recreate the example_script but for this new file
- Add code to the example script to save the molecules as smiles strings in a new file
- Write a function which replaces oxygen atoms with nitrogen atoms
- Write a function which finds carbon atoms with single bonds