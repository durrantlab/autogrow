"""AutoGrow 3.1.2 is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AutoGrow is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2015 Jacob D. Durrant. If you have any questions, comments, or
suggestions, please don't hesitate to contact me at jdurrant [at] ucsd [dot] edu.

The latest version of AutoGrow can be downloaded from 
http://autogrow.ucsd.edu.

If you use AutoGrow in your work, please cite:

J. D. Durrant, R. E. Amaro, J. A. McCammon. Autogrow: A Novel Algorithm for Protein Inhibitor Design. Chem. Biol. Drug Des. 2009, 73, 168-178.

J. D. Durrant, S. Lindert, J. A. McCammon. Autogrow 3.0: An Improved Algorithm for Chemically Tractable, Semi-Automated Protein Inhibitor Design. J. Mol. Graphics Modell. 2013, 44, 104-112.

"""

def define_defaults(): # this function is first to make it easy for the user to modify the default parameters (see, especially, the boxes below)
	"""Sets the command-line parameters to their default values."""
	
	vars = {}
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	
	
	################################################################# FILE-LOCATION VARIABLES ###################################################################################
	vars['mgltools_directory'] = ""							   # Example: vars['mgltools_directory'] = "/home/myname/MGLTools-1.5.4/"
	vars['openbabel_bin_directory'] = ""						   # Example: vars['openbabel_bin_directory'] = "/home/myname/openbabel-2.2.0/bin/"
	vars['vina_executable'] = ""							   # Example: vars['vina_executable'] ="/home/myname/autodock_vina_1_1_2_linux_x86/bin/vina"
	vars['autoclickchem_script'] = script_dir + os.sep + "support" \
			+ os.sep + "autoclickchem" + os.sep + "autoclickchem.py"	   # Example: vars['autoclickchem_script'] ="/home/myname/autoclickchem/autoclickchem.py"
	vars['nn1_script'] = script_dir + os.sep + "support" + os.sep + "NNScore" \
			+ os.sep + "NNScore_1.0" + os.sep + "NNScore.py"		   # Example: vars['nn1_script'] ="/home/myname/NNScore/NNScore_1.0/NNScore.py"
	vars['nn2_script'] = script_dir + os.sep + "support" + os.sep + "NNScore" \
			+ os.sep + "NNScore_2.01" + os.sep + "NNScore2.01.py"		   # Example: vars['nn2_script'] ="/home/myname/NNScore/NNScore_2.01/NNScore2.01.py"
	#############################################################################################################################################################################
	

	############################################################# OPTIONAL FILE-LOCATION VARIABLES ################################################################
	#					     (RECOMMEND SETTING TO "" SO AUTOGROW CAN AUTOLOCATE THESE FILES)												 #
	###############################################################################################################################################################
	vars['babel_executable'] = "" 	   # vars['babel_executable'] = "/home/myname/openbabel-2.2.0/bin/babel"
	vars['obprop_executable'] = "" 	   # vars['obprop_executable'] = "/home/myname/openbabel-2.2.0/bin/obprop"
	vars['prepare_ligand4.py'] = ""    # vars['prepare_ligand4.py'] = "/home/myname/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
	vars['prepare_receptor4.py'] = ""  # vars['prepare_receptor4.py'] = "/home/myname/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
	vars['mgl_python'] = "" 	   # vars['mgl_python'] = "/home/myname/MGLTools-1.5.4/bin/pythonsh"
	###############################################################################################################################################################


	vars['number_of_mutants'] = 40
	vars['number_of_crossovers_first_generation'] = 0
	vars['number_of_mutants_first_generation'] = 0
	vars['number_of_crossovers'] = 0
	vars['top_ones_to_advance_to_next_generation'] = 10
	vars['num_generations'] = 10
	vars['max_seconds_per_generation'] = 1000000 # about a week and a half
	
	vars['center_x'] = 0.0
	vars['center_y'] = 0.0
	vars['center_z'] = 0.0
	
	vars['size_x'] = 30.0
	vars['size_y'] = 30.0
	vars['size_z'] = 30.0
	
	vars['num_processors'] = 0 # with this setting, the program will use all processors available
	
	vars['directory_of_source_compounds'] = "." + os.sep
	vars['directory_of_fragments'] = os.path.realpath(__file__) + os.sep + "fragments" + os.sep + "MW_250" + os.sep + ""
	vars['filename_of_receptor'] = "" 
	vars['output_dir'] = "." + os.sep + "output" + os.sep
	
	vars['use_lipinski_filter'] = "TRUE"
	vars['use_strict_lipinski_filter'] = "FALSE"
	vars['use_ghose_filter'] = "FALSE"
	
	vars['allow_modification_without_frag_addition'] = "FALSE"

	vars['maintain_core'] = "FALSE"
	vars['minimum_core_atoms_required'] = 1
	
	vars['score_by_ligand_efficiency'] = "FALSE"
	
	vars['scoring_function'] = "VINA"
	
	vars['additional_autoclickchem_parameters'] = ""
	
	return vars

# import modules
import sys
import os
import glob
import random
import operator
import shutil
import scipy.optimize # used for its nice minimization functions
import numpy
import math
import time
import cPickle
import subprocess
import datetime
import textwrap
import copy
import pymolecule
import multiprocessing
import string

class Molecule2(pymolecule.Molecule):
	"""This class extends the pymolecule.Molecule class with additional functions needed by AutoGrow"""
	
	def rotate_molecule_around_a_fixed_center(self, theta_x, theta_y, theta_z, x, y, z):
		"""Rotates this pymolecule.Molecule about a fixed center.
		
		Arguments:
		theta_x -- Rotation about x axis, in radians (float).
		theta_y -- Rotation about y axis, in radians (float).
		theta_z -- Rotation about z axis, in radians (float).
		x -- x coordinate of fixed pivot point (float)
		y -- y coordinate of fixed pivot point (float)
		z -- z coordinate of fixed pivot point (float)
		
		"""
		
		# now, move the pdb to the given rotation point
		self.translate_molecule(pymolecule.Point(-x, -y, -z))

		# now rotate
		sinx = math.sin(theta_x)
		siny = math.sin(theta_y)
		sinz = math.sin(theta_z)
		cosx = math.cos(theta_x)
		cosy = math.cos(theta_y)
		cosz = math.cos(theta_z)

		cosy_cosz = cosy * cosz
		sinx_siny_cosz_plus_cosx_sinz = sinx * siny * cosz + cosx * sinz
		sinx_sinz_minus_cosx_siny_cosz = sinx * sinz - cosx * siny * cosz
		cosy_sinz = cosy * sinz
		cosx_cosz_minus_sinx_siny_sinz = cosx * cosz - sinx * siny * sinz
		cosx_siny_sinz_plus_sinx_cosz = cosx * siny * sinz + sinx * cosz
		sinx_cosy = sinx * cosy
		cosx_cosy = cosx * cosy

		for atomindex in self.all_atoms: 
			vector = self.all_atoms[atomindex].coordinates
			
			new_x = vector.x * cosy_cosz + vector.y * sinx_siny_cosz_plus_cosx_sinz + vector.z * sinx_sinz_minus_cosx_siny_cosz
			new_y = -vector.x * cosy_sinz + vector.y * cosx_cosz_minus_sinx_siny_sinz + vector.z * cosx_siny_sinz_plus_sinx_cosz
			new_z = vector.x * siny - vector.y * sinx_cosy + vector.z * cosx_cosy
	
			self.all_atoms[atomindex].coordinates = pymolecule.Point(new_x, new_y, new_z)

		# now move it back from the origin
		self.translate_molecule(pymolecule.Point(x, y, z))

	def rotate_molecule(self, theta_x, theta_y, theta_z):
		"""Rotates this pymolecule.Molecule about it's center.
		
		Arguments:
		theta_x -- Rotation about x axis, in radians (float).
		theta_y -- Rotation about y axis, in radians (float).
		theta_z -- Rotation about z axis, in radians (float).
		
		"""

		# first, identify the geometric center
		x = 0.0
		y = 0.0
		z = 0.0
		count = 0
		for index in self.all_atoms:
			if self.all_atoms[index].element != "H":
				count = count + 1
				x = x + self.all_atoms[index].coordinates.x
				y = y + self.all_atoms[index].coordinates.y
				z = z + self.all_atoms[index].coordinates.z
		x = x / count
		y = y / count
		z = z / count

		# now, move the pdb to the origin
		self.translate_molecule(pymolecule.Point(-x, -y, -z))

		# now rotate
		sinx = math.sin(theta_x)
		siny = math.sin(theta_y)
		sinz = math.sin(theta_z)
		cosx = math.cos(theta_x)
		cosy = math.cos(theta_y)
		cosz = math.cos(theta_z)

		cosy_cosz = cosy * cosz
		sinx_siny_cosz_plus_cosx_sinz = sinx * siny * cosz + cosx * sinz
		sinx_sinz_minus_cosx_siny_cosz = sinx * sinz - cosx * siny * cosz
		cosy_sinz = cosy * sinz
		cosx_cosz_minus_sinx_siny_sinz = cosx * cosz - sinx * siny * sinz
		cosx_siny_sinz_plus_sinx_cosz = cosx * siny * sinz + sinx * cosz
		sinx_cosy = sinx * cosy
		cosx_cosy = cosx * cosy

		for atomindex in self.all_atoms: 
			vector = self.all_atoms[atomindex].coordinates
			
			new_x = vector.x * cosy_cosz + vector.y * sinx_siny_cosz_plus_cosx_sinz + vector.z * sinx_sinz_minus_cosx_siny_cosz
			new_y = -vector.x * cosy_sinz + vector.y * cosx_cosz_minus_sinx_siny_sinz + vector.z * cosx_siny_sinz_plus_sinx_cosz
			new_z = vector.x * siny - vector.y * sinx_cosy + vector.z * cosx_cosy
	
			self.all_atoms[atomindex].coordinates = pymolecule.Point(new_x, new_y, new_z)

		# now move it back from the origin
		self.translate_molecule(pymolecule.Point(x, y, z))

	def get_center(self):
		"""Calculates the geometric center of this pymolecule.Molecule
		
		Returns:
		A list of floats, representing the x, y, and z coordinates of the geometric center
		
		"""

		center = [0.0, 0.0, 0.0]
		x = 0.0
		y = 0.0
		z = 0.0
		count = 0
		for index in self.all_atoms:
			if self.all_atoms[index].element != "H":
				count = count + 1
				x = x + self.all_atoms[index].coordinates.x
				y = y + self.all_atoms[index].coordinates.y
				z = z + self.all_atoms[index].coordinates.z
		center[0] = x / count
		center[1] = y / count
		center[2] = z / count
		return center

	def max_index(self):
		"""Identifies the maximum index of any atom in this pymolecule.Molecule object
		
		Returns:
		An integer, the maximum index of any atom in this pymolecule.Molecle object.
		
		"""

		maxindex = -99999999
		for index in self.all_atoms.keys():
			if index > maxindex: maxindex = index
		return maxindex

	def save_pdb2(self, filename, remarks = []):
		"""Saves the python.Molecule object to a PDB file.
		
		Arguments:
		filename -- The name of the file to save (String)
		remarks -- Any remarks to include in the file header (a list of strings)
		
		"""

		if len(self.all_atoms) > 0: # so the pdb is not empty (if it is empty, don't save)

			file = open(filename,"w")

			for remark in remarks:
				file.write("REMARK " + remark + "\n")

			# write coordinates
			for atomindex in self.all_atoms:
				file.write(self.all_atoms[atomindex].create_pdb_line(atomindex) + "\n")

			file.close()

##########################
# multithreading classes #
##########################

class multi_threading(): 
	"""Launch jobs on multiple processors"""
	
	def __init__(self, inputs, task_class_name, variables_to_pass):
		"""Launches jobs on multiple processors
		
		Arguments:
		inputs -- the data to be processed, in a list
		task_class_name -- the name of the class charged with managing the jobs assigned to each processor, a string
		variables_to_pass -- additional variables (usually command-line parameters) to be passed to each processor, a dictionary

		"""
		
		if len(inputs) == 0: return

		global vars
		num_processors = vars['num_processors']
		
		if num_processors == 1: # so it's running on just one processor, don't even try to parallelize it
			run_class = task_class_name()
			for inp in inputs: run_class.value_func(inp, variables_to_pass)

		else: # so run in parallel
			# first, if num_processors <= 0, determine the number of processors to use programatically
			if num_processors <= 0: num_processors = multiprocessing.cpu_count()
	
			# reduce the number of processors if too many have been specified
			if len(inputs) < num_processors: num_processors = len(inputs)
	
			# if the appropriate filename is present, write the contents for record-keeping purposes
			if 'log_filename' in variables_to_pass.keys():
				f = open(variables_to_pass['log_filename'],'w')
				for an_input in inputs: f.write(an_input + "\n")
				f.close()
	
			# now, divide the inputs into the appropriate number of processors
			inputs_divided = {}
			for t in range(num_processors): inputs_divided[t] = []
	
			for t in range(0, len(inputs), num_processors):
				for t2 in range(num_processors):
					index = t + t2
					if index < len(inputs): inputs_divided[t2].append(inputs[index])
	
			# now, run each division on its own processor
			running = multiprocessing.Value('i', num_processors)
			mutex = multiprocessing.Lock()
	
			arrays = []
			threads = []
			for i in range(num_processors):
				threads.append(task_class_name())
				arrays.append(multiprocessing.Array('i',[0, 1]))
	
			processes = []
			for i in range(num_processors):
				p = multiprocessing.Process(target=threads[i].runit, args=(running, mutex, inputs_divided[i], variables_to_pass))
				p.start()
				processes.append(p)
	
			while running.value > 0: is_running = 0 # wait for everything to finish
			
class general_task: # other, more specific classes with inherit this one
	"""Run jobs on a single processor"""
	
	def runit(self, running, mutex, items, variables_to_pass):
		"""Run jobs on a single processor
		
		Arguments:
		running -- a multiprocessing.Value() object
		mutex -- a multiprocessing.Lock() object
		items -- the data to be processed, in a list
		variables_to_pass -- additional variables (usually command-line parameters) to be passed to each processor, a dictionary

		"""

		for item in items: self.value_func(item, variables_to_pass)
		mutex.acquire()
		running.value -= 1
		mutex.release()

class execute_command(general_task):
	"""Run shell commands on a single processor"""
	
	def value_func(self, command, variables_to_pass):
		"""Run a single shell command
		
		Arguments:
		command -- the command to run, a string
		variables_to_pass -- additional variables (usually command-line parameters)

		"""

		log("Executing: " + command)
		proc = subprocess.Popen(command, shell=True)
		for idx in range(variables_to_pass['seconds_per_job']):
			if proc.poll() is not None: # so it is not running
				break
			time.sleep(1)
		if proc.poll() is None: # so it is still running
			proc.terminate()
			log("Had to terminate a job early: " + command)

class run_autoclick_multithread(general_task):
	"""Run autoclickchem on a single processor"""
	
	def value_func(self, command, variables_to_pass):
		"""Run a autoclickchem in silico reaction
		
		Arguments:
		command -- parameter not actually used, but needed for compatibility with the multi_threading class
		variables_to_pass -- additional variables (usually command-line parameters)

		"""

		global vars

		check_time_limit_reached() # has the program been running too long? If so, stop the program.

		# pick a ligand
		ligand = random.choice(variables_to_pass['ligands'])

		# pick a fragment
		frag_dir = random.choice(variables_to_pass['frag_dirs']) # pick the kind of fragment
		fragment = random.choice(variables_to_pass['fragments'][frag_dir]) # pick the specific fragment

		# run autoclick chem to generate new ligands		
		args = ['autoclickchem.py', '-reactants1', ligand, '-reactants2', fragment, '-output_dir', variables_to_pass['directory']] # so all possible reactions will be performed, then you'll prune down later
		
		if vars['allow_modification_without_frag_addition'] == "FALSE": # so prevent autoclickchem from running the reactions that just modify the ligand without adding fragments
			args.append('-alkene_to_epoxide')
			args.append('-halide_to_cyanide')
			args.append('-alcohol_to_cyanide')
			args.append('-carboxylate_to_cyanide')
			args.append('-acyl_halide_to_cyanide')
			args.append('-acid_anhydride_to_cyanide')
			args.append('-halide_to_azide')
			args.append('-alcohol_to_azide')
			args.append('-carboxylate_to_azide')
			args.append('-acyl_halide_to_azide')
			args.append('-acid_anhydride_to_azide')
			args.append('-amine_to_azide')
			args.append('-amine_to_isocyanate')
			args.append('-amine_to_isothiocyanate')
			args.append('-azide_to_amine')

		# now see if the user has passed any additional AutoClickChem parameters
		params = vars['additional_autoclickchem_parameters']
		while "  " in params: params = params.replace('  ',' ')
		params = params.split(' ')
		args.extend(params)
		
		new_ligs_filenames = run_autoclick_chem(args)

		for lig in new_ligs_filenames: log("Making mutant: " + os.path.basename(ligand) + " + " + os.path.basename(fragment) + " => " + os.path.basename(lig))


class run_ligmerge_multithread(general_task):
	"""Run ligmerge on a single processor"""
	
	def value_func(self, command, variables_to_pass):
		"""Run a autoclickchem in silico reaction
		
		Arguments:
		command -- parameter not actually used, but needed for compatibility with the multi_threading class
		variables_to_pass -- additional variables (usually command-line parameters)

		"""

		global vars

		ligands = variables_to_pass['ligands']
		directory = variables_to_pass['directory']
		
		check_time_limit_reached()

		# pick two ligands
		count2 = 0
		while True:
			count2 = count2 + 1
			ligand1 = random.choice(ligands)
			ligand2 = random.choice(ligands)
			if ligand1 != ligand2: break # you've got different names
			if count2 > 10000:
				log("ERROR: I've tried 10,000 times to identify two unique ligands to serve as \"parents\" in the production of a crossover ligand. Aborting program...")
				sys.exit(0)
			
		pdb1 = Molecule2()
		pdb1.load_pdb(ligand1)
		pdb2 = Molecule2()
		pdb2.load_pdb(ligand2)

		aLigMerge = LigMerge()
		merged,message = aLigMerge.main(pdb1, pdb2, 'false', 'false', 'false', 20) # so if it doesn't generate a hybrid molecule in 20 seconds, abort.

		if merged is not None:
			if len(merged) > 0: merged = merged[0]
			else: merged is None

		if merged is not None: 
			newindex = random.randint(0,1000000)
			filename = directory + str(newindex) + ".crossover.pdb" 

			try:
				merged.save_pdb2(filename, ['Source Files:',ligand1, ligand2])
			except AttributeError:
				nthing = ""

			log("Making crossover: " + os.path.basename(ligand1) + " + " + os.path.basename(ligand2) + " => " + os.path.basename(filename))

class druglike_multithread(general_task):
	"""Evaluate druglike properties of moleules on a single processor"""
	
	def value_func(self, filename, variables_to_pass):
		"""Run a autoclickchem in silico reaction
		
		Arguments:
		filename -- the filename of the compounds, a string
		variables_to_pass -- additional variables (usually command-line parameters)

		"""

		check_time_limit_reached()

		if check_drug_likeness(filename) == False:
			log("Not druglike: deleting ligand " + filename)
			delete_all_associated_files(filename)
		else:
			log("Passes drug-like filter (if any): " + filename)


####################
# ligmerge classes #
####################

class LigMerge:
	"""The AutoGrow crossover operator, adapted from a beta version of LigMerge."""

	def __transform_molecule(self,x,*args):
		"""Transforms (translates and rotates) the atoms of a pymolecule.Molecule object.
		
		Arguments:
		x -- A numpy array of the parameters to be optimized (numpy.array).
		args -- A tuple containing the fixed parameters for the optimizations
		
		"""
		
		# rename variables to make things easier.
		trans_x = x[0]
		trans_y = x[1]
		trans_z = x[2]

		theta_x = x[3]
		theta_y = x[4]
		theta_z = x[5]

		pdb_to_move = args[0]
		if len(args) > 1:
			center = args[1]
	
		pdb_to_move.undo() # we move that PDB back to its default value.
		pdb_to_move.translate_molecule(pymolecule.Point(trans_x, trans_y, trans_z)) # we then translate the smaller molecule

		if len(args) > 1:
			pdb_to_move.rotate_molecule_around_a_fixed_center(theta_x, theta_y, theta_z, center[0], center[1], center[2]) # we then rotate it
		else:
			pdb_to_move.rotate_molecule(theta_x, theta_y, theta_z) # we then rotate it
		
	def __calculate_rmsd_assignments(self,x,*args):
		"""Aligns and calculates the RMSD between two pymolecule.Molecule objects.
		
		Arguments:
		x -- A numpy array of the parameters to be optimized (numpy.array).
		args -- A tuple containing the fixed parameters for the optimizations
		
		Returns:
		A float, the RMSD between the two ligands.
		
		"""

		# rename variables to make things easier.
		trans_x = x[0]
		trans_y = x[1]
		trans_z = x[2]
	
		theta_x = x[3]
		theta_y = x[4]
		theta_z = x[5]

		pdb_to_move = args[0]
		pdb_orig = args[1]
		assignment = args[2]
		
		# the smaller one is the one that will move. 
		pdb_to_move.undo() # we move that PDB back to its default value.
		delta = pymolecule.Point(trans_x, trans_y, trans_z)
		pdb_to_move.translate_molecule(delta) # we then translate the smaller molecule
	
		pdb_to_move.rotate_molecule(theta_x, theta_y, theta_z) # we then rotate it

		# calculate the RMSD (between the atoms tethered in the assignment)
		total_dist_squared = 0.0 # sum the square distances between atoms, the sum is initially 0
	
	
		count = 0 # to keep track of the number of atom pairs found
		for index in range(len(assignment[0])): 
			to_move_atom = pdb_to_move.all_atoms[assignment[0][index]] # define the atom in pdb_to_move corresponding to that index
			orig_atom = pdb_orig.all_atoms[assignment[1][index]] # define that atom in pdb_orig corresponding to that index
			dist = to_move_atom.coordinates.distance_to_another_point(orig_atom.coordinates) # calculate the distance between the atom from the smaller molecule and the one from the larger
			total_dist_squared = total_dist_squared + math.pow(dist,2) # add the distance squared to the growing sum
			count = count + 1 
	
		if count > 0:
			rmsd = math.pow(total_dist_squared / count,0.5) # here, calculate the rmsd
		else:
			rmsd = 99999999.9
	
		return rmsd # return
	
	def __add_handles(self,molecule_object): # molecule_object is the pdb that the handles will be added to
		"""Adds "handles" to a pymolecule.Molecule object.
		
		Arguments:
		molecule_object -- The molecular model to which handles will be added (pymolecule.Molecule).
		
		Returns:
		A pymolecule.Molecule object, with handles added.
		
		"""
		
		global common_pdb # this is the pdb with common (superimposed) atoms, duplicated removed.
	
		already_used = [] # a list of atoms from the common that have been identified as handles. This is just for a check.
		toadd = [] # a list of the handle atoms that will be added to the molecule_object
		for atom_index in molecule_object.all_atoms: # go through each atom in molecule_object
			atom = molecule_object.all_atoms[atom_index]
			for common_index in common_pdb.all_atoms: # go through each atom in common_pdb
				common_atom = common_pdb.all_atoms[common_index]
				dist = atom.coordinates.distance_to_another_point(common_atom.coordinates) # get the distance between these two
				
				ideal_bondlength = molecule_object.bond_length(atom.element, common_atom.element) # calculate the ideal length a bond between these atoms should have
				if (dist < ideal_bondlength * 1.2) and not common_index in already_used: # so an atom has been identified that should be bound. Also, this atom has not already been bound to some other atom
					toadd.append(common_atom) # add this atom to the list of atoms that will eventually be added to molecule_object
					already_used.append(common_index) # make sure this atom is never again selected as a handle
	
		for atom in toadd: molecule_object.all_atoms[molecule_object.max_index() + 1] = atom # go through each of the identified handle atoms, and add to the molecule_object
		# return list of handle atoms for this pdb object	
		return toadd
	
	def __frags_connected(self,frag1, frag2, molecule_object, handle_list): 
		"""Determines if two fragments of a pymolecule.Molecule object are connected
		
		Arguments:
		frag1 -- A list of atom indices corresponding to the atoms of the first fragment.
		frag2 -- A list of atom indices corresponding to the atoms of the second fragment.
		molecule_object -- The molecular model containing the two fragments (pymolecule.Molecule).
		handle_list -- A list of handle atoms corresponding to the molecule_object.
		
		Returns:
		A boolean, true if the fragments are connected, false otherwise.
		
		"""
		
		for index1 in frag1: # for through the indices in frag1
			atom1 = molecule_object.all_atoms[index1] # this is the atom corresponding to index1
			for index2 in frag2: # go through the indices in frag2
				atom2 = molecule_object.all_atoms[index2] # assign atom to that index
				dist = atom1.coordinates.distance_to_another_point(atom2.coordinates) # calculate the distance
				if (dist < molecule_object.bond_length(atom1.element, atom2.element) * 1.2) and not (atom1 in handle_list and atom2 in handle_list): # so the fragments have atoms that are bonded to each other, but make sure every handle atom is in a separate list
					return True # you could probably change this to return True and eliminate some lines above. Check.
		return False
	
	def __will_fragments_clash(self, fragments, handle_list1, handle_list2):
		"""Detect if the fragment combination that has been chosen will lead to steric clashes. Avoid creating merged molecules with steric clashes.
		
		Arguments:
		fragments -- A list of fragments (pymolecules.Molecule) to add to the common substructure. Each fragment should come from a different handle.
		handle_list1 -- List of handle atoms.
		handle_list2 -- List of handle atoms.
		
		Returns:
		A boolean and string, true and a text message if steric clashes are present, false and "" otherwise.
		
		"""

		for frag1_index in range(0,len(fragments) - 1):
			for frag2_index in range(frag1_index + 1,len(fragments)):
				#extract the fragments
				frag1 = fragments[frag1_index]
				frag2 = fragments[frag2_index]
				#start to iterate through all atom pairs within these two fragments
				for index1 in frag1.all_atoms: # for through the indices in frag1
					atom1 = frag1.all_atoms[index1] # this is the atom corresponding to index1
					for index2 in frag2.all_atoms: # go through the indices in frag2
						atom2 = frag2.all_atoms[index2] # assign atom to that index
						# check whether the two atoms come close
						if(atom1.coordinates.distance_to_another_point(atom2.coordinates)) < 1.5:
							flag_handles = False
							# so we know they are close, but they could still both be handle atoms, check that they are not!!
							#iterate through all the atoms in the handle list and check whether they are really close to one of the atoms in question						
							for handle_atom1 in handle_list1:
								if atom1.coordinates.distance_to_another_point(handle_atom1.coordinates) < 0.01 or atom2.coordinates.distance_to_another_point(handle_atom1.coordinates) < 0.01:
									flag_handles = True
							for handle_atom2 in handle_list2:
								if atom1.coordinates.distance_to_another_point(handle_atom2.coordinates) < 0.01 or atom2.coordinates.distance_to_another_point(handle_atom2.coordinates) < 0.01:
									flag_handles = True
							# so if the atoms are close and none of them are handle atoms, you will have a steric clash
							if not flag_handles:
								return True, "Merging this fragment combination would lead to steric clashes, this particular pdb will not be generated"
							
		# otherwise just return false
		return False, ""
	
	def __get_frags(self,the_pdb_object, handle_list):
		"""Separate a pymolecule.Molecule object with multiple, probably unconnected fragments into their own pymolecule.Molecule objects.
		
		Arguments:
		the_pdb_object -- A pymolecule.Molecule object.
		handle_list -- List of handle atoms.
		
		Returns:
		A list of python.Molecule objects, corresponding to each individual fragment.
		
		"""

		frag_indicies = [] # a list of lists, where each list contains the indices of the atoms in a single fragment
		
		# first put each atom into it's own list
		for index in the_pdb_object.all_atoms: frag_indicies.append([index])
		
		# so at this point, the list looks like this: [[0], [1], [2]...]
			
		# now go through all pairs of groups and combine ones that are close
		merge_made = True # this is true as long as it's possible to merge lists. Starts off true.
		while merge_made == True: # keep going as long as it's possible to merge lists.
			merge_made = False # assume you can no longer merge lists. The program has to prove this supposition wrong.
			for index1 in range(len(frag_indicies)-1): # go through each of the lists in the frag_indicies list, except the last one
				frag1 = frag_indicies[index1] # rename the list to make it easier to reference
				for index2 in range(index1 + 1, len(frag_indicies)): # go through all the lists whose index is greater than index1. To not double count.
					frag2 = frag_indicies[index2] # rename for easy reference
					if self.__frags_connected(frag1,frag2,the_pdb_object, handle_list) == True: # check to see if these two fragments are somehow connected
						frag_indicies[index1].extend(frag2) # if they are, add the indices of one to the other
						frag_indicies[index2] = [] # delete the indices of the other. Don't delete the list, but keep it empty, so as to not mess up the counting.
						merge_made = True # a merge was made, so the loop should continue
		
		while [] in frag_indicies: frag_indicies.remove([]) # delete all empty lists frmo the list of lists
		
		# Just convert the list of indices into pdb objects.
		frags = [] # a list of pdb objects
		for frag in frag_indicies:
			pdb_temp = Molecule2()
			for index in frag: pdb_temp.all_atoms[index] = the_pdb_object.all_atoms[index]
			frags.append(pdb_temp) # add the molecule_object to the list frags
			
		return frags # return the list of pdb objects
	
	def __closest_pdb_dist(self, pdb1, pdb2):
		"""Find the closest distance two pymolecule.Molecule objects come to each other.
		
		Arguments:
		pdb1 -- The first pymolecule.Molecule object.
		pdb2 -- The second pymolecule.Molecule object.
		
		Returns:
		A float, the distance between the two pymolecule.Molecule objects.
		
		"""
		
		closest = 999999999.9
		for index1 in pdb1.all_atoms:
			atom1 = pdb1.all_atoms[index1]
			for index2 in pdb2.all_atoms:
				atom2 = pdb2.all_atoms[index2]
				dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
				if dist < closest: closest = dist
		return closest
	
	def __find_combinations(self, combination, new_list, return_value):
		"""A recursive function to generate all possible combinations of elements within a list of lists.
		
		Arguments:
		combination -- A growing list of combinations
		new_list -- A list of lists from which the combinations are generated
		return_value -- A growing list of possible combinations.
		
		"""
		
		# define flag to keep track of when you backrack in the recursive call
		global flag_jump_back 
		flag_jump_back = False
		#remove first element from new_list and store it in temp
		temp = new_list.pop(0)
		#iterate through all the elements of temp
		for index in range(len(temp)):
			# when you just jumped back (in recursive call), you need to remove the last element from combination
			if( flag_jump_back == True):
				combination.pop()
				flag_jump_back = False
			# add the current element from temp to combination
			combination.append(temp[index])
			# if new_list still contains elements you need to go deeper into the recursion
			if len(new_list) > 0:	
				__find_combinations(copy.deepcopy(combination), copy.deepcopy(new_list), return_value)
			# if not, you have reached the end and can store the result
			else:
				return_value.append(copy.deepcopy(combination))
				#if the deepest elemennt still contains more entries, keep going, remove the last element from combination	
				if index < len(temp) - 1:
					combination.pop()
				# if not, you need to jump back up (in the recursive call history)
				else:
					flag_jump_back = True				
	
	
	def main(self, pdb1, pdb2, all_symmetry_relations, all_substituent_combinations, output_mcs, max_runtime):
		"""The main function that runs the LigMerge algorithm, combining two pymolecule.Molecule into one.
		
		Arguments:
		pdb1 -- The first pymolecule.Molecule.
		pdb2 -- The second pymolecule.Molecule.
		all_symmetry_relations -- A string, "TRUE" if the algorithm is to take into account molecular symmetry.
		all_substituent_combinations -- A string, "TRUE" if the algorithm should produce all possible substituent combinations.
		output_mcs -- A string, "TRUE" if the algorithm should output the maximum common substructure of the two pymolecule.Molecule objects.
		max_runtime -- An integer, the number of seconds 
		
		Returns:
		A list of pymolecule.Molecule objects, the possible combined products.
		
		"""
		
		global common_pdb
		
		message = ""
		# set start time of program call
		start_time = time.time()
	
		# define original position
		pdb1.set_undo_point()
		pdb2.set_undo_point()
	
		# calculate the largest common substructure from the two pdbs
		acommon_structures = CommonSubstructures()
		substructure_indices = acommon_structures.get_common_structures(pdb1,pdb2,start_time,max_runtime)
	
		if substructure_indices == []: return None, "No candidate for MCS with at least 3 atoms has been identified!!" 
		
		# at this point there may be several index combinations that are all symmetry related (like: ([1, 2, 3], [1, 2, 3]) vs ([1, 2, 3], [2, 3, 1])) but can constitute different combinations (if there are different non-H atom substituents)
	
		# if all_symmetry_relations flag is 
		if (all_symmetry_relations == "TRUE"):
			lower_bound = 0
			upper_bound = len(substructure_indices)
		else:
			random_index = random.randrange(0, len(substructure_indices), 1)
			lower_bound = random_index
			upper_bound = random_index + 1
		
		pdbs_list = []
		
		for substructure_combinations_index in range(lower_bound, upper_bound):
		
			pdb_to_move = Molecule2()
			pdb_orig = Molecule2()
	
			for index1 in pdb1.all_atoms:
				pdb_to_move.all_atoms[index1] = pdb1.all_atoms[index1].copy_of()
			for index2 in pdb2.all_atoms:
				pdb_orig.all_atoms[index2] = pdb2.all_atoms[index2].copy_of()
	
			# define original position
			pdb_to_move.set_undo_point()
			pdb_orig.set_undo_point()
	
			substructure_pdb1_indices = substructure_indices[substructure_combinations_index][0]
			substructure_pdb2_indices = substructure_indices[substructure_combinations_index][1]
	
			# create two pdbs that just contain the common substructure atoms in them
			common_substructure_pdb1 = Molecule2()
			common_substructure_pdb2 = Molecule2()
	
			for index1 in substructure_pdb1_indices:
				# the .copy_of() is very(!!!) important here; otherwise the new pdb object will contain the SAME atoms as the original pdb; whatever you do to one, you will also do to the other...
				common_substructure_pdb1.all_atoms[index1] = pdb_to_move.all_atoms[index1].copy_of()
			for index2 in substructure_pdb2_indices:
				common_substructure_pdb2.all_atoms[index2] = pdb_orig.all_atoms[index2].copy_of()
	
			# if flag for output of maximum common substructure is set, print out the pdbs
			if (output_mcs == "TRUE"):
				common_substructure_pdb1.save_pdb2(vars['output_dir'] + 'common_substructure_pdb1_' + repr(substructure_combinations_index) + '.pdb')
				common_substructure_pdb2.save_pdb2(vars['output_dir'] + 'common_substructure_pdb2_' + repr(substructure_combinations_index) + '.pdb')
	
			# define original position of substructures
			common_substructure_pdb1.set_undo_point()
			common_substructure_pdb2.set_undo_point()
	
			# the general idea here will be to find the best transformation that will superimpose the two common substructure fragments
			# this transformation will then be applied to the full structures
	
			# do several minimizations of the rmsd of the two common substructures and save best transformation
			best_rmsd = 999999999.9 # to keep track of what the best minimization is
			initial = numpy.ndarray(shape=(6), dtype=float)
			best_array = numpy.ndarray(shape=(6), dtype=float)
	
			# so minimize 25 times (this can be modified)
			for t in range(1):
		
				# move the two ligands randomly relative to each other
				initial[0] = 20.0 * random.random() - 10.0 
				initial[1] = 20.0 * random.random() - 10.0
				initial[2] = 20.0 * random.random() - 10.0
				initial[3] = (4.0 * random.random() - 2.0) * math.pi # no need to add tomove_array value, since already covering full range of rotations
				initial[4] = (4.0 * random.random() - 2.0) * math.pi
				initial[5] = (4.0 * random.random() - 2.0) * math.pi
	
				# each minimization is done in two steps
	
				# First do a minimization of the tethered atoms starting with a random set of transformations
				test_array = scipy.optimize.fmin_slsqp(self.__calculate_rmsd_assignments,initial,args=(common_substructure_pdb1,common_substructure_pdb2,substructure_indices[substructure_combinations_index]),iprint=-1) # fmin_slsqp was faster than powell, fmin, fmin_cg, fmin_bfgs, and anneal in a test case
	
				# Then do a second minimization of the tethered atoms starting with the transformation determined in the first step
				# note: I don't think this is really necessary!!!
				test_array = scipy.optimize.fmin_slsqp(self.__calculate_rmsd_assignments,test_array,args=(common_substructure_pdb1,common_substructure_pdb2,substructure_indices[substructure_combinations_index]),iprint=-1) # fmin_slsqp was faster than powell, fmin, fmin_cg, fmin_bfgs, and anneal in a test case
		
				# calculate the rmsd of the minimization
				rmsd = self.__calculate_rmsd_assignments(test_array,common_substructure_pdb1,common_substructure_pdb2,substructure_indices[substructure_combinations_index])
	
				# now check to see if this minimization has a better rmsd than the one on file
				if rmsd < best_rmsd:
					best_rmsd = rmsd
					best_array = test_array
	
			# now we have identified the transformation that best superimposes the two substructures
			# we will now determine the transformation that is necessary to superimpose the entire structures
	
			# for this we need the geometric center of the original common_substructure_pdb1
			# so set it back to its original position
			common_substructure_pdb1.undo()
	
			# calculate its center
			center = common_substructure_pdb1.get_center()
	
			# and add back the translational information from the best transformation
			center[0] += best_array[0]
			center[1] += best_array[1]
			center[2] += best_array[2]
	
			# now finally transform the pdb_to_move with the transformation best_array around the calculated center
			self.__transform_molecule(best_array,pdb_to_move, center)
	
			# just for test purposes also transform the common_substructure_pdb1 back to its best superimposed position
			#self.__transform_molecule(best_array,common_substructure_pdb1)
	
			# so now you have aligned structures.
	
			# now identify all the atoms that are superimposed
			pairs = [] # this is a list containg couples corresponding to the indices of atoms that are superimposed
			for index1 in pdb_to_move.all_atoms: # go through the indices in the first pdb
				atom1 = pdb_to_move.all_atoms[index1] # define the atom from the index
				for index2 in pdb_orig.all_atoms: # go through the indices in the second pdb
					atom2 = pdb_orig.all_atoms[index2] # define the atom from the index
					dist = atom1.coordinates.distance_to_another_point(atom2.coordinates) # get the distance between these atoms
					if dist < 0.2 and atom1.element.strip() == atom2.element.strip() : # they are very close together (had to increase that value to make sure we identify all superimposed atoms correctly)
						#superimposed[str(index1) + "_to_move"] = index2
						#superimposed[str(index2) + "_orig"] = index1
						pairs.append((index1, index2)) # add the couple to the list of pairs of atoms that are superimposed
	
			# make a single pdb that only contains the superimposed atoms; delete those atoms from the other pdbs
			common_pdb = Molecule2() # common_pdb is the name of the PDB object containing the atoms they have in common (superimposed)
			for pair in pairs: # go through each of the superimposed atoms in the list.
		
				# renaming for convenience
				index1 = pair[0]
				index2 = pair[1]
		
				# Add to the list of atoms in common_pdb the common atom from pdb_to_move
				common_pdb.all_atoms[len(common_pdb.all_atoms)+1] = pdb_to_move.all_atoms[index1]
		
				# delete this common atom from pdb_to_move and pdb_orig
				pdb_to_move.delete_atom(index1)
				pdb_orig.delete_atom(index2)
	
			# at this stage check that the common pdb contain at least 3 atoms
			# if it doesn't this means that the superimposition did not overlay at least 3 atoms directly, i.e. the structures were too dissimilar to be superimposed (didn't have a common core)
			if (len(common_pdb.all_atoms) < 3):
				return None, "WARNING: Molecules are not superimposable. There are no geometrically similar common substructures containing more than three atoms."
	
			# now go through each of the atoms in pdb_to_move and pdb_orig and add back any atom from common_pdb that is immediately bound to it (essentially "handles")
			# add hanndles to each of the non-common ligands.
	
			handles_pdb_to_move = self.__add_handles(pdb_to_move)
			handles_pdb_orig = self.__add_handles(pdb_orig)
	
			# at this point, pdb_to_move and pdb_orig contain multiple moieties with they associated handle atoms. These fragments are not separate, but all together.
	
			# separate each of the fragments into their own pdb objects
			frags_to_move = self.__get_frags(pdb_to_move, handles_pdb_to_move)
			frags_orig = self.__get_frags(pdb_orig, handles_pdb_orig)
	
			#######################################################################
			## check for fragments that connect multiple handle atoms externally ##
			#######################################################################
	
			## identify fragments that have more than one handle atoms (due to all the splitting up that we did in __get_frags and __frags_connected this can only happen if the handle atoms are connected through the fragment itself (not through the common core))
	
			multiple_handle_frags_to_move = [] # a list of pdb objects that contain fragments that have more than one handle atom
			multiple_handle_frags_orig = [] # a list of pdb objects that contain fragments that have more than one handle atom
	
			# iterate through fragments associated with the pdb_to_move object
			for to_move_index in range(len(frags_to_move)): # go through the fragments that were originally associated with the pdb_to_move object
				frag_to_move = frags_to_move[to_move_index] #frag_to_move is the actual pdb object of the fragment
				handle_counter = 0
				for atom_index in frag_to_move.all_atoms: # go through each atom in that molecule_object
					atom = frag_to_move.all_atoms[atom_index] # get the actual atom
					if atom in handles_pdb_to_move: # check whether the atom is a handle
						handle_counter += 1 #increment handle counter
					if handle_counter == 2: #as soon as you hit 2 handles per fragment
						multiple_handle_frags_to_move.append(frag_to_move) # append the fragment to the list
						break # break out of the loop (no need to iterate over more atom, we already found 2 handle atoms)
	
			# iterate through fragments associated with the pdb_orig object
			for orig_index in range(len(frags_orig)): # go through the fragments that were originally associated with the pdb_orig object
				frag_orig = frags_orig[orig_index] #frag_orig is the actual pdb object of the fragment
				handle_counter = 0
				for atom_index in frag_orig.all_atoms: # go through each atom in that molecule_object
					atom = frag_orig.all_atoms[atom_index] # get the actual atom
					if atom in handles_pdb_orig: # check whether the atom is a handle
						handle_counter += 1 #increment handle counter
					if handle_counter == 2: #as soon as you hit 2 handles per fragment
						multiple_handle_frags_orig.append(frag_orig) # append the fragment to the list
						break # break out of the loop (no need to iterate over more atom, we already found 2 handle atoms)
	
			problematic_handle_atoms = [] # list of handle atoms that are in multiple-handle fragment in both compounds
	
			# the bug with the missing hydrogens should only occur if there are overlapping multiple-handle fragment from both compounds
			if len(multiple_handle_frags_to_move) >= 1 and len(multiple_handle_frags_orig) >= 1: #you do have multiple-handle fragment from both compounds
				for to_move_index in range(len(multiple_handle_frags_to_move)): # go through the multiple-handle fragments that were originally associated with the pdb_to_move object
					multiple_handle_frag_to_move = multiple_handle_frags_to_move[to_move_index] #multiple_handle_frag_to_move is the actual pdb object of the fragment
					for to_move_atom_index in multiple_handle_frag_to_move.all_atoms: # go through each atom in that molecule_object
						to_move_atom = multiple_handle_frag_to_move.all_atoms[to_move_atom_index] # get the actual atom
						for orig_index in range(len(multiple_handle_frags_orig)): # go through the multiple-handle fragments that were originally associated with the pdb_orig object
							multiple_handle_frag_orig = multiple_handle_frags_orig[orig_index] #multiple_handle_frag_to_move is the actual pdb object of the fragment
							for orig_atom_index in multiple_handle_frag_orig.all_atoms: # go through each atom in that molecule_object
								orig_atom = multiple_handle_frag_orig.all_atoms[orig_atom_index] # get the actual atom
								dist = to_move_atom.coordinates.distance_to_another_point(orig_atom.coordinates) # calculate the distance between the atoms
								if dist <= 0.01: #if they are the same atom
									problematic_handle_atoms.append(to_move_atom)
	
			# the idea now is to randomly pick one of the two compounds and have the merged structure get its problem fragments (the externally connected ones)
			# delete problem fragments from other compound
			# in addition also delete all the fragments that overlap with the chosen problem fragment
			if len(problematic_handle_atoms) >= 1: #if there is the case that you have multiple handle atoms that are connected externally in both compounds
				if random.choice([1,2]) == 1: # randomly choose which compounds problem fragments will be deleted
					for frag_to_be_removed in multiple_handle_frags_to_move: #delete problem fragments associted with pdb_to_move object
						frags_to_move.remove(frag_to_be_removed)
					for frag_no_overlap_possible in multiple_handle_frags_orig: #delete all fragments associated with pdb_to_move object that overlap with the retained problem fragment associted with pdb_orig object
						for potentially_overlapping_frag in frags_to_move:
							if self.__closest_pdb_dist(frag_no_overlap_possible, potentially_overlapping_frag) < 0.25:
								frags_to_move.remove(potentially_overlapping_frag)
				else:
					for frag_to_be_removed in multiple_handle_frags_orig: #delete problem fragments associted with pdb_orig object
						frags_orig.remove(frag_to_be_removed)
					for frag_no_overlap_possible in multiple_handle_frags_to_move: #delete all fragments associated with pdb_orig object that overlap with the retained problem fragment associted with pdb_to_move object
						for potentially_overlapping_frag in frags_orig:
							if self.__closest_pdb_dist(frag_no_overlap_possible, potentially_overlapping_frag) < 0.25:
								frags_orig.remove(potentially_overlapping_frag)
	
			###########################################################################
			## end check for fragments that connect multiple handle atoms externally ##
			###########################################################################
	
			# create composite compounds with all possible combinations of fragments
			if (all_substituent_combinations == "TRUE"):
				# now figure out which fragments go with each other (overlap); this is the same as seeing which ones have the same handle 
				# whenever you have a multiple-handle fragment and it is picked in the random process, make sure that it is picked at its second handle automatically (otherwise you end up with lost fragments)
				fragments_that_have_to_be_picked = []
				fragments_already_picked = []
				fragments_with_same_handles = []
	
				# first of all take care of all the multiple handle atoms
				# find all the single handle fragments from orig that overlap with multiple handle fragment from to_move
				for multiple_handle_to_move_index in range(len(multiple_handle_frags_to_move)):
					multiple_handle_to_move = multiple_handle_frags_to_move[multiple_handle_to_move_index]
					frag_go_together = [[],[]]
					frag_go_together[0] = multiple_handle_to_move
					for orig_index in range(len(frags_orig)): # go through the fragments that were originally associated with the pdb_orig object
						frag_orig = frags_orig[orig_index]
						if self.__closest_pdb_dist(multiple_handle_to_move, frag_orig) < 0.25: # so fragments come within 0.25 angstroms of each other
							frag_go_together[1].append(frag_orig)
					fragments_with_same_handles.append(frag_go_together)
	
				# find all the single handle fragments from to_move that overlap with multiple handle fragment from orig
				for multiple_handle_orig_index in range(len(multiple_handle_frags_orig)):
					multiple_handle_orig = multiple_handle_frags_orig[multiple_handle_orig_index]
					frag_go_together = [[],[]]
					frag_go_together[0] = multiple_handle_orig
					for to_move_index in range(len(frags_to_move)): # go through the fragments that were originally associated with the pdb_orig object
						frag_to_move = frags_to_move[to_move_index]
						if self.__closest_pdb_dist(multiple_handle_orig, frag_to_move) < 0.25: # so fragments come within 0.25 angstroms of each other
							frag_go_together[1].append(frag_to_move)
					fragments_with_same_handles.append(frag_go_together)
	
	
				#now go through the to_move fragments and find the corresponding overlapping fragment in orig
				#skip all of this if to_move fragment is multiple handle or a corresponding orig fragments is multiple handle; they are already in the list, see above
				for to_move_index in range(len(frags_to_move)): # go through the fragments that were originally associated with the pdb_to_move object
					frag_to_move = frags_to_move[to_move_index]
					if frag_to_move in multiple_handle_frags_to_move: continue
					frag_go_together = []
					frag_go_together.insert(0,frag_to_move)
					for orig_index in range(len(frags_orig)): # go through the fragments that were originally associated with the pdb_orig object
						frag_orig = frags_orig[orig_index]
						if self.__closest_pdb_dist(frag_to_move, frag_orig) < 0.25: # so fragments come within 0.25 angstroms of each other, delete one.
							if frag_orig in multiple_handle_frags_orig: continue
							frag_go_together.insert(0,frag_orig)
					if len(frag_go_together) > 1:
						fragments_with_same_handles.append(frag_go_together)
	
				#define local variable that will hold the results of the search for all possible combinations
				return_value_find_combinations = []
			
	
				# at this stage you should have a list where each entry are fragments with the same handle atoms (fragments_with_same_handles)
				start_combinations =[]
				if not len(fragments_with_same_handles) == 0:
					__find_combinations(start_combinations, fragments_with_same_handles, return_value_find_combinations)
	
				# in case you had multiple handle fragments your list of possible combinations will not only contain pure fragments but sometimes lists of fragments (this will be the case when a multiple handle fragment was discarded and multiple single handle fragments have to be chosen simultaniously)
				# we need to consolidate all the entries into single fragment entries, for this split up the lists of fragments
	
				for fragment_combination in return_value_find_combinations:
					for element_index in range(len(fragment_combination) - 1, -1, -1):
						element_in_fragment_combination = fragment_combination[element_index]
						if type(element_in_fragment_combination) == list and len(element_in_fragment_combination) >= 1:
							for individual_fragment_in_list in fragment_combination[element_index]:
								fragment_combination.append(individual_fragment_in_list)
							fragment_combination.remove(fragment_combination[element_index])
						elif type(element_in_fragment_combination) == list and len(element_in_fragment_combination) == 0:
							fragment_combination.remove(fragment_combination[element_index])
	
	
				substituent_combinations_counter = 0;
	
				#iterate through all the fragment combinations
				for fragment_combination in return_value_find_combinations:
	
					# check that the merged pdb cannot contain atoms that are too close in space
					result,message = self.__will_fragments_clash(fragment_combination, handles_pdb_to_move, handles_pdb_orig)
					# if this fragment combination leads to a steric clash, don't build it and give out warning, don't assert as we still need to test the other combinations
					if result:
						print message
						continue
	
					all_frags = fragment_combination
					substituent_combinations_counter = substituent_combinations_counter + 1
	
					# merge the pdbs into one
					i = 0
					merged_pdb = Molecule2() # define a new PDB object called merged_pdb
	
					for index in common_pdb.all_atoms: # add in the atoms from the common pdb into the merged_pdb
						i = i + 1 # keep track of pdb-output atom indices
						atom = common_pdb.all_atoms[index].copy_of()
						atom.PDBIndex = str(i)
						merged_pdb.all_atoms[i] = atom
	
					for frag in all_frags: # now go through all the fragments and add them into the merged pdb
						for frag_atom_index in frag.all_atoms: # look at each atom index in each of the fragments
							frag_atom = frag.all_atoms[frag_atom_index] # here's where the atom is actually defined from the index above
							too_close = False # assume the atoms are too close. Disprove this supposition.
							for common_index in common_pdb.all_atoms: # go through each of the indices in common
								common_atom = common_pdb.all_atoms[common_index] # here's the atom of that index
								dist = frag_atom.coordinates.distance_to_another_point(common_atom.coordinates) # calculate the distance between fragment atom and the common atom
								if dist < 0.25: # they are very close, this is probably a handle.
									too_close = True # they are too close, the fragment atom is probably a handle.
									break # don't proceed with adding this atom to the merged pdb
							if too_close == False: # they are not too close, add to merge pdb
								i = i + 1 # keep track of pdb-output atom indices
								atom = frag_atom.copy_of()
								atom.PDBIndex = str(i)
								merged_pdb.all_atoms[i] = atom
				
					# save the merged pdb!!!!
					pdbs_list.append(merged_pdb)
	
	
			else:
				# now figure out which fragments go with each other (overlap) and randomly select one or the other where appropriate to create a composite compound.
				# this is the same as seeing which ones have the same handle.
	
				# frags_to_move = a list of pdb objects
	
				# whenever you have a multiple-handle fragment and it is picked in the random process, make sure that it is picked at its second handle automatically (otherwise you end up with lost fragments)
				fragments_that_have_to_be_picked = []
				fragments_already_picked = []
	
				for to_move_index in range(len(frags_to_move)): # go through the fragments that were originally associated with the pdb_to_move object
					frag_to_move = frags_to_move[to_move_index]
					for orig_index in range(len(frags_orig)): # go through the fragments that were originally associated with the pdb_orig object
						frag_orig = frags_orig[orig_index]
						if self.__closest_pdb_dist(frag_to_move, frag_orig) < 0.25: # so fragments come within 0.25 angstroms of each other, delete one.
							if frag_to_move in fragments_that_have_to_be_picked: #if frag_to_move has been previously identifies as multiple handle fragment that will always be picked from here on
								frags_orig[orig_index].all_atoms = {} #delete frag_orig
								continue #skip the random picking part for this fragment pair
							if frag_orig in fragments_that_have_to_be_picked: #if frag_orig has been previously identifies as multiple handle fragment that will always be picked from here on
								frags_to_move[to_move_index].all_atoms = {} #delete frag_to_move
								continue #skip the random picking part for this fragment pair
							if random.choice([1,2]) == 1: # randomly choose which fragment will be deleted
								frags_to_move[to_move_index].all_atoms = {} #delete frag_to_move
								if frag_orig in multiple_handle_frags_orig: #check whether frag_to_move is a multiple handle fragment
									fragments_that_have_to_be_picked.append(frag_orig) #if it is, make sure it will always be picked from here on (if it wouldn't this handle atom had no fragment)
							else:
								frags_orig[orig_index].all_atoms = {} #delete frag_orig
								if frag_to_move in multiple_handle_frags_to_move: #check whether frag_orig is a multiple handle fragment
									fragments_that_have_to_be_picked.append(frag_to_move) #if it is, make sure it will always be picked from here on (if it wouldn't this handle atom had no fragment)
	
				# move all the fragments, regardless of the pdb object from which they were derived, into a single list
				all_frags = frags_to_move[:]
				all_frags.extend(frags_orig)
	
				# check that the merged pdb cannot contain atoms that are too close in space
				result,message = self.__will_fragments_clash(all_frags, handles_pdb_to_move, handles_pdb_orig)
				# if this fragment combination leads to a steric clash, don't build it and give out warning, don't assert as we still need to test the other combinations
				if result:
					print message
					continue
	
	
				# merge the pdbs into one
				i = 0
				merged_pdb = Molecule2() # define a new PDB object called merged_pdb
	
				for index in common_pdb.all_atoms: # add in the atoms from the common pdb into the merged_pdb
					i = i + 1 # keep track of pdb-output atom indices
					atom = common_pdb.all_atoms[index].copy_of()
					atom.PDBIndex = str(i)
					merged_pdb.all_atoms[i] = atom
	
				for frag in all_frags: # now go through all the fragments and add them into the merged pdb
					for frag_atom_index in frag.all_atoms: # look at each atom index in each of the fragments
						frag_atom = frag.all_atoms[frag_atom_index] # here's where the atom is actually defined from the index above
						too_close = False # assume the atoms are too close. Disprove this supposition.
						for common_index in common_pdb.all_atoms: # go through each of the indices in common
							common_atom = common_pdb.all_atoms[common_index] # here's the atom of that index
							dist = frag_atom.coordinates.distance_to_another_point(common_atom.coordinates) # calculate the distance between fragment atom and the common atom
							if dist < 0.25: # they are very close, this is probably a handle.
								too_close = True # they are too close, the fragment atom is probably a handle.
								break # doin't proceed with adding this atom to the merged pdb
						if too_close == False: # they are not too close, add to merge pdb
							i = i + 1 # keep track of pdb-output atom indices
							atom = frag_atom.copy_of()
							atom.PDBIndex = str(i)
							merged_pdb.all_atoms[i] = atom
				
				# save the merged pdb!!!!
				pdbs_list.append(merged_pdb)

		return pdbs_list, ""

class CommonSubstructures:
	"""A class for analyzing two pymolecule.Molecule objects to identify common substructures."""

	def __do_assignment_atoms_overlay(self, assignment, pdb1, pdb2): #input is a list of 2 lists, and the two original pdbs to have the coordinates
		"""Checks whether or not an assignment actually corresponds to overlying atoms.
		
		Arguments:
		assignment -- A list.
		pdb1 -- A pymolecule.Molecule object.
		pdb2 -- A pymolecule.Molecule object.
		
		Returns:
		A boolean, True if there is an overlay, False otherwise.
		
		"""
		
		#create pdbs for the two structures in assignment
		pdb1_substructure = Molecule2()
		pdb2_substructure = Molecule2()
		for index in range(len(assignment[0])): 
			pdb1_substructure.all_atoms[assignment[0][index]] = pdb1.all_atoms[assignment[0][index]]
			pdb2_substructure.all_atoms[assignment[1][index]] = pdb2.all_atoms[assignment[1][index]]
		
		#iterate through all possible atom pairs in pdb1 and the corresponding atom pairs in pdb2 (correpondance is made by assignment)
		for index_first_atom in range(len(assignment[0])):
			for index_second_atom in range(index_first_atom + 1, len(assignment[0])):
				index1_first = assignment[0][index_first_atom]
				index2_first = assignment[1][index_first_atom]
				index1_second = assignment[0][index_second_atom]
				index2_second = assignment[1][index_second_atom]
				atom1_first = pdb1.all_atoms[index1_first]
				atom2_first = pdb2.all_atoms[index2_first]
				atom1_second = pdb1.all_atoms[index1_second]
				atom2_second = pdb2.all_atoms[index2_second]
	
				# for each of the atom pairs in the assignment test the following: 
				# 1) is the distance between the two atoms in pdb1 the same as in pdb2 (allow for small deviation)
				# 2) do the pairs of atoms correspond to the same elements in both pdbs
				if math.fabs(atom1_first.coordinates.distance_to_another_point(atom1_second.coordinates) - atom2_first.coordinates.distance_to_another_point(atom2_second.coordinates)) > 0.02 or atom1_first.element != atom2_first.element or atom1_second.element != atom2_second.element :
					#if only one of these instances is violated, then the assignment doesn't correspond to overlaying atoms
					return False
		#if none of the violations are triggered, the assignment does correspond to overlaying atoms
		return True
	
	def get_common_structures(self, pdb1, pdb2, start_time, time_limit):
		"""Get the common substructures of two pymolecule.Molecule objects.
		
		Arguments:
		pdb1 -- A pymolecule.Molecule object.
		pdb2 -- A pymolecule.Molecule object.
		start_time -- The number of seconds since the unix epoch when the LigMerge algorithm was started (Integer).
		time_limit -- The maximum number of seconds to allow the LigMerge algorithm to run (Integer).
		
		Returns:
		A list of lists, specifying which atoms in each of the pymolecule.Molecule objects are in common substructures.
		
		"""

		# first, generate a list of atoms (indices) and elements
		pdb1_lists = []
		for index in pdb1.all_atoms:
			atom = pdb1.all_atoms[index]
			if atom.element != "H": pdb1_lists.append([[index],[atom.element]]) # ignoring hydrogens
	
		pdb2_lists = []
		for index in pdb2.all_atoms:
			atom = pdb2.all_atoms[index]
			if atom.element != "H": pdb2_lists.append([[index],[atom.element]]) # ignoring hydrogens
		
		#at this point you have a list just like this one for each of the two molecules
		#[[[1], ['C']], [[2], ['C']], [[3], ['C']], [[4], ['C']], [[5], ['C']], [[6], ['C']], [[7], ['C']], [[8], ['C']], [[9], ['C']], [[10], ['C']], [[11], ['C']], [[12], ['C']], [[13], ['N']]]
		# this of this as a list of pairs, where each pair contains index/element info
	
		no_repeats = [] # to make sure no repeated substructures are included. 1, 2, 3 is the same as 3, 2, 1
		similar_substructures = [] # to keep track of substructures that are similar
		
		# So, the growing list similar_substructures keeps track of all potentially similar substructures, not just the biggest ones, in case the biggest ones among the two ligands end up not being geometrically similar
		
		for t in range(50): # so the maximum size of common fragments is 50 heavy atoms. probably way to big.
			# now go through the lists and remove elements that do only occur in one of the two pdbs
			# this is just to speed things up, since there's no chance these atoms could participate in common structures
			for index1 in range(len(pdb1_lists)-1,-1,-1): # going through the list backwards so I can delete items as I go...
				list1 = pdb1_lists[index1]
				
				# now look through the other list and see if it exists there.
				found = False
				for index2 in range(len(pdb2_lists)-1,-1,-1):
					list2 = pdb2_lists[index2]
					if list1[1] == list2[1]: #if there is an atom of the same element type in the other pdb
						found = True
						break
					
				if found == False: # so it's not found
					pdb1_lists.pop(index1) # remove it
					
			# now go the other way as well.
			for index2 in range(len(pdb2_lists)-1,-1,-1):
				list2 = pdb2_lists[index2]
				
				# now look through the other list and see if it exists there.
				found = False
				for index1 in range(len(pdb1_lists)-1,-1,-1):
					list1 = pdb1_lists[index1]
					if list2[1] == list1[1]: #if there is an atom of the same element type in the other pdb
						found = True
						break
					
				if found == False: # so it's not found
					pdb2_lists.pop(index2)

			#this starts off with two lists (one for each pdb) that contains all the index-element pairs of elements that are present in both structures
			#in every step of the iteration (for t in range(50)) the list members grow by one element; i.e. you will have sets of t number of atoms that are connected to each other, that occur in both pdbs
			# now match stuff up. What substructres on one pdb are similar to the other
			
			for item1 in pdb1_lists: # item1 is serial/element pair #pdb1_lists looks like this: [[[1], ['C']], [[2], ['C']], [[3], ['C']], [[4], ['C']], [[5], ['C']], [[6], ['C']], [[7], ['C']], [[8], ['C']], [[9], ['C']], [[10], ['C']], [[11], ['C']], [[12], ['C']], [[13], ['N']]]
				if len(item1[0]) > 3: # so you don't care about substructures that have only 3 atoms or less
					for item2 in pdb2_lists:
						if len(item2[0]) > 3: # so you don't care about substructures that have only 3 atoms or less
							if item1[1] == item2[1]: # so element list is the same (you are comparing lists, e.g. [C, S, C] against [C, S, C])
								
								# make a key that uniquely identifies these substructures regardless of index order
								tmp1 = item1[0][:] #make a real copy of the element
								tmp1.sort() # this sorting is just to make a unique key. The item's themselves are not sorted.
								tmp2 = item2[0][:] #make a real copy of the element
								tmp2.sort() # this sorting is just to make a unique key. The item's themselves are not sorted.
								key = (tmp1, tmp2) # so a duplet of the indices associated with each
	
								if key not in no_repeats: # so this structure has not been detected before
									no_repeats.append(key)
									similar_substructures.append((item1[0], item2[0]))
									# now you have this long list (similar_substructures) that has ordered pairs of sets of indicies in it that correspond to matching atoms
									
			# okay, now add the neighboring atoms to the growing list
			# so this is just generating lists of atoms that are bound together, first single atoms, then pairs, then tripplets, etc.
			# first, the first list
			temp_pdb1_lists = pdb1_lists[:]
			pdb1_lists = []
			for item in temp_pdb1_lists:
				
				# get a list of all atoms connected to the last atom in the pdb1_list
				indices = item[0]
				last_index = indices[-1]
				elements = item[1]
				connected_atoms = pdb1.all_atoms[last_index].indecies_of_atoms_connecting
				
				for connected_atom in connected_atoms:
					if connected_atom not in indices: # so no backtracking
						atom = pdb1.all_atoms[connected_atom]
						if atom.element != "H": # so no hydrogens
							new_indices = indices[:]
							new_indices.append(connected_atom)
							new_elements = elements[:]
							new_elements.append(atom.element)
							pdb1_lists.append([new_indices, new_elements])
				
			# now do the same for the other list
			temp_pdb2_lists = pdb2_lists[:]
			pdb2_lists = []
			for item in temp_pdb2_lists:
				indices = item[0]
				last_index = indices[-1]
				elements = item[1]
				connected_atoms = pdb2.all_atoms[last_index].indecies_of_atoms_connecting
				
				for connected_atom in connected_atoms:
					if connected_atom not in indices: # so no backtracking
						atom = pdb2.all_atoms[connected_atom]
						if atom.element != "H": # so no hydrogens
							new_indices = indices[:]
							new_indices.append(connected_atom)
							new_elements = elements[:]
							new_elements.append(atom.element)
							pdb2_lists.append([new_indices, new_elements])
			# check whether the excecution time until this point comes within 10s of time_limit
			if( int(time.time()) - int(start_time) > int(time_limit)):
				log("WARNING: the execution time of the subroutine identifying the largest common substructure has exceeded the maximum run time provided in the command line; a set of smaller substructures is used for the superposition.")
				break
	
		# it is very important to note that at this stage similar_substructures does not contain assignments, it just contains stretches of indeces of neighboring atoms that occur in both pdbs (the only thing that is compared there is the element type)
		# example: ([1, 2, 3, 4], [10, 9, 8, 7]) means that in pdb1 atoms 1,2,3,4 are bound to each other (10,9,8,7 are connected in pdb2) and they represent the same sequence of elements (like for instance C,C,S,O); at this stage this does NOT associate atom 1 (pdb1) with atom 10 (pdb2) and so on
	
		# the next step is the most crucial one, i.e. converting the lists into assignments (where atoms are actually associated with each other)
	
		# create list that will hold the actual valid assignments
		most_similar_substructures = []
	
		#determine the size of the longest list in similar_substructures
		size_longest_common_substructure = len(similar_substructures[-1][0])
	
		#iterate from longest list to subsequently shorter lists 
		for size_common_substructure in range(size_longest_common_substructure, 3, -1):
	
			# iterate through all the assignments in the list
			for assignment_pair in similar_substructures:
				#assignment_pair = copy.deepcopy(similar_substructures[k])
				
				# look at all the assignment pairs with this length
				if len(assignment_pair[0]) == size_common_substructure:
					# for each of these assignment lists permutate the second list entry to get all possible combinations (this corresponds to checking all the permutations related by rotations)
					for i in range(0,len(assignment_pair[1])):
						#this is the actual permutation
						assignment_pair[1].append(assignment_pair[1].pop(0))
						#test whether this assignment corresponds to overlaying compounds
						if self.__do_assignment_atoms_overlay( assignment_pair, pdb1, pdb2):
							#if it does append the assignment pair to the list
							most_similar_substructures.append(copy.deepcopy(assignment_pair))
	
					# reverse the order of the second list entry (corresponds to flipping it) and then permutate the second list entry to get all possible combinations (corresponding to checking all the permutations of the flipped one by rotations)
					# inversion!!
					assignment_pair[1].reverse()
					for i in range(0,len(assignment_pair[1])):
						#this is the actual permutation
						assignment_pair[1].append(assignment_pair[1].pop(0))
						#test whether this assignment corresponds to overlaying compounds
						if self.__do_assignment_atoms_overlay( assignment_pair, pdb1, pdb2):
							#if it does append the assignment pair to the list
							most_similar_substructures.append(copy.deepcopy(assignment_pair))
	
		# now, discard substructures that don't have the maximum common atoms.
		max_possible = 0
		for similar in most_similar_substructures:
			if len(similar[0]) > max_possible:
				max_possible = len(similar[0])
		for index in range(len(most_similar_substructures)-1,-1,-1):
			similar = most_similar_substructures[index]
			if len(similar[0]) < max_possible:
				most_similar_substructures.pop(index)
				x = 0
	
		return most_similar_substructures

###################################################
# SECTION OF CODE TO INTERFACE WITH AUTOCLICKCHEM #
###################################################

def get_pdb_files(loc):
	"""Get a list of all the PDB files in a given directory
	
	Arguments:
	loc -- A string, the directory to search.
	
	Returns:
	A list, the names of the identified PDB files.
	
	"""
	
	files = []
	if os.path.isdir(loc): # so it's a directory, go through the directory and find all the pdb files
		if loc[-1:]!=os.sep: loc = loc + os.sep # so add a / to the end of the directory
		files.extend(glob.glob(loc + '*.pdb'))
		files.extend(glob.glob(loc + '*.PDB'))
	else: # so it's a file
		if loc[-3:]=="PDB" or loc[-3:]=="pdb":
			files.append(loc)
		else:
			log_file("The file " + loc + " does not have the PDB file extention and so cannot be used.",log)
	return files

def run_autoclick_chem(args):
	"""Run the autoclickchem module.
	
	Arguments:
	args -- a list containing the autoclickchem command-line parameters to pass.

	Returns:
	A list, the names of the PDB files that were created.
	
	"""

	kinds_of_reactions = ["azide_and_alkyne_to_azole", "epoxide_alcohol_opening", "epoxide_thiol_opening", "chloroformate_and_amine_to_carbamate", "sulfonyl_azide_and_thio_acid", "carboxylate_and_alcohol_to_ester", "carboxylate_and_thiol_to_thioester", "acyl_halide_and_alcohol_to_ester", "acyl_halide_and_thiol_to_thioester", "ester_and_alcohol_to_ester", "ester_and_thiol_to_thioester", "acid_anhydride_and_alcohol_to_ester", "acid_anhydride_and_thiol_to_thioester", "carboxylate_and_amine_to_amide", "acyl_halide_and_amine_to_amide", "ester_and_amine_to_amide", "acid_anhydride_and_amine_to_amide", "isocyanate_and_amine_to_urea", "isothiocyanate_and_amine_to_thiourea", "isocyanate_and_alcohol_to_carbamate", "isothiocyanate_and_alcohol_to_carbamothioate", "isocyanate_and_thiol_to_carbamothioate", "isothiocyanate_and_thiol_to_carbamodithioate", "alkene_to_epoxide", "halide_to_cyanide", "alcohol_to_cyanide", "carboxylate_to_cyanide", "acyl_halide_to_cyanide", "acid_anhydride_to_cyanide", "halide_to_azide", "alcohol_to_azide", "carboxylate_to_azide", "acyl_halide_to_azide", "acid_anhydride_to_azide", "amine_to_azide", "amine_to_isocyanate", "amine_to_isothiocyanate", "azide_to_amine"]
	
	# build autoclickchem intermediates directory by intilizing an AutoClickChem object
	useless_object = autoclickchem.AutoClickChem()
	
	# start by enabling all reactions
	reactions_to_perform = kinds_of_reactions[:]
	
	# Get reactants
	reactants1 = []
	reactants2 = []
	
	output_dir = "." + os.sep
	
	for index in range(1,len(args)):
		arg = args[index]
		if arg.lower()=="-reactants1":
			loc = args[index + 1]
			reactants1.extend(get_pdb_files(loc)) # so it's a directory
			args[index] = ""
			args[index + 1] = ""
		elif arg.lower()=="-reactants2":
			loc = args[index + 1]
			reactants2.extend(get_pdb_files(loc)) # so it's a directory
			args[index] = ""
			args[index + 1] = ""
		elif arg.lower()=="-output_dir":
			
			output_dir = args[index + 1]
			if output_dir[-1:] != os.sep: output_dir = output_dir + os.sep
			
			# if the output directory doesn't exist, create it
			if not os.path.exists(output_dir):
				os.mkdir(output_dir)
			
			args[index] = ""
			args[index + 1] = ""
		elif arg.lower()=="-max_reactions":
			max_reactions = int(args[index + 1])
			args[index] = ""
			args[index + 1] = ""
		elif arg.lower()=="-all_reactions":
			reactions_to_perform = []
			args[index] = ""
		elif arg.lower()=="+all_reactions":
			reactions_to_perform = kinds_of_reactions[:]
			args[index] = ""
		elif arg[:1]=="+": # so it starts with a plus
			tag = arg[1:].lower()
			if tag in kinds_of_reactions: # so it's a valid tag
				if tag not in reactions_to_perform: # so it's not already in the list of reactions that will be performed
					reactions_to_perform.append(tag) # add it to the list
					args[index] = ""
		elif arg[:1]=="-": # so it starts with a minus
			tag = arg[1:].lower()
			if tag in reactions_to_perform:
				while tag in reactions_to_perform: # just in case for some reason it appears more than once in the list
					reactions_to_perform.remove(tag)
					args[index] = ""
	
	# build a list of the possible reactants
	reactions = []
	# first, consider all the reactions where there's just one reactant
	for react1 in reactants1: reactions.append([react1])
	for react1 in reactants2: reactions.append([react1])
	# now get all combinations of reactant1 x reactant2
	for react1 in reactants1:
		for react2 in reactants2: reactions.append([react1,react2])
		
	filenames_created = []
	
	for reaction in reactions: # go through each of the reactions
		
		# load in the first reactant file
		file1 = reaction[0]
		pdb1 = Molecule2()
		pdb1.load_pdb(file1)
		file1 = os.path.basename(file1)
	
		# by default, assume the second file is empty
		file2 = "empty"
		pdb2 = Molecule2()
		
		if len(reaction) > 1: # but if there are more than one items in the reactants list, load in the second reactant file
			file2 = reaction[1]
			pdb2.load_pdb(file2)
			file2 = os.path.basename(file2)

		# perform the reaction
		react = autoclickchem.OperatorsReact()

		try :
			pdb_list = react.react_molecules(pdb1, pdb2, reactions_to_perform)
			
			if len(pdb_list) > 0: # if there some products, need to save them
				pdb = random.choice(pdb_list)
				
				while True: # this just to make sure you don't repeat the filename
					newindex = random.randint(0,1000000)
					filename = output_dir + str(newindex) + ".mutant.pdb" 
					if not os.path.exists(filename): break
		
				printout = pdb.save_pdb(filename)
				
				if os.path.exists(filename): filenames_created.append(filename)
		except :
			log("\t\t" + 'Error trying to react ' + file1 + ' and ' + file2 + '... skipping... ')

	return filenames_created

###################################################
# Functions to check for drug-like properties	 #
###################################################

def is_numeric(s):
	"""Check to see if a string is numeric
	
	Arguments:
	s -- a string.

	Returns:
	A boolean, True if it is numeric, False otherwise.
	
	"""
	
	try:
		float(s)
		return True
	except ValueError:
		return False

def check_lipinsky_strict(h_bond_donor, h_bond_acceptor, mol_mass, part_coeff):
	"""Check to see if strict Lipinski's Rule of Five is satisfied. Strict means no violations are permitted.
	
	Arguments:
	h_bond_donor -- the number of hydrogen-bond donors, an int.
	h_bond_acceptor -- the number of hydrogen-bond acceptors, an int.
	mol_mass -- the molecular mass, a float.
	part_coeff -- the partition coefficient, a float.

	Returns:
	A boolean, True if the criteria are satisfied, False otherwise.
	
	"""
	
	# number of H-bond donors has to be 5 or less
	limit_h_bond_donor = 5
	
	# number of H-bond acceptors has to be 10 or less
	limit_h_bond_acceptor = 10
	
	# molecular mass has to be 500 daltons or less
	limit_mol_mass = 500
	
	# partion coefficient has to be 5 or less
	limit_part_coeff = 5
	
	if (h_bond_donor > limit_h_bond_donor): return False
	if (h_bond_acceptor > limit_h_bond_acceptor): return False
	if (mol_mass > limit_mol_mass): return False
	if (part_coeff > limit_part_coeff): return False
	return True

def check_lipinsky(h_bond_donor, h_bond_acceptor, mol_mass, part_coeff):
	"""Check to see if strict Lipinski's Rule of Five is satisfied. One violation is permitted.
	
	Arguments:
	h_bond_donor -- the number of hydrogen-bond donors, an int.
	h_bond_acceptor -- the number of hydrogen-bond acceptors, an int.
	mol_mass -- the molecular mass, a float.
	part_coeff -- the partition coefficient, a float.

	Returns:
	A boolean, True if the criteria are satisfied, False otherwise.
	
	"""
	
	# number of H-bond donors has to be 5 or less
	limit_h_bond_donor = 5
	
	# number of H-bond acceptors has to be 10 or less
	limit_h_bond_acceptor = 10
	
	# molecular mass has to be 500 daltons or less
	limit_mol_mass = 500
	
	# partion coefficient has to be 5 or less
	limit_part_coeff = 5
	
	num_violations = 0
	
	if (h_bond_donor > limit_h_bond_donor): num_violations = num_violations + 1
	if (h_bond_acceptor > limit_h_bond_acceptor): num_violations = num_violations + 1
	if (mol_mass > limit_mol_mass): num_violations = num_violations + 1
	if (part_coeff > limit_part_coeff): num_violations = num_violations + 1
	
	if num_violations > 1:
		return False
	else:
		return True

def check_ghose(part_coeff, mol_refr, mol_weight, num_atoms, surf_area): # http://pubs.acs.org/doi/abs/10.1021/cc9800071
	"""Check to see if strict Ghose's Rule is satisfied.
	
	Arguments:
	part_coeff -- the partition coefficient, a float.
	mol_refr -- the molar refractivity, an float.
	mol_weight -- the molecular weight, a float.
	num_atoms -- the number of atoms, an int.
	surf_area -- the surface area, a float.
	
	Returns:
	A boolean, True if the criteria are satisfied, False otherwise.
	
	"""
	
	# partion coefficient has to be between -0.4 and 5.6
	lower_limit_part_coeff = -0.4
	upper_limit_part_coeff = 5.6
	
	# molar refractivity has to be between 40 and 130
	lower_limit_mol_refr = 40
	upper_limit_mol_refr = 130
	
	# molecular weight has to be between 160 and 500 daltons
	lower_limit_mol_weight = 160
	upper_limit_mol_weight = 500
	
	# number of atoms has to be between 20 and 70
	lower_limit_num_atoms = 20
	upper_limit_num_atoms = 70
	
	# polar surface area has to be 140 A2 or less
	limit_surf_area = 140
	if ((part_coeff > upper_limit_part_coeff) or (part_coeff < lower_limit_part_coeff)): return False
	if ((mol_refr > upper_limit_mol_refr) or (mol_refr < lower_limit_mol_refr)): return False
	if ((mol_weight > upper_limit_mol_weight) or (mol_weight < lower_limit_mol_weight)): return False
	if ((num_atoms > upper_limit_num_atoms) or (num_atoms < lower_limit_num_atoms)): return False
	if (surf_area > limit_surf_area): return False
	return True

def check_drug_likeness(pdb_filename):
	"""Check to see if a molecule is drug-like.
	
	Arguments:
	pdb_filename -- the filename of the associated molecular model, a string.
	
	Returns:
	A boolean, True if the molecule is druglike, False otherwise.
	
	"""
	
	global vars
	
	pdb_basename = pdb_filename[:-4]
	pdb_basename = os.path.dirname(pdb_basename) + os.sep + "support" + os.sep + os.path.basename(pdb_basename)

	lipinski = None
	lipinski_strict = None
	ghose = None

	# make the default assumption that the molecule is druglike
	proceed = True

	if os.path.exists(pdb_basename + ".lipinski_pass"): lipinski = True
	if os.path.exists(pdb_basename + ".lipinski_fail"): lipinski = False
	if os.path.exists(pdb_basename + ".lipinski_strict_pass"): lipinski_strict = True
	if os.path.exists(pdb_basename + ".lipinski_strict_fail"): lipinski_strict = False
	if os.path.exists(pdb_basename + ".ghose_pass"): ghose = True
	if os.path.exists(pdb_basename + ".ghose_fail"): ghose = False

	if lipinski == None or lipinski_strict == None or ghose == None:
	
		# run obprop to get some molecular properties
		if not os.path.exists(pdb_basename + ".info"): os.popen(vars['obprop_executable'] + " " + pdb_filename + " > " + pdb_basename + ".info 2>&1", 'r')
		logP = None
		f=open(pdb_basename + ".info", "r")
		for a in f.readlines():
			if "mol_weight" in a: mol_weight = float(a.split()[1])
			if "num_atoms" in a: num_atoms = float(a.split()[1])
			if "logP" in a: logP = float(a.split()[1])
			if "PSA" in a: polar_surface_area = float(a.split()[1])
			if "MR" in a: molecular_refractivity = float(a.split()[1])
		f.close()
		if logP is None:
			log("WARNING: Could not calculate logP for " + pdb_basename + ". Assuming 0.0.")
			logP = 0.0

		# count the number of C and N atoms in the molecule, as well as the C and N atoms that are connected to at least one H atom (these numbers are used to estimate the number of H-bond acceptors and donors)
		O_or_N_counter = int(0)
		O_or_N_with_H_counter = int(0)
		
		molecule_object = Molecule2()
		molecule_object.load_pdb(pdb_filename)
		
		#iterate over all atoms in the molecule
		for atom_index in molecule_object.all_atoms.keys():
			# store object atom in variable atom
			atom = molecule_object.all_atoms[atom_index]
			# deterime if element type of this atom is O or N
			if (atom.element == "O") or (atom.element == "N"):

				# increment counter
				O_or_N_counter += 1

				#iterate over the indices of the atoms that are connected to this atom
				for atomnumber in atom.indecies_of_atoms_connecting:

					# determine element type of connecting atom
					elementtype = molecule_object.all_atoms[atomnumber].element

					# if it is a H atom, increment counter and stop iterating
					if elementtype == "H":
						O_or_N_with_H_counter += 1
						break
	
		lipinski = check_lipinsky(O_or_N_with_H_counter, O_or_N_counter, mol_weight, logP)
		lipinski_strict = check_lipinsky_strict(O_or_N_with_H_counter, O_or_N_counter, mol_weight, logP)
		ghose = check_ghose(logP, molecular_refractivity, mol_weight, num_atoms, polar_surface_area)
	
		if lipinski == True:
			f = open(pdb_basename + ".lipinski_pass", 'w')
			f.close()
		else:
			f = open(pdb_basename + ".lipinski_fail", 'w')
			f.close()
	
		if lipinski_strict == True:
			f = open(pdb_basename + ".lipinski_strict_pass", 'w')
			f.close()
		else:
			f = open(pdb_basename + ".lipinski_strict_fail", 'w')
			f.close()
	
		if ghose == True:
			f = open(pdb_basename + ".ghose_pass", 'w')
			f.close()
		else:
			f = open(pdb_basename + ".ghose_fail", 'w')
			f.close()

	if vars['use_lipinski_filter'] == "TRUE" and lipinski == False: proceed = False
	if vars['use_strict_lipinski_filter'] == "TRUE" and lipinski_strict == False: proceed = False
	if vars['use_ghose_filter'] =="TRUE" and ghose == False: proceed = False
	
	return proceed

#############################################
# Functions for managing the user interface #
#############################################

def program_info():
	"""Get the program version number, etc."""

	program_output = "\nAutoGrow Version 3.1.2\n"
	program_output = program_output + "==================\n"
	program_output = program_output + "If you use AutoGrow 3.1.2 in your research, please cite the following reference:\n"
	program_output = program_output + "Durrant, J. D., ....\n\n"
	
	return program_output

def underline(text):
	"""Print out some text and underline it.
	
	Arguments:
	text -- the text to underline, a string.
	
	"""
	
	log(text)
	toadd = ''
	for t in range(len(text)): toadd = toadd + '='
	log(toadd + "\n")

def show_flag(text, flag, default = "(no default)", prefix = "-"):
	"""Print out information about a command-line parameter.
	
	Arguments:
	text -- the description of the parameter, a string.
	flag -- the name of the parameter, a string.
	default -- the default value of the parameter, a string.
	prefix -- the character to add to the front of the parameter name, a string.
	
	"""
	
	flag = prefix + flag
	
	total_size = 140 # the size of thet whole line
	tag_size = 45 # the width of the initial column
	default_size = 62
	
	while len(flag) < tag_size: flag = flag + " "
	
	# add the default, centered
	flag = flag + default.center(default_size,' ')
	
	#add the description
	flag = flag + text
		
	indent2 = (" "*(tag_size + default_size)) + "   " #(" "*41) + " " + (" "*22) + "   "

	wrapper = textwrap.TextWrapper(total_size,'',indent2)

	log(wrapper.fill(flag))

def help_output():
	"""Print out the help file"""
	
	underline('INTRODUCTION')
	
	wrapper = textwrap.TextWrapper(140)
	log(wrapper.fill('AutoGrow 3.1.2 is a genetic algorithm that attempts to automate the small-molecule inhibitor identification/optimization process. Though no substitute for the medicinal chemist, AutoGrow 3.1.2 can produce chemically feasible drug-like molecules that may supplement the chemist\'s efforts. Version 3 is significantly improved over previous versions. AutoGrow 1.0 and 2.0 simply replace selected hydrogen atoms with molecular fragments, often leading to molecules that are not easy to synthesize; AutoGrow 3 adds fragments using the AutoClickChem algorithm, thus mimicking the many reactions of click chemistry in silico and producing molecules that can actually be synthesized. Additionally, AutoGrow 1.0 and 2.0 add fragments without regard for the drug-like properties of the growing ligands, often leading to compounds that are not necessarily lead like; AutoGrow 3, in contrast, assesses the generated ligands for drug-like properties at each generation, discarding any compounds that do not meet the appropriate criteria.') + "\n")
	
	underline('COMMAND LINE PARAMETERS')

	log("PARAMETER".center(45) + " " + "DEFAULT".center(62) + " " + "DESCRIPTION".center(140-45-62))
	log("---------".center(45) + " " + "-------".center(62) + " " + "-----------".center(140-45-62))

	show_flag("A directory containing the initial population of ligands to optimize, in PDB format.", 'directory_of_source_compounds',"." + os.sep)
	show_flag("A properly formatted directory containing fragments that can be added to the evolving ligands.", 'directory_of_fragments',"{AUTOGROW_HOME}" + os.sep + "fragments" + os.sep + "MW_250" + os.sep + "")
	show_flag("The filename of the receptor into which the evolving ligands will be docked, in PDB format.", 'filename_of_receptor',"." + os.sep)
	show_flag("The number of 'mutant' ligands to create per generation.", 'number_of_mutants', '40')
	show_flag("The number of 'crossover' ligands to create in the first generation. If only one or a few ligands are present in the directory specified by the -directory_of_source_compounds tag, it may not be possible to generate many crossover ligands in the first generation.", 'number_of_crossovers_first_generation','0')
	show_flag("The number of 'mutant' ligands to create in the first generation. If only one or a few ligands are present in the directory specified by the -directory_of_source_compounds tag, it may not be possible to generate many crossover ligands in the first generation.", 'number_of_mutants_first_generation','0')
	show_flag("The number of 'crossover' ligands to create in subsequent generations.", 'number_of_crossovers','0')
	show_flag("The number of best-scoring ligands selected to serve as the founders of the next generation.", 'top_ones_to_advance_to_next_generation','10')
	show_flag("The number of generations to run.", 'num_generations','10')
	show_flag("The program is terminated if a generation runs for longer than the number of seconds specified by this tag.", 'max_seconds_per_generation','1000000')
	show_flag("The x-coordinate of the center of the AutoDock-Vina box used for docking, typically corresponding to the location of the the active site.", 'center_x','0.0')
	show_flag("The y-coordinate of the center of the AutoDock-Vina box used for docking, typically corresponding to the location of the the active site.", 'center_y','0.0')
	show_flag("The z-coordinate of the center of the AutoDock-Vina box used for docking, typically corresponding to the location of the the active site.", 'center_z','0.0')
	show_flag("The number of processors to use, for use on multi-core machines.", 'num_processors','all_processors')
	show_flag("The size of the AutoDock-Vina docking box along the x axis.", 'size_x','30.0')
	show_flag("The size of the AutoDock-Vina docking box along the y axis.", 'size_y','30.0')
	show_flag("The size of the AutoDock-Vina docking box along the z axis.", 'size_z','30.0')
	show_flag("The output directory to which all docking and scoring output files are written.", 'output_dir',"." + os.sep + "output" + os.sep)
	show_flag("The absolute path to the AutoDock Vina executable. The location of this executable can also be specified by modifying the variable definition near the beginning of this script. Example: /home/myname/autodock_vina_1_1_2_linux_x86/bin/vina", 'vina_executable')
	show_flag("The absolute path to the Open Babel bin directory. The location of this directory can also be specified by modifying the variable definition near the beginning of this script. Example: /home/myname/openbabel-2.2.0/bin/", 'openbabel_bin_directory')
	show_flag("AutoGrow attempts to locate the obprop executable based on the directory specified by the -openbabel_bin_directory tag, but, if desired, the location of this executable can be specified explicitly.", 'obprop_executable')
	show_flag("AutoGrow attempts to locate the babel executable based on the directory specified by the -openbabel_bin_directory tag, but, if desired, the location of this executable can be specified explicitly.", 'babel_executable')
	show_flag("The absolute path to the main directory of MGLTools. The location of this directory can also be specified by modifying the variable definition near the beginning of this script. Example: /home/myname/MGLTools-1.5.4/", 'mgltools_directory')
	show_flag("AutoGrow attempts to locate the prepare_ligand4.py script based on the directory specified by the -mgltools_directory tag, but, if desired, the location of this script can be specified explicitly.", 'prepare_ligand4.py')
	show_flag("AutoGrow attempts to locate the prepare_receptor4.py script based on the directory specified by the -mgltools_directory tag, but, if desired, the location of this script can be specified explicitly.", 'prepare_receptor4.py')
	show_flag("AutoGrow attempts to locate the python executable distributed with MGLTools based on the directory specified by the -mgltools_directory tag, but, if desired, the location of this executable can be specified explicitly.", 'mgl_python')
	show_flag("The absolute path to the main AutoClickChem python file, in case you don't want to use the version that ships with AutoGrow.", 'autoclickchem_script', '{AUTOGROW_HOME}' + os.sep + 'support' + os.sep + 'autoclickchem' + os.sep + 'autoclickchem.py')
	show_flag("The absolute path to the main NNScore 1.0 python file, in case you don't want to use the version of NNScore 1.0 that ships with AutoGrow.", 'nn1_script', '{AUTOGROW_HOME}' + os.sep + 'support' + os.sep + 'NNScore' + os.sep + 'NNScore_1.0' + os.sep + 'NNScore.py')
	show_flag("The absolute path to the main NNScore 2.0 python file, in case you don't want to use the version of NNScore 2.0 that ships with AutoGrow. Note that if NNScore 2.0 is used, the -vina_executable command-line parameter must point to a copy of Vina 1.1.2 specifically.", 'nn2_script', '{AUTOGROW_HOME}' + os.sep + 'support' + os.sep + 'NNScore' + os.sep + 'NNScore_2.01' + os.sep + 'NNScore2.01.py')
	show_flag("Which scoring function to use when determining the best-binding ligands of each generation. Acceptable values are \"VINA\", \"NN1\", and \"NN2\" for the AutoDock Vina, NNScore 1.0, and NNScore 2.0 scoring functions, respectively. Note that if NNScore 2.0 is used, the -vina_executable command-line parameter must point to a copy of Vina 1.1.2 specifically.", 'scoring_function', 'VINA')
	show_flag("Only advance ligands judged druglike because they satisfy Lipinski's Rule of Five (Lipinski, C. A., F. Lombardo, et al. 2001. \"Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings.\" Adv Drug Deliv Rev 46: 3-26). Allow for one violation, per Lipinski's recommendations. Acceptable values are \"TRUE\" and \"FALSE\".", 'use_lipinski_filter',"TRUE")
	show_flag("Only advance ligands judged druglike because they satisfy Lipinski's Rule of Five (Lipinski, C. A., F. Lombardo, et al. 2001. \"Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings.\" Adv Drug Deliv Rev 46: 3-26). Ignoring Lipinski's recommendations (which allow for one violation), permit no violations. Acceptable values are \"TRUE\" and \"FALSE\".", 'use_strict_lipinski_filter',"TRUE")
	show_flag("Only advance ligands judged druglike based on the criteria proposed by Ghose et al. (Ghose, A. K., V. N. Viswanadhan, et al. 1999. \"A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. 1. A qualitative and quantitative characterization of known drug databases.\" J Comb Chem 1: 55-68). Acceptable values are \"TRUE\" and \"FALSE\".", 'use_ghose_filter',"FALSE")
	show_flag("Use ligand efficiency (predicted binding affinity divided by the number of heavy atoms in the ligand) instead of the predicted binding affinity alone to assess the goodness of binding. Acceptable values are \"TRUE\" and \"FALSE\".", 'score_by_ligand_efficiency',"FALSE")
	show_flag("Do not use those AutoClickChem reactions that modify compounds without adding new fragments. For example, converting a halide to an azide is considered a modification. On the other hand, merging two compounds via an azide-alkyne cycloaddition is a fragment addition. Acceptable values are \"TRUE\" and \"FALSE\".", 'allow_modification_without_frag_addition',"FALSE")
	show_flag("Sometimes users may wish to ensure that selected atoms from the starting compounds are always present in all ligands subsequently generated. For example, perhaps the starting compounds are known inhibitors, and a certain chemical moiety that is common among them is known to be critical for binding and so should not be modified. These \"core\" atoms are marked by placing a single exclamation mark in the atom name of the associated PDB file. If the -maintain_core flag is set to TRUE, AutoGrow will only generate compounds that contain a certain number of marked, \"core\" atoms. The minimum number of marked atoms permissible is set by the -minimum_core_atoms_required flag. Acceptable values are \"TRUE\" and \"FALSE\".", 'maintain_core',"FALSE")	
	show_flag("This flag is only used if the -maintain_core flag is set to TRUE. AutoGrow will only generate compounds that contain at the very least this number of atoms marked by placing a single exclamation point in the atom name. To protect a deprotonated carboxylate group from being converted to an amide, for example, the user could mark the associated carbon and oxygen atoms by putting an \"!\" in their atom names. The minimum_core_atoms_required tag could then be set to 3, thus requiring that all these carboxylate atoms be present in every generated ligand, thereby preventing AutoGrow from converting this carboxylate group into an amide.", 'minimum_core_atoms_required',"1")	
	show_flag("If you'd like to pass specific commandline parameters to AutoClickChem (e.g., to limit the types of click-chemistry reactions performed), used this parameter. Enclose AutoClickChem parameters in quotes.", 'additional_autoclickchem_parameters',"")	

	underline('EXAMPLE')
	
	wrapper = textwrap.TextWrapper(140)
	log(wrapper.fill('python autogrow_3_1_2.py -number_of_mutants 20 -number_of_crossovers 10 -top_ones_to_advance_to_next_generation 10 -center_x 37.641998 -center_y 25.400999 -center_z 15.80400 -size_x 30.0 -size_y 30.0 -size_z 30. -directory_of_source_compounds ./benzenes/ -directory_of_fragments ./fragments/MW_150/ -filename_of_receptor ./receptors/myreceptor.pdb -output_dir ./output/ -mgltools_directory /.../mgltools_x86_64Linux2_1.5.4 -vina_executable /.../autodock_vina_1_1_2/bin/vina -obprop_executable /.../obprop -openbabel_bin_directory /.../bin -num_generations 10 -use_strict_lipinski_filter TRUE -use_ghose_filter FALSE -allow_modification_without_frag_addition TRUE -maintain_core FALSE') + "\n")
	
	sys.exit(0)

def load_in_commandline_parameters(argv):
	"""Load in the command-line parameters"""

	printout = program_info()

	printout = printout + "(RE)STARTING AUTOGROW 2.0: " + str(datetime.datetime.now())
	printout = printout + "\nUse the -help tag to get detailed help regarding program usage.\n"
	
	vars = define_defaults()
	
	# now get the user-defined variables from the commandline
	allowable_paramters = ['directory_of_source_compounds', 'directory_of_fragments', 'filename_of_receptor', 'number_of_mutants', 'number_of_crossovers_first_generation', 'number_of_mutants_first_generation', 'number_of_crossovers', 'top_ones_to_advance_to_next_generation', 'num_generations', 'max_seconds_per_generation', 'allow_modification_without_frag_addition', 'maintain_core', 'minimum_core_atoms_required', 'center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z', 'output_dir', 'vina_executable', 'openbabel_bin_directory', 'obprop_executable', 'babel_executable', 'mgltools_directory', 'prepare_ligand4.py', 'prepare_receptor4.py', 'mgl_python', 'autoclickchem_script', 'nn1_script', 'nn2_script', 'scoring_function', 'use_lipinski_filter', 'use_strict_lipinski_filter', 'use_ghose_filter', 'score_by_ligand_efficiency', 'num_processors', 'additional_autoclickchem_parameters']
	
	# look for the -help tag
	for t in argv:
		key = t.replace("-","").upper()
		if key.upper() == "HELP": help_output()
	
	parameters_used = []
	
	for t in range(len(argv)-1):
		key = argv[t].replace("-","")
		if key in allowable_paramters:
			parameters_used.append(key)
			value = argv[t+1]
			if is_numeric(value):
				value = float(value)
			vars[key] = value
			argv[t] = ""
			argv[t+1] = ""
	
	# some of these parameters need to be integers, not float
	vars['number_of_mutants'] = int(vars['number_of_mutants'])
	vars['number_of_crossovers_first_generation'] = int(vars['number_of_crossovers_first_generation'])
	vars['number_of_mutants_first_generation'] = int(vars['number_of_mutants_first_generation'])
	
	vars['number_of_crossovers'] = int(vars['number_of_crossovers'])
	vars['top_ones_to_advance_to_next_generation'] = int(vars['top_ones_to_advance_to_next_generation'])
	vars['num_generations'] = int(vars['num_generations'])
	vars['max_seconds_per_generation'] = int(vars['max_seconds_per_generation'])
	vars['minimum_core_atoms_required'] = int(vars['minimum_core_atoms_required'])
	
	vars['num_processors'] = int(vars['num_processors'])
	
	if os.name == "nt" or os.name == "ce": # so it's running under windows. multiprocessing disabled
		vars['num_processors'] = 1
		printout = printout + "\nWARNING: Multiprocessing is disabled on windows machines.\n"


	if not 'number_of_crossovers_first_generation' in parameters_used: vars['number_of_crossovers_first_generation'] = vars['number_of_crossovers'] # so if not specified, set these equal
	if not 'number_of_mutants_first_generation' in parameters_used: vars['number_of_mutants_first_generation'] = vars['number_of_mutants'] # so if not specified, set these equal
	
	# some of the parameters need to be in caps
	vars['use_strict_lipinski_filter'] = vars['use_strict_lipinski_filter'].upper().strip()
	vars['use_lipinski_filter'] = vars['use_lipinski_filter'].upper().strip()
	vars['use_ghose_filter'] = vars['use_ghose_filter'].upper().strip()
	vars['allow_modification_without_frag_addition'] = vars['allow_modification_without_frag_addition'].upper().strip()
	vars['maintain_core'] = vars['maintain_core'].upper().strip()
	vars['score_by_ligand_efficiency'] = vars['score_by_ligand_efficiency'].upper().strip()
	vars['scoring_function'] = vars['scoring_function'].upper().strip()
	
	# Make sure command-line parameter values are acceptable
	if vars['use_strict_lipinski_filter'] != "TRUE" and vars['use_strict_lipinski_filter'] != "FALSE":
		printout = printout + "\nERROR: The -use_strict_lipinski_filter command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['use_strict_lipinski_filter'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['use_lipinski_filter'] != "TRUE" and vars['use_lipinski_filter'] != "FALSE":
		printout = printout + "\nERROR: The -use_lipinski_filter command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['use_lipinski_filter'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['use_ghose_filter'] != "TRUE" and vars['use_ghose_filter'] != "FALSE":
		printout = printout + "\nERROR: The -use_ghose_filter command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['use_ghose_filter'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['allow_modification_without_frag_addition'] != "TRUE" and vars['allow_modification_without_frag_addition'] != "FALSE":
		printout = printout + "\nERROR: The -allow_modification_without_frag_addition command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['allow_modification_without_frag_addition'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['maintain_core'] != "TRUE" and vars['maintain_core'] != "FALSE":
		printout = printout + "\nERROR: The -maintain_core command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['maintain_core'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['score_by_ligand_efficiency'] != "TRUE" and vars['score_by_ligand_efficiency'] != "FALSE":
		printout = printout + "\nERROR: The -score_by_ligand_efficiency command-line parameter must be given a value of \"TRUE\" or \"FALSE\". Value given: \"" + vars['score_by_ligand_efficiency'] + "\"\n"
		log(printout)
		sys.exit(0)
	if vars['scoring_function'] != "VINA" and vars['scoring_function'] != "NN1" and vars['scoring_function'] != "NN2":
		printout = printout + "\nERROR: The -scoring_function command-line parameter must be given a value of \"VINA\", \"NN1\", or \"NN2\". Value given: \"" + vars['scoring_function'] + "\"\n"
		log(printout)
		sys.exit(0)
	
	# Tell the user if any commandline paramters were not used
	keywords_not_used = ""
	for item in argv:
		if item != "": keywords_not_used = keywords_not_used + item + " "
	
	if keywords_not_used.strip() == "": printout = printout + "The following commandline paramters were not used: " + keywords_not_used + "\n" + "\n"
	
	# convert paths to abspath, in case necessary
	vars['directory_of_source_compounds'] = os.path.abspath(vars['directory_of_source_compounds'])
	vars['directory_of_fragments'] = os.path.abspath(vars['directory_of_fragments'])
	vars['filename_of_receptor'] = os.path.abspath(vars['filename_of_receptor'])
	vars['output_dir'] = os.path.abspath(vars['output_dir'])
	vars['mgltools_directory'] = os.path.abspath(vars['mgltools_directory'])
	vars['openbabel_bin_directory'] = os.path.abspath(vars['openbabel_bin_directory'])
	vars['autoclickchem_script'] = os.path.abspath(vars['autoclickchem_script'])
	vars['nn1_script'] = os.path.abspath(vars['nn1_script'])
	vars['nn2_script'] = os.path.abspath(vars['nn2_script'])
	
	# make sure directories end in /
	if vars['directory_of_source_compounds'][-1] != os.sep: vars['directory_of_source_compounds'] = vars['directory_of_source_compounds'] + os.sep
	if vars['directory_of_fragments'][-1] != os.sep: vars['directory_of_fragments'] = vars['directory_of_fragments'] + os.sep
	if vars['output_dir'][-1] != os.sep: vars['output_dir'] = vars['output_dir'] + os.sep
	if vars['mgltools_directory'][-1] != os.sep: vars['mgltools_directory'] = vars['mgltools_directory'] + os.sep
	if vars['openbabel_bin_directory'][-1] != os.sep: vars['openbabel_bin_directory'] = vars['openbabel_bin_directory'] + os.sep
	
	# find other mgltools-related scripts
	if vars['prepare_ligand4.py'] == "": vars['prepare_ligand4.py'] = vars['mgltools_directory'] + 'MGLToolsPckgs' + os.sep + 'AutoDockTools' + os.sep + 'Utilities24' + os.sep + 'prepare_ligand4.py'
	if vars['prepare_receptor4.py'] == "": vars['prepare_receptor4.py'] = vars['mgltools_directory'] + 'MGLToolsPckgs' + os.sep + 'AutoDockTools' + os.sep + 'Utilities24' + os.sep + 'prepare_receptor4.py'
	if vars['mgl_python'] == "": vars['mgl_python'] = vars['mgltools_directory'] + 'bin' + os.sep + 'pythonsh'
	
	# find other babel executables
	if vars['obprop_executable'] == "": vars['obprop_executable'] = vars['openbabel_bin_directory'] + "obprop"
	if vars['babel_executable'] == "": vars['babel_executable'] = vars['openbabel_bin_directory'] + "babel"
	
	# check to see if above files don't exist. If they don't, check if '.exe' versions exist (Windows)
	if not os.path.exists(vars['obprop_executable']) and os.path.exists(vars['obprop_executable'] + ".exe"): vars['obprop_executable'] = vars['obprop_executable'] + ".exe"
	if not os.path.exists(vars['babel_executable']) and os.path.exists(vars['babel_executable'] + ".exe"): vars['babel_executable'] = vars['babel_executable'] + ".exe"
	
	# convert path names with spaces if this is windows
	if os.name == "nt" or os.name == "ce": # so it's running under windows. multiprocessing disabled
		if " " in vars['directory_of_source_compounds']: vars['directory_of_source_compounds'] = '"' + vars['directory_of_source_compounds'] + '"'
		if " " in vars['directory_of_fragments']: vars['directory_of_fragments'] = '"' + vars['directory_of_fragments'] + '"'
		if " " in vars['filename_of_receptor']: vars['filename_of_receptor'] = '"' + vars['filename_of_receptor'] + '"'
		if " " in vars['output_dir']: vars['output_dir'] = '"' + vars['output_dir'] + '"'
		if " " in vars['mgltools_directory']: vars['mgltools_directory'] = '"' + vars['mgltools_directory'] + '"'
		if " " in vars['openbabel_bin_directory']: vars['openbabel_bin_directory'] = '"' + vars['openbabel_bin_directory'] + '"'
		if " " in vars['autoclickchem_script']: vars['autoclickchem_script'] = '"' + vars['autoclickchem_script'] + '"'
		if " " in vars['nn1_script']: vars['nn1_script'] = '"' + vars['nn1_script'] + '"'
		if " " in vars['nn2_script']: vars['nn2_script'] = '"' + vars['nn2_script'] + '"'
		if " " in vars['prepare_ligand4.py']: vars['prepare_ligand4.py'] = '"' + vars['prepare_ligand4.py'] + '"'
		if " " in vars['prepare_receptor4.py']: vars['prepare_receptor4.py'] = '"' + vars['prepare_receptor4.py'] + '"'
		if " " in vars['mgl_python']: vars['mgl_python'] = '"' + vars['mgl_python'] + '"'
		if " " in vars['obprop_executable']: vars['obprop_executable'] = '"' + vars['obprop_executable'] + '"'
		if " " in vars['babel_executable']: vars['babel_executable'] = '"' + vars['babel_executable'] + '"'

	
	# output the paramters used
	printout = printout + "\nPARAMETERS" + "\n"
	printout = printout + "==========" + "\n"
	for key in allowable_paramters: printout = printout + key + ": " + str(vars[key]) + "\n"
	printout = printout + "" + "\n"
	
	# make sure scripts and executables exist
	if not os.path.exists(vars['prepare_ligand4.py']) and not os.path.exists(vars['prepare_ligand4.py'].replace('"','')):
		printout = printout + "\nERROR: Could not find prepare_ligand4.py at " + vars['prepare_ligand4.py'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['prepare_receptor4.py']) and not os.path.exists(vars['prepare_receptor4.py'].replace('"','')):
		printout = printout + "\nERROR: Could not find prepare_receptor4.py at " + vars['prepare_receptor4.py'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['mgl_python']) and not os.path.exists(vars['mgl_python'].replace('"','')):
		printout = printout + "\nERROR: Could not find pythonsh at " + vars['mgl_python'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['obprop_executable']) and not os.path.exists(vars['obprop_executable'].replace('"','')):
		printout = printout + "\nERROR: Could not find pythonsh at " + vars['obprop_executable'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['babel_executable']) and not os.path.exists(vars['babel_executable'].replace('"','')):
		printout = printout + "\nERROR: Could not find pythonsh at " + vars['babel_executable'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['autoclickchem_script']) and not os.path.exists(vars['autoclickchem_script'].replace('"','')):
		printout = printout + "\nERROR: Could not find " + os.path.basename(vars['autoclickchem_script']) + " at " + vars['autoclickchem_script'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['nn1_script']) and not os.path.exists(vars['nn1_script'].replace('"','')):
		printout = printout + "\nERROR: Could not find " + os.path.basename(vars['nn1_script']) + " at " + vars['nn1_script'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['nn2_script']) and not os.path.exists(vars['nn2_script'].replace('"','')):
		printout = printout + "\nERROR: Could not find " + os.path.basename(vars['nn2_script']) + " at " + vars['nn2_script'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['obprop_executable']) and not os.path.exists(vars['obprop_executable'].replace('"','')):
		printout = printout + "\nERROR: Could not find " + os.path.basename(vars['obprop_executable']) + " at " + vars['obprop_executable'] + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['babel_executable']) and not os.path.exists(vars['babel_executable'].replace('"','')):
		printout = printout + "\nERROR: Could not find " + os.path.basename(vars['babel_executable']) + " at " + vars['babel_executable'] + "\n"
		log(printout)
		sys.exit(0)
	
	# now check to make sure receptors and ligands exist
	ligands_list = glob.glob(vars['directory_of_source_compounds'] + "*.pdb")
	ligands_list.extend(glob.glob(vars['directory_of_source_compounds'] + "*.PDB"))
	
	if len(ligands_list) == 0:
		printout = printout + "\nERROR: There are no PDB-formatted ligands in the directory \"" + vars['directory_of_source_compounds'] + "\"." + "\n"
		log(printout)
		sys.exit(0)
	if not os.path.exists(vars['filename_of_receptor']):
		printout = printout + "\nERROR: There receptor file does not exist: \"" + vars['filename_of_receptor'] + "\"." + "\n"
		log(printout)
		sys.exit(0)
	
	# if the output directory doesn't exist, then make it
	if not os.path.exists(vars['output_dir']):
		os.makedirs(vars['output_dir'])
		os.makedirs(vars['output_dir'] + "support")
	
	return vars, printout

def log(thestring):
	"""Print text to the screen and to a log file

	Arguments:
	thestring -- the text to print, a string.
	
	"""

	global log_text_file

	print thestring

	# Ensure log_text_file variable is defined. If so, write to file.
	try:
		log_text_file
	except NameError:
		log_text_file = None
	
	# Test whether variable is defined to be None
	if not log_text_file is None: log_text_file.write(str(thestring) + "\n")


###########################################
# Functions for creating molecular models #
###########################################

def make_mutants(directory):
	"""Make mutant compounds in the specified directory
	
	Arguments:
	directory -- the directory name, a string.
	
	"""
	
	log("\nCreating new ligands...")
	
	global vars
	
	# first check if all the mutants have been made
	if len(glob.glob(directory + "*.mutant.pdb")) < vars['number_of_mutants_current_generation']: # so you need to make more mutants
	
		# get a list of the fragments and ligands
		frag_dirs = glob.glob(vars['directory_of_fragments'] + "*")
		fragments = {}
		for frag_dir in frag_dirs:
			fragments[frag_dir] = glob.glob(frag_dir + os.sep + "*.pdb")

		ligands = glob.glob(vars['directory_of_source_compounds'] + "*.pdb")
		
		show_initial_message = False
		
		vars_to_pass = {}
		vars_to_pass['ligands'] = ligands
		vars_to_pass['frag_dirs'] = frag_dirs
		vars_to_pass['fragments'] = fragments
		vars_to_pass['directory'] = directory

		while len(glob.glob(directory + "*.mutant.pdb")) < vars['number_of_mutants_current_generation']: # keep generating new mutants until you've passed the allowable level
			multi_threading(range(100), run_autoclick_multithread, vars_to_pass) # so, performing 100 ligands at a time

		# note that this is going to overshoot, and generate too many ligands. These will be deleted later

def make_crossovers(directory):
	"""Make crossover compounds in the specified directory
	
	Arguments:
	directory -- the directory name, a string.
	
	"""
	
	global vars
	
	if len(glob.glob(directory + "*.crossover.pdb")) < vars['number_of_crossovers_current_generation']:
	
		# get a list of the ligands
		ligands = glob.glob(vars['directory_of_source_compounds'] + "*.pdb")
		
		vars_to_pass = {}
		vars_to_pass['ligands'] = ligands
		vars_to_pass['directory'] = directory

		while len(glob.glob(directory + "*.crossover.pdb")) < vars['number_of_crossovers_current_generation']: # keep generating new crossovers until you've passed the allowable level
			multi_threading(range(100), run_ligmerge_multithread, vars_to_pass) # so, performing 100 ligands at a time
			# this is a lousy solution!!! Need to pass variables between threads!!! ******
			
		# note that this is going to overshoot, and generate too many ligands. These will be deleted later

##########################################
# Functions for docking molecular models #
##########################################

def convert_ligand_pdb_files_to_pdbqt(directory):
	"""Convert the ligands of a given directory from pdb to pdbqt format
	
	Arguments:
	directory -- the directory name, a string.
	
	"""
	
	count = 0

	count = count + 1
	if count > 10000:
		log("ERROR: I've tried 10,000 times to convert all the ligand PDB files in " + directory + " to the PDBQT format. Aborting program...")
		sys.exit(0)
	
	# check to make sure each of the pdb files has been converted to pdbqt
	need_to_make_pdbqt = []
	for filename in glob.glob(directory + "*.pdb"):
		if not os.path.exists(filename + "qt"): need_to_make_pdbqt.append(filename)
	
	# create a file to run the pdbqt making (from command line or through qsub... should eventually be two options...)
	jobs = []
	for filename in need_to_make_pdbqt:
		convert_pdb_to_pdbqt_acceptable_format(filename)
		jobs.append(vars['mgl_python'] + " " + vars['prepare_ligand4.py'] + " -g -l " + filename + " -o " + filename + "qt")
	multi_threading(jobs, execute_command, {'seconds_per_job': 100000, 'log_filename':directory + "support" + os.sep + "1.make_pdbqt"}) # give it all the time it needs to make these receptor pdbqt files


def convert_pdb_to_pdbqt_acceptable_format(filename): 
	"""Make sure a PDB file is properly formatted for conversion to pdbqt
	
	Arguments:
	directory -- the filename, a string.
	
	"""
	
	# read in the file
	output_lines = []
	f = open(filename,'r')
	for line in f.readlines():
		line = line.strip()
		if line[:5] == "ATOM " or line[:7] == "HETATM ":
			#fix things like atom names with two letters
			first = line[:11]
			middle = line[11:17].upper().strip()
			last = line[17:]
			
			middle_firstpart = ""
			middle_lastpart = middle
			
			for i in range(len(middle_lastpart)):
				if middle_lastpart[:1] in string.uppercase:
					middle_firstpart = middle_firstpart + middle_lastpart[:1]
					middle_lastpart = middle_lastpart[1:]
				else: break # you reached the first number
			
			# now if there are more than two letters in middle_firstpart, keep just two
			if len(middle_firstpart) > 2:
				middle_lastpart = middle_firstpart[2:] + middle_lastpart
				middle_firstpart = middle_firstpart[:2]
			
			if not (middle_firstpart == "BR" or middle_firstpart == "ZN" or middle_firstpart == "FE" or middle_firstpart == "MN" or middle_firstpart == "CL" or middle_firstpart == "MG"):
				# so just keep the first letter for the element part of the atom name
				middle_lastpart = middle_firstpart[1:] + middle_lastpart
				middle_firstpart = middle_firstpart[:1]
			
			middle = middle_firstpart.rjust(3) + middle_lastpart.ljust(3)
			
			line = first + middle + last
			
			# make sure all parts of the molecule belong to the same chain and resid
			line = line[:17] + "LIG X 999" + line[26:]
			
			output_lines.append(line)
		else:
			output_lines.append(line)
			
	f.close()
	f = open(filename,'w')

	for line in output_lines:
		f.write(line + "\n")

	f.close()

def convert_receptor_pdb_files_to_pdbqt():
	"""Make sure a PDB file is properly formatted for conversion to pdbqt
	
	Arguments:
	directory -- the filename, a string.
	
	"""
	
	global vars
	
	count = 0
	
	while not os.path.exists(vars['filename_of_receptor'] + "qt"): 
	
		count = count + 1
		if count > 10000:
			log("ERROR: I've tried 10,000 times to convert the file \"" + vars['filename_of_receptor'] + "\" to the PDBQT format. Aborting program...")
			sys.exit(0)
		
		# make sure the receptors have been converted to PDBQT. If not, do the conversion.
		receptors = [vars['filename_of_receptor']] 
		need_to_covert_receptor_to_pdbqt = []
		for filename in receptors: 
			if not os.path.exists(filename + "qt"): need_to_covert_receptor_to_pdbqt.append(filename)
		
		# create a file to run the pdbqt making (from command line or through qsub... should eventually be two options...)
		jobs = []
		for filename in need_to_covert_receptor_to_pdbqt: jobs.append(vars['mgl_python'] + " " + vars['prepare_receptor4.py'] + " -r " + filename + " -o " + filename + "qt")
		
		multi_threading(jobs, execute_command, {'seconds_per_job': 30, 'log_filename': os.path.dirname(vars['filename_of_receptor']) + os.sep + "make_pdbqt"}) # should only take 30 seconds to make a pdbqt file

def dock_compounds(directory):
	"""Dock the ligand pdbqt files in a given directory using AutoDock Vina
	
	Arguments:
	directory -- the filename, a string.
	
	"""
	
	count = 0
	
	count = count + 1
	if count > 10000:
		log("ERROR: I've tried 10,000 times to dock the PDBQT files of " + directory + ". Aborting program...")
		sys.exit(0)

	receptors = [vars['filename_of_receptor'] + "qt"] 

	# find ligands that have not been docked
	need_to_dock = []
	for filename in glob.glob(directory + "*.pdbqt"):
		if not os.path.exists(filename + ".vina"): need_to_dock.append(filename)
	
	# do the docking of the needed ligands
	jobs = []
	for receptor_filename in receptors:
		for lig_filename in need_to_dock:
			torun = vars['vina_executable'] + " --center_x " + str(vars['center_x']) + " --center_y " + str(vars['center_y']) + " --center_z " + str(vars['center_z']) + " --size_x " + str(vars['size_x']) + " --size_y " + str(vars['size_y']) + " --size_z " + str(vars['size_z']) + " --receptor " + receptor_filename + " --ligand " + lig_filename + " --out " + lig_filename + ".vina --cpu 1"
			if vars['scoring_function'] == "NN1": # so also evaluate with NNScore 1.0 scoring function
				torun = torun + "; " + sys.executable + " " + vars['nn1_script'] + " -receptor " + receptor_filename + " -vina_output " + lig_filename + ".vina -networks_dir " + os.path.dirname(vars['nn1_script']) + os.sep + "networks" + os.sep + "top_3_networks" + os.sep + " > " + lig_filename + ".vina.nn1"
			if vars['scoring_function'] == "NN2": # so also evaluate with NNScore 2.01 scoring function
				torun = torun + "; " + sys.executable + " " + vars['nn2_script'] + " -receptor " + receptor_filename + " -ligand " + lig_filename + ".vina -vina_executable " + vars['vina_executable'] + " > " + lig_filename + ".vina.nn2"
			jobs.append(torun)

	multi_threading(jobs, execute_command, {'seconds_per_job': 600, 'log_filename': directory + "support" + os.sep + "2.do_docking"}) # give it only ten minutes per docking

########################
# Delete extra ligands #
########################

def delete_bad_ligands(current_generation_dir):
	"""Delete ligands in a directory that don't meet user specified criteria
	
	Arguments:
	current_generation_dir -- the name of the directory, a string.
	
	"""
	
	global vars

	log("\nThere are " + str(len(glob.glob(current_generation_dir + "*.pdb"))) + " ligands in the current directory.")
	
	# now delete mutants that have steric problems
	steric_show_message = False
	for filename in glob.glob(current_generation_dir + "*steric*.pdb"):
		delete_all_associated_files(filename)
		log("Steric clash: deleting ligand " + filename)

	# delete ligands that correspond to empty files
	for filename in glob.glob(current_generation_dir + "*.pdb"):
		if os.path.getsize(filename) == 0:
			log("Ligand file contains no data: deleting ligand " + filename)
			delete_all_associated_files(filename)

	for filename in glob.glob(current_generation_dir + "*.mutant.pdb"):
		check_time_limit_reached()
		
		tmp_f = open(filename,'r')
		line = tmp_f.readline()
		line = tmp_f.readline() # the second line contains info about the files used to create this mutant
		tmp_f.close()
		
		line = line.split(":")
		line = line[1].strip()
		lines = line.split(";")
		
		if len(lines) == 1: # because otherwise it must be drawing parents form both the fragment library and the other source
			frag_dir = vars['directory_of_fragments']
			if frag_dir[-1:] == os.sep: frag_dir = frag_dir[:-1]
			
			if os.path.dirname(lines[0].strip())[:len(frag_dir.strip())] == frag_dir.strip():
				log("Ligand derived only from fragments: deleting ligand " + filename)
				delete_all_associated_files(filename)
				
	shown_delete_general_message = False
	if vars['allow_modification_without_frag_addition'] == "FALSE": # so throw out ones that only modified the ligand, without adding a fragment to it. This will probably never happen.
		bad_ligand_show_message = False
		for filename in glob.glob(current_generation_dir + "*.mutant.pdb"):
			check_time_limit_reached()
			
			tmp_f = open(filename,'r')
			line = tmp_f.readline()
			line = tmp_f.readline() # the second line contains info about the files used to create this mutant
			tmp_f.close()
			
			# the way autoclickchem is currently formatted, if a PDB was derived from two PDBS, it will have a ";" in its second line. Be sure not to change this formatting!
			if not ";" in line: # so it was just modifying the ligand, since fragment-only modifications have already been deleted
				if shown_delete_general_message == False:
					log("\nRemoving undesirable ligands from the current generation...")
					shown_delete_general_message = True
				if bad_ligand_show_message == False:
					log("\tRemoving ligands that were modified without fragment addition...")
					bad_ligand_show_message = True
				log("\t\tDeleting ligand " + filename)
				delete_all_associated_files(filename)
	
	# if user requested, make sure all generated ligands have the initially specified ligand
	if vars['maintain_core'] == 'TRUE':
		bad_ligand_show_message = False
		for filename in glob.glob(current_generation_dir + "*.pdb"):
			check_time_limit_reached()
			
			tmp_f = open(filename,'r')
			count = 0
			for line in tmp_f.readlines():
				if line[:7] == "HETATM " or line[:5] == "ATOM ":
					if "!" in line[11:16]: # so there's an atom marked with ! in its name
						count = count + 1
			tmp_f.close()
			
			if count < vars['minimum_core_atoms_required']: # so there is no atom marked with ! in its name
				if shown_delete_general_message == False:
					log("\nRemoving undesirable ligands from the current generation...")
					shown_delete_general_message = True
				if bad_ligand_show_message == False:
					log("\tRemoving ligands that do not contain at least " + str(vars['minimum_core_atoms_required']) + " marked (\"core\") atoms...")
					bad_ligand_show_message = True
				log("\t\tDeleting ligand " + filename + " because it contains only " + str(count) + " marked (\"core\") atoms.")
				delete_all_associated_files(filename)
		if bad_ligand_show_message == True: log("\t\tThere are now " + str(len(glob.glob(current_generation_dir + "*.pdb"))) + " ligands in the current directory.")
	
	shown_initial_message = False
	
	# check to see if any of the created ligands are not druglike. If so, delete them.
	# this can be very time consuming, so it's good to use multi-threading
	multi_threading(glob.glob(current_generation_dir + "*.pdb"), druglike_multithread, {})

	remove_duplicates(current_generation_dir)
	
	# if after all this filtering there are still more than the required number of mutants and crossovers, delete the extra ones
	while len(glob.glob(current_generation_dir + "*.crossover.pdb")) > vars['number_of_crossovers_current_generation']:
		# pick a random ligand
		filename = random.choice(glob.glob(current_generation_dir + "*.crossover.pdb"))
		
		if len(glob.glob(filename[:-3] + "*")) == 1: # just the pdb file
			log("Too many crossovers: deleting file " + filename + "...")
			delete_all_associated_files(filename)

	while len(glob.glob(current_generation_dir + "*.mutant.pdb")) > vars['number_of_mutants_current_generation']:
		# pick a random ligand
		filename = random.choice(glob.glob(current_generation_dir + "*.mutant.pdb"))
		
		if len(glob.glob(filename[:-3] + "*")) == 1: # just the pdb file
			log("Too many mutants: deleting file " + filename + "...")
			delete_all_associated_files(filename)

def remove_duplicates(current_directory):
	"""Delete ligands in a directory that are duplicates
	
	Arguments:
	current_directory -- the name of the directory, a string.
	shown_delete_general_message -- controls whether or not messages are displayed
	
	"""
	
	# first first step is to identify the smiles strings of any missing ligands
	# best to do this with multithreading
	jobs = []
	for filename in glob.glob(current_directory + '*.pdb'):
		can_filename = os.path.dirname(filename) + os.sep + "support" + os.sep + os.path.basename(filename) + ".can"
		if not os.path.exists(can_filename): jobs.append(vars['babel_executable'] + ' -ipdb ' + filename + ' -ocan ' + can_filename)
	if len(jobs) > 0: multi_threading(jobs, execute_command, {'seconds_per_job': 10})

	# now get a list of all the smiles strings
	smiles = {}
	for filename in glob.glob(current_directory + '*.pdb'):
		can_filename = os.path.dirname(filename) + os.sep + "support" + os.sep + os.path.basename(filename) + ".can"
		f = open(can_filename, 'r')
		smile = f.read().split("\t")[0]
		f.close()
		
		if not smile in smiles.keys(): smiles[smile] = []
		smiles[smile].append(filename)
	
	# now go through each and remove duplicates
	for smile in smiles.keys():
		if len(smiles[smile]) > 1: # so there are multiple ones
			
			one_has_been_removed = False
			
			# if any of these have already been processed, remove it from the list
			for filename in smiles[smile]:
				if len(glob.glob(filename + "*")) > 1:
					log("Keeping " + filename + "...")
					smiles[smile].remove(filename)
					one_has_been_removed = True
					break # break so this sparing can happen only one
			
			if one_has_been_removed == False:
				log("Keeping " + smiles[smile][0] + "...")
				smiles[smile].pop(0) # just remove the first one
			
			# now delete the remaining files
			for filename in smiles[smile]:
				log("Deleting " + filename + ", which is a duplicate...")
				delete_all_associated_files(filename)

def delete_no_pdbqt(current_generation_dir):
	"""Delete ligands that cannot be converted to pdbqt
	
	Arguments:
	current_generation_dir -- the name of the directory containing the ligands, a string.
	
	"""
	
	if len(glob.glob(current_generation_dir + "*.pdbqt")) == 0:
		log("ERROR: No pdbqt files were generated, suggesting a significant problem. Most likely, AutoGrow is unable to successfully execute the prepare_ligand4.py MGLTools script.")
		sys.exit(0)
	for filename in glob.glob(current_generation_dir + "*.pdb"):
		if not os.path.exists(filename + "qt"): # so this pdbqt file didn't exist
			log("PDBQT not generated: Deleting " + os.path.basename(filename) + "...")
			delete_all_associated_files(filename)
	
def delete_no_docked(current_generation_dir):
	"""Delete ligands that could not be docked
	
	Arguments:
	current_generation_dir -- the name of the directory containing the ligands, a string.
	
	"""

	if len(glob.glob(current_generation_dir + "*.vina")) == 0:
		log("ERROR: None of the dockings were successful, suggesting a significant problem. Most likely, AutoGrow is unable to successfully run the AutoDock vina executable.")
		sys.exit(0)
	for filename in glob.glob(current_generation_dir + "*.pdb"):
		if not os.path.exists(filename + "qt.vina"): # so this pdbqt file didn't exist
			log("Docking unsuccessful: Deleting " + os.path.basename(filename) + "...")
			delete_all_associated_files(filename)

def delete_all_associated_files(pdb_filename):
	"""Delete files associated with a compound
	
	Arguments:
	pdb_filename -- the filename of the compounds, a string
	
	"""

	toremove = [pdb_filename]
	toremove.extend(glob.glob(pdb_filename[:-3] + "*"))
	toremove.extend(glob.glob(os.path.dirname(pdb_filename) + os.sep + "support" + os.sep + os.path.basename(pdb_filename)[:-3] + "*"))
	
	for todel in toremove:
		if os.path.exists(todel): os.remove(todel)

#############################
# save and load smiles data #
#############################

def save_smiles_and_druglike():
	"""Save the smiles and druglike-properties information that has been calculated so it need not be recalculated if the program restarts."""
	
	global vars, smiles, druglike, smiles_used_in_previous_generations

	if not os.path.exists(vars['output_dir'] + 'support' + os.sep): os.mkdir(vars['output_dir'] + 'support' + os.sep)

	fh = open(vars['output_dir'] + 'support' + os.sep + 'smiles_and_druglike.pickle', 'w')
	cPickle.dump([smiles, druglike, smiles_used_in_previous_generations],fh)
	fh.close()

def load_smiles_and_druglike():
	"""Load smiles and druglike-properties data saved from a previous run"""
	
	global smiles, druglike, smiles_used_in_previous_generations
	
	# first check to see if the pickle file exists
	if os.path.exists(vars['output_dir'] + 'support' + os.sep + 'smiles_and_druglike.pickle'):
		log("Loading previously calculated SMILES-string and druglike data for ligands already generated...\n")
		try:
			fh = open(vars['output_dir'] + 'support' + os.sep + 'smiles_and_druglike.pickle', 'r')
			tmp = cPickle.load(fh)
			fh.close()
			smiles = tmp[0]
			druglike = tmp[1]
			smiles_used_in_previous_generations = tmp[2]
		except:
			log("\tERROR: Could not load data from the smiles_and_druglike.pickle file!")


######################
# manage directories #
######################

def get_current_generation_directory():
	"""Get the current generation directory
	
	Returns:
	A string, the name of the current generation directory.

	"""

	current_generation_dir = ""
	if not os.path.exists(vars['output_dir'] + "generation1" + os.sep): # check to see if the first generation directory exists
		current_generation_dir = vars['output_dir'] + "generation1" + os.sep
		os.makedirs(current_generation_dir)
		os.makedirs(current_generation_dir + "support")
	else: # it does exist, so find the latest generation
		# gen all the generations that have been written
		gens = glob.glob(vars['output_dir'] + "generation*")
		for i in range(len(gens)):
			gen = gens[i]
			gen = gen.replace(vars['output_dir'],'')
			gens[i] = gen
		
		# identify the latest generation
		do_continue = True
		t = 0
		while do_continue == True:
			t = t + 1
			if "generation" + str(t) not in gens:
				do_continue = False
		
		current_generation_dir = vars['output_dir'] + "generation" + str(t-1) + os.sep
		if t-2 != 0: vars['directory_of_source_compounds'] = vars['output_dir'] + "generation" + str(t-2) + os.sep + "best_ligands" + os.sep
	return current_generation_dir

def identify_top_compounds(directory):
	"""Identify the ligands from the current generational directory that will be selected as the founding members of the next generation
	
	Arguments:
	directory -- the name of the current generational directory containing the ligands, a string.
	
	Returns:
	A list of strings, the names of the files describing the top-scoring molecules.
	
	"""
	
	vina_outputs = glob.glob(directory + "*.vina")
	results = {}
	vina_pose_to_use = {} # {filename: pose}
	for vina_output in vina_outputs:
		if vars['scoring_function'] == "VINA":
			f = open(vina_output,'r')
			score = f.readline()
			score = f.readline()
			while "  " in score: score = score.replace("  "," ")
			score = score.split(" ")
			score = float(score[3])
			f.close()
			vina_pose_to_use[vina_output] = 1
		elif vars['scoring_function'] == "NN1":
			nn1_output = vina_output + ".nn1"
			f = open(nn1_output,'r')
			for line in f.readlines():
				if "Best score:" in line:
					tmp = line.split(" ")
					score = float(tmp[2])
					tmp = line.replace(")",'')
					tmp = tmp.split(', MODEL ')
					vina_pose_to_use[vina_output] = int(tmp[1])
			f.close()
		elif vars['scoring_function'] == "NN2":
			nn1_output = vina_output + ".nn2"
			f = open(nn1_output,'r')
			
			while True:
				line = f.readline()
				if len(line) == 0: break
				if "When the poses were ranked by the average of the 20 network scores" in line:
					line = f.readline()
					tmp = line.split(" ")
					vina_pose_to_use[vina_output] = int(tmp[9])
					line = f.readline()
					tmp = line.split(" ")
					score = float(tmp[0])
			f.close()
		
		if vars['score_by_ligand_efficiency'] == "TRUE": # so you need to count the number of non-hydrogen atoms and divide the score by that.
			pdb = Molecule2()
			pdb.load_pdb(vina_output[:-7])
			
			count = 0
			
			for index in pdb.all_atoms.keys():
				atom = pdb.all_atoms[index]
				if atom.element != "H": count = count + 1
			
			score = score / float(count)
		
		results[vina_output] = score
	
	results = sorted(results.iteritems(), key=operator.itemgetter(1))
		
	if vars['scoring_function'] == "VINA": results = results[:vars['top_ones_to_advance_to_next_generation']]
	elif vars['scoring_function'] != "VINA": # because for both NNScore functions, a higher score is better
		results = results[-vars['top_ones_to_advance_to_next_generation']:]
		results.reverse()
	
	new_results = []
	new_scores = []
	for result in results:
		new_results.append(result[0][:-7])
		new_scores.append(result[1])

	# save scores
	thegen = os.path.basename(directory[:-1])
	underline = ""
	for t in range(len(thegen)): underline = underline + "="
	
	f = open(vars['output_dir'] + "scores", 'a')
	f.write(thegen + "\n")
	f.write(underline + "\n\n")
	
	if vars['scoring_function'] == "VINA": docked_pose_filename_tag = 'best_vina_pose'
	elif vars['scoring_function'] == "NN1": docked_pose_filename_tag = 'best_nn1_pose'
	elif vars['scoring_function'] == "NN2": docked_pose_filename_tag = 'best_nn2_pose'
	
	for t in range(len(new_results)):
		f.write(str(new_scores[t]) + "\t" + new_results[t].replace(vars['output_dir'],"") + "\n")
	f.write("\n")
	f.close()
	
	if not os.path.exists(directory + "best_ligands" + os.sep): os.makedirs(directory + "best_ligands" + os.sep)
	if not os.path.exists(directory + "best_ligands_docked" + os.sep): os.makedirs(directory + "best_ligands_docked" + os.sep)

	# Save the compound to a separate directory
	for result in new_results:
		shutil.copy2(result,directory + "best_ligands" + os.sep)
		extract_a_vina_pose(result + "qt.vina", vina_pose_to_use[result + "qt.vina"], directory + "best_ligands_docked" + os.sep + os.path.basename(result)[:-4] + "." + docked_pose_filename_tag + ".pdb")
		log("\tTop ligand: " + os.path.basename(result))
	
	return new_results

def extract_a_vina_pose(vina_output_filename, pose_number, pdb_filename):
	"""Loads in a vina output file and saves the specified frame to a PDB file.
	
	Arguments:
	vina_output_filename -- the name of the vina output file, a string.
	pose_number -- the pose to save, and integer.
	pdb_filename -- the name of the pdb file that will be written, a string.
	
	"""
	
	f = open(vina_output_filename, 'r')
	index = 1
	pdbs = {}
	this_pdb = []
	for line in f.readlines():
		line = line.strip()
		if line[:4] == "ATOM" or line[:6] == "HETATM": this_pdb.append(line)
		if line[:6] == "ENDMDL":
			pdbs[index] = this_pdb[:]
			index = index + 1
			this_pdb = []
	f.close()
	
	one_to_use = pdbs[pose_number]
	
	f = open(pdb_filename, 'w')
	for line in one_to_use:
		f.write(line + "\n")
	f.close()
	
def make_next_generation_directory(current_directory):
	"""Create the next generational directory
	
	Arguments:
	current_directory -- the name of the current generational directory, a string.
	
	Returns:
	A list of strings, the names of the files describing the top-scoring molecules.
	
	"""
	
	global vars, top_compounds
	
	# make the next generation
	current_index = int(os.path.basename(current_directory[:-1]).replace("generation",""))
	newdir = vars['output_dir'] + "generation" + str(current_index+1) + os.sep
	if not os.path.exists(newdir):
		os.makedirs(newdir)
		os.makedirs(newdir + "support")
	
	# get a list of the best compounds from the previous generation
	if top_compounds != []:
		for cmp in top_compounds:
			nme = os.path.basename(cmp)
			nme = nme[:-3]
			for filename in glob.glob(current_directory + nme + "*"):
				shutil.copy2(filename, newdir)

			for filename in glob.glob(current_directory + os.sep + "support" + os.sep + nme + "*"):
				shutil.copy2(filename, newdir + os.sep + "support" + os.sep)
			
def check_time_limit_reached():
	"""Check if it has taken too long to process the current generation"""
	
	global vars, start_time_current_generation
	
	time_passed = time.time() - start_time_current_generation
	
	if time_passed > vars['max_seconds_per_generation']:
		log("\nERROR: It has been " + str(time_passed) + " seconds since the current generation began. You requested program termination after " + str(vars['max_seconds_per_generation']) + " seconds.\n")
		log("Program terminated\n")
		sys.exit(0)

def set_start_time_current_generation(logit = True):
	"""Mark the time the current generation starts."""
	
	global start_time_current_generation
	start_time_current_generation = time.time()
	if logit == True: log("The current generation began at " + time.asctime(time.localtime(start_time_current_generation)))


################
# Run AutoGrow #
################

class run_main():
	"""The class to run the program"""
	
	def __init__(self, args):
		"""Run the program
		
		Arguments:
		args -- an array containing the system parameters. sys.argv format.
		
		"""
		
		# load the commandline parameters
		global vars, printout
		vars,printout = load_in_commandline_parameters(args)
		
		# load the autoclickchem modules
		cmd_folder = os.path.dirname(vars['autoclickchem_script']) 
		if cmd_folder not in sys.path: sys.path.insert(0, cmd_folder)
		global autoclickchem
		import autoclickchem
		
		# open a log file
		global log_text_file
		log_text_file = open(vars['output_dir'] + "log.txt", 'a') 
		log(printout)
		
		# keep track of the calculated smiles strings and druglike status of the ligands so these need only be calculated once
		global smiles, druglike, smiles_used_in_previous_generations, top_compounds
		smiles = {}
		druglike = {}
		smiles_used_in_previous_generations = []
		top_compounds = []
		
		# run the algorithm for the correct number of generations
		generation_num = 0
		
		# keep track of the time each generation takes
		global start_time_current_generation
		start_time_current_generation = 0
		last_generation_since_time_recorded = 0
		
		# turn the receptor.pdb into a receptor.pdbqt
		set_start_time_current_generation(False)
		log("Converting receptor pdb files to pdbqt format if needed...")
		convert_receptor_pdb_files_to_pdbqt()
		
		while generation_num < vars['num_generations']:
			
			# first, identify the directory of the current generation, set up generation-sepcific parameters
			
			current_generation_dir = get_current_generation_directory()
			
			if "generation1" + os.sep in current_generation_dir: vars['number_of_crossovers_current_generation'] = vars['number_of_crossovers_first_generation']
			else: vars['number_of_crossovers_current_generation'] = vars['number_of_crossovers']
			
			if "generation1" + os.sep in current_generation_dir: vars['number_of_mutants_current_generation'] = vars['number_of_mutants_first_generation']
			else: vars['number_of_mutants_current_generation'] = vars['number_of_mutants']
			
			thedirinfo = "\n\nThe current generational directory is " + current_generation_dir
			log(thedirinfo)
			underline = ""
			for t in range(len(thedirinfo)): underline = underline + "="
			log(underline)
			log("")
			
			if current_generation_dir != last_generation_since_time_recorded:
				last_generation_since_time_recorded = current_generation_dir
				set_start_time_current_generation()
		
			check_time_limit_reached()
		
			# now, remove bad ligands
			delete_bad_ligands(current_generation_dir)
			check_time_limit_reached()
			
			# make all the required ligands
			while len(glob.glob(current_generation_dir + "*.pdb")) < vars['number_of_mutants_current_generation'] + vars['number_of_crossovers_current_generation']:
			
				# now see if the latest generation is finished
				make_mutants(current_generation_dir)
				check_time_limit_reached()
		
				# second check if all the crossovers have been made
				make_crossovers(current_generation_dir)
				check_time_limit_reached()
		
				# remove bad ligands
				delete_bad_ligands(current_generation_dir)
				check_time_limit_reached()
			
			# convert ligands to pdbqt format
			log("\nConverting ligand PDB files to PDBQT format...")
			convert_ligand_pdb_files_to_pdbqt(current_generation_dir)
			check_time_limit_reached()
			
			# delete ones that were not converted to PDBQT sucessfully
			delete_no_pdbqt(current_generation_dir)
			check_time_limit_reached()
			
			log("Docking compounds using AutoDock Vina...")
			dock_compounds(current_generation_dir)
			check_time_limit_reached()
			
			# delete ones that were not docked successfully
			delete_no_docked(current_generation_dir)
			check_time_limit_reached()
			
			# if the proper number of ligands have been generated and docked, advance to the next generation. Otherwise, try this generation all over again.
			total_num_compoounds = vars['number_of_mutants_current_generation'] + vars['number_of_crossovers_current_generation']
		
			if len(glob.glob(current_generation_dir + "*.pdb")) == total_num_compoounds and len(glob.glob(current_generation_dir + "*.pdbqt")) == total_num_compoounds and len(glob.glob(current_generation_dir + "*.vina")) == total_num_compoounds:
			
				log("Identifying top compounds from recent dockings...")
				top_compounds = identify_top_compounds(current_generation_dir)
		
				make_next_generation_directory(current_generation_dir)
				check_time_limit_reached()
		
				vars['directory_of_source_compounds'] = current_generation_dir + "best_ligands" + os.sep
				
				generation_num = generation_num + 1
		
		log_text_file.close()

if __name__=="__main__": dorun = run_main(sys.argv)

