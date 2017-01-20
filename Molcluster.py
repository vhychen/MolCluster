#!/usr/bin/python
import copy
import os.path as op
import os
import MolClass
import re
import shutil
from decimal import Decimal
from GlobalVariables import *
import sys
import ChemShellInputFunctions as CSIF
import datetime
import textwrap
import math
import BondCorrelationFunctions as BCF
import CommonFunctions as CF
import TrajectoryGeometryAnalysis
import GenerateEmbeddedClusterInput as GECI
import ReadInput as RI

#-------------------------------------------------------------------------------
# Function main
#
# Operation: Main user interface of MolCluster
#-------------------------------------------------------------------------------
def main():

  temp_dir = 'temp'
  RI.genTempDirs()

  log_file = op.join(temp_dir,'molcluster.log')
  
  version='2.40'
  
  ######################################
  ### Placeholder for cluster radius ###
  ######################################
  radius = None
  
  ######################################
  ### Some common outputs to screen. ###
  ######################################
  bugs_line = 'Please report any bugs to Vincent Chen (vc0489@gmail.com). The log file can be found in \'' + str(log_file) + '\'.'
  exit_line = '\n' + 'E'*80 + '\nEEE  ' + 'Exiting molcluster'.center(70) + '  EEE\n' + 'E'*80 + '\n'
  opening_text = double_dash_line + '\n===' + ('MolCluster V' + version + ', a molecular cluster cutting tool.').center(74) + '===\n===' + 'Vincent Chen (vc0489@gmail.com)'.center(74) + '===\n' + double_dash_line
  print (opening_text)
  
  ########################
  ### Create log file. ###
  ########################
  with open(log_file,'w') as f:
    f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + '\n')
    f.write(opening_text + '\n') 
  
  ###############################
  ### Begin user interface... ###
  ###############################
  user_input = None
  while not user_input:
    user_input = raw_input('Please give directory containing the FIELD and HISTORY files: ').strip()
    user_input = user_input.strip()
    field_file = op.join(user_input,'FIELD')
    history_file = op.join(user_input,'HISTORY')
    if not op.isdir(user_input):
      print('The directory \'' + user_input + '\' does not exist.')
      print (dash_line)
      with open(log_file,'a') as f:
        f.write('The directory \'' + user_input + '\' does not exist.\n')
        f.write(dash_line + '\n')
      user_input = None
    elif not op.exists(field_file):
      print ('\'FIELD\' is not present in \'' + user_input + '\'.')
      print (dash_line)
      with open(log_file,'a') as f:
        f.write('\'FIELD\' is not present in \'' + user_input + '\'.\n')
        f.write(dash_line + '\n')
      user_input = None
    elif not op.exists(history_file):
      print ('\'HISTORY\' is not present in \'' + user_input + '\'.')
      print (dash_line)
      with open(log_file,'a') as f:
        f.write('\'HISTORY\' is not present in \'' + user_input + '\'.\n')
        f.write(dash_line + '\n')
      user_input = None
    else: 
      input_dir = user_input
      with open(log_file,'a') as f:
        f.write('FIELD and HISTORY files found in \'' + input_dir + '\'.\n')
  print (dash_line)
  print ('FIELD and HISTORY files found in \'' + input_dir + '\'. Reading files...')
  mol_info = None
  try:
    mol_info,vdw_info,n_atoms = RI.field(field_file)
  except:
    #mol_info,vdw_info = ReadInput.field(field_file)
    print ('FIELD file format not as expected. ' + bugs_line)
    print (exit_line)
    with open(log_file,'a') as f:
      f.write('FIELD file format not as expected.\n')
      f.write(exit_line)
    return
  print ('')
  #n_atoms = 0
  n_mols = 0
  for t in mol_info:
    #n_atoms += t['nmols']*t['natoms']
    n_mols += t['nmols']
  #n_configs,timestep_size = ReadInput.configurationsInHistory(history_file,n_atoms)
  try:
    n_configs,timestep_size,ns_between_configs = RI.configurationsInHistory(history_file,n_atoms)
    if n_configs == -1:
      print ('Mismatch in the number of atoms expected (from FIELD) and number of atoms in HISTORY. ' + bugs_line)
      print (exit_line)
      with open(log_file,'a') as f:
        f.write('Mismatch in the number of atoms expected (from FIELD) and number of atoms in HISTORY.\n')
        f.write(exit_line)
      return
    elif n_configs == 0:
      print ('No configurations detected in HISTORY. ' + bugs_line)
      print (exit_line)
      with open(log_file,'a') as f:
        f.write('No configurations detected in HISTORY.\n')
        f.write(exit_line)
      return
  except:
    print ('HISTORY file format not as expected. ' + bugs_line)
    print (exit_line)
    with open(log_file,'a') as f:
      f.write('HISTORY file format not as expected.\n')
      f.write(exit_line)
    return

  with open(log_file,'a') as f:
    f.write('A total of ' + str(n_configs) + ' configuration(s) found in HISTORY.\n')

  #first_config,unit_cell_length,first_timestep = ReadInput.readHistory(history_file,mol_info,None,0)
  config_list,timestep_list,unit_cell_length_list = RI.genCoordinatesListFromHistory(history_file,mol_info,[1],n_atoms,timestep_size)
  first_system = MolClass.BulkSystem(config_list[0],mol_info,unit_cell_length_list[0],vdw_info)
  first_system.print_system_info()
  
  bulk_system_info_file = op.join(temp_dir,'bulk_system_info.txt')
  first_system.print_system_info('BULK SYSTEM INFORMATION',bulk_system_info_file,False)
  
  user_input = None
  while not user_input:
    print('Would you like to:') 
    print('-->[1] Cut bare cluster(s) for (ChemShell) QM/MM simulation(s)')
    print('-->[2a] Generate ChemShell input for cutting embedded clusters')
    print('-->[2b] Cut cluster(s) and add point charges (generated after running [2a])')
    print('-->[3] Geometry analysis of a DL_POLY trajectory')
    print('-->[4] Calculate the bond correlation function between two specified atoms')
    user_input = raw_input('Please choose option [1/2a/2b/3/4]: ').strip()
    if user_input == '4':
      with open(log_file,'a') as f:
        f.write('User chose to calculate the bond correlation function between two specified atoms.\n')
      BCF.userInterface(history_file,mol_info,n_atoms,first_system,n_configs,timestep_size,ns_between_configs,input_dir)
      return
    elif user_input == '1':
      with open(log_file,'a') as f:
        f.write('User chose to cut bare cluster(s) for (ChemShell) QM/MM simulation(s).\n')
      print (dash_line)
    elif user_input == '3':
      with open(log_file,'a') as f:
        f.write('User chose to run geometry analysis on a DL_POLY trajectory.\n')
      print (dash_line)
      TrajectoryGeometryAnalysis.main(input_dir)
      return
    elif user_input == '2a':
      with open(log_file,'a') as f:
        f.write('User chose to generate ChemShell input for cutting embedded cluster(s).\n')
      print (dash_line)
      GECI.main(input_dir,history_file,first_system,log_file,n_configs,mol_info,timestep_size,vdw_info)
      return
    elif user_input == '2b':
      with open(log_file,'a') as f:
        f.write('User chose to cut cluster(s) and merge with ChemShell embedded cluster(s) fragment.\n')
      print (dash_line)
      user_input_2 = None
      while user_input_2 == None:
        user_input_2 =  raw_input(textwrap.fill('Please enter parent directory (containing \'embedded_cluster_specs.obj\') of the embedded cluster(s) generated by ChemShell:',80) + ' ').strip()  
        if not op.isdir(user_input_2):
          print ('The directory ' + user_input_2 + ' does not exist.')
          print (dash_line)
          user_input_2 = None
        else:
          specs_file = op.join(user_input_2,'embedded_cluster_specs.obj')
          if not op.exists(specs_file):
            print (textwrap.fill('The file embedded_cluster_specs.obj cannot be found. Please copy it into ' + user_input_2,80))
          embedded_cluster_specs_obj = RI.loadSystemObj(specs_file) 
          try:
            embedded_cluster_specs_obj = RI.loadSystemObj(specs_file) 
          except:
            print ('Error in reading embedded_cluster_specs.obj. Exiting.')
            return
          print ('Expecting ' + str(embedded_cluster_specs_obj.n_clusters) + ' embedded clusters...')
          for i in range(embedded_cluster_specs_obj.n_clusters):
            cluster_dir = op.join(user_input_2,embedded_cluster_specs_obj.embedded_cluster_prefix + '_' + str(i+1))
            if not op.isdir(cluster_dir):
              print ('The directory ' +  cluster_dir + ' does not exist. Exiting.')
              return 
            fragment_file = op.join(cluster_dir,'embedded_cluster.pun')
            if not op.exists(fragment_file):
              print ('The fragment file ' + fragment_file + ' does not exist. Exiting.')
              return
      GECI.genFullEmbeddedClusters(user_input_2,embedded_cluster_specs_obj,log_file)
      return
    else:
      user_input = None
      print ('Please enter a valid option.')
      print (dash_line)

  #############################################
  # Choose configurations to extract clusters #
  #############################################
  prompt_message = 'Would you like to extract a cluster from:'
  config_no_list = CF.chooseConfigurations(n_configs,prompt_message,log_file)
  print ('config_no_list: ' + str(config_no_list))
  print (dash_line)

  try: 
    configuration_list,timestep_list,unit_cell_length_list = RI.genCoordinatesListFromHistory(history_file,mol_info,config_no_list,n_atoms,timestep_size)
  except:
    print (textwrap.fill('Error when extracting configuration(s) from HISTORY. ' + bugs_line,80))
    print (exit_line)
    with open(log_file,'a') as f:
      f.write('Error when extracting configuration(s) from HISTORY.\n')
      f.write(exit_line)
    return
  #print ('...complete.')
  print (dash_line)
  
  ############################
  # Build BulkSystem objects #
  ############################
  print ('Creating BulkSystem objects from configuration(s)...')
  system_list = []
  try:
    for i,c in enumerate(configuration_list):
      s = MolClass.BulkSystem(c,mol_info,unit_cell_length_list[i],vdw_info,timestep_list[i],(timestep_list[i]-1)*timestep_size)
      system_list.append(s)
      with open(log_file,'a') as f:
        f.write('-->Configuration extracted from timestep ' + str(timestep_list[i]) + ' (time = ' + str(s.time/1000) + ' ns).\n')
  except:
    print (textwrap.fill('Error when building BulkSystem object(s). ' + bugs_line,80))
    print (exit_line)
    with open(log_file,'a') as f:
      f.write('Error when building BulkSystem object(s).\n')
      f.write(exit_line)
    return
  print ('...complete.')
  print (dash_line)
  with open(log_file,'a') as f:
    f.write('BulkSystem object created for all configuration(s) extracted.\n')

  #########################
  # Define cluster origin #
  #########################
  origin_type,origin_type_var = CF.chooseOrigin(None,log_file,first_system)
  print (dash_line)
  
  ###################################################
  # User chooses the way the cluster is constructed #
  ###################################################
  user_input = None
  while not user_input:
    print ('Available cluster construction definition:')
    print ('-->[R] Fixed radius')
    #print ('-->[T] Specify total number of molecules - NOT SUPPORTED YET')
    print ('-->[M] Specify number of each molecule type') 
    user_input = raw_input('Please choose cluster size definition [R/M] --> ').strip()
    upper_user_input = user_input.upper()
    if upper_user_input == 'R':
      # User defines cluster size (radius/number of atoms)
      cluster_construction_type = 'R'
      second_user_input = None
      first_system.print_estimate_n_atoms_in_cluster()
      while not second_user_input:
        second_user_input = raw_input('Enter the radius (in Angstroms) of the cluster(s) to be cut: ').strip()
        try:
          radius = float(second_user_input)
          if radius < 0:
            second_user_input = None
            print ('Please enter a positive number.')
            print (dash_line)
          else:
            pass
        except:
          print (str(second_user_input) + ' is not a valid radius length.')
          print (dash_line)
          second_user_input = None
      with open(log_file,'a') as f:
        f.write('Cluster(s) with a radius of ' + str(radius) + ' ' + ANGSTROM_STR + ' will be cut.\n')
      
      ########################
      # Define boundary type #
      ########################
      print (dash_line)
      secound_user_input = None
      while not secound_user_input:
        print ('The boundary can be:')
        print (textwrap.fill('-->[H]ard: all atoms have to within the cluster boundary defined by the radius.',80))
        print (textwrap.fill('-->[S]oft: atoms can be outside of the cluster boundary as long as the centre of mass of its parent molecule is within the cluster boundary',80))
        secound_user_input = raw_input('Choose type of boundary [H/S]: ').strip().upper()
        if secound_user_input == 'H':
          radius_type = 'H'
          radius_type_text = 'Hard'
        elif secound_user_input == 'S':
          radius_type = 'S'
          radius_type_text = 'Soft'
        else:
          secound_user_input = None
          print ('Please enter \'H\' or \'S\'')
          print (dash_line)

      print (radius_type_text + ' boundary selected.')
      with open(log_file,'a') as f:
        f.write(radius_type_text + ' boundary type selected.\n')

      print (dash_line)
      secound_user_input = None
      while not secound_user_input:
        secound_user_input = raw_input('Are the cluster(s) required to be neutral [Y/N]? ').strip().upper()
        if secound_user_input == 'Y':
          neutral = True
          with open(log_file,'a') as f:
            f.write('User requested cluster(s) to be neutral.\n')
        elif secound_user_input == 'N':
          neutral = False
        else:
          secound_user_input = None
          print ('Please enter \'Y\' or \'N\'')
          print (dash_line)
      #print (dash_line)
    elif upper_user_input == 'T':
      #cluster_construction_type = 'T'
      print (user_input + ' is not a valid option')
      user_input = None
      print (dash_line)
    elif upper_user_input == 'M':
      cluster_mol_numbers_list,est_cluster_radius = CF.define_region_of_number_of_each_molecule_type('cluster',first_system,'BULK SYSTEM INFO',True,log_file)
      cluster_construction_type = 'M'
      with open(log_file,'a') as f:
        f.write('Cluster(s) with specified numbers of each molecule type to be cut.\n')
    else:
      print (user_input + ' is not a valid option')
      user_input = None
      print (dash_line)
  
  #################
  # Choose output #
  #################
  print (dash_line)
  user_input = None
  while not user_input:
    print ('The following output formats are available:')
    print ('-->[C]hemShell fragment input file (with accompanying .xyz file)')
    print ('-->[X]yz co-ordinates file')
    user_input = raw_input('Please choose the output format [C/X]: ').strip().upper()
    if user_input in ('C','X'):
      output_format = user_input
    else:
      print ('Please enter \'C\' or \'X\'')
      print (dash_line)
      user_input = None

  ####################
  # ChemShell inputs #
  ####################
  print (dash_line) 
  if output_format == 'C':
    user_input = None
    with open(log_file,'a') as f:
      f.write('User requested ChemShell fragment input files to be generated.\n')
    while not user_input:
      user_input = raw_input(textwrap.fill('Would you also like to generate Chemshell input files for QM/MM optimisations [Y/N]? ',80) + ' ' ).strip()
      if user_input in ('Y','y'):
        with open(log_file,'a') as f:
          f.write('User also requested ChemShell input files for QM/MM optimisations to be generated.\n')
        gen_chemshell_opt_input_flag = True

        ########################
        # Define the QM region #
        ########################
        if radius:
          qm_region_type,qm_region_var = CF.chooseQMRegionDef(radius,first_system,log_file)
        else:
          qm_region_type,qm_region_var = CF.chooseQMRegionDef(est_cluster_radius,first_system,log_file)
        if qm_region_type == 'S':
          qm_radius = qm_region_var
        elif qm_region_type == 'M':
          n_qm_mol = qm_region_var
        elif qm_region_type == 'A':
          n_qm_atoms = qm_region_var
        
        ############################
        # Define the active region #
        ############################
        if radius:
          active_region_type,active_region_var = CF.chooseActiveRegionDef(radius,first_system,log_file)
        else:
          active_region_type,active_region_var = CF.chooseActiveRegionDef(est_cluster_radius,first_system,log_file)
        active_radius = None
        if active_region_type == 'S':
          active_radius = active_region_var
        elif active_region_type == 'M':
          n_active_mol = active_region_var
        elif active_region_type == 'A':
          n_active_atoms = active_region_var

        ##############################
        # Choose co-ordinates system #
        ##############################
        print (dash_line) 
        user_input2 = None
        while not user_input2:
          print ('Available co-ordinate systems:')
          print ('[C]artesian (constraints not permitted)')
          print ('[D]elocalised internal coordinates (DLC)')
          print ('[H]ybrid delocalised internal coordinates (HDLC)')
          user_input2 = raw_input(textwrap.fill('Please select co-ordinate system to be used during optimisation [C/D/H]:',80) + ' ').strip().upper()
          if user_input2 == 'C':
            coord_system = 'C'
          elif user_input2 == 'D':
            coord_system = 'D'
          elif user_input2 == 'H':
            coord_system = 'H'
          else:
            user_input2 = None
            print ('Please enter \'C\', \'D\' or \'H\'.')
            print (dash_line)

      elif user_input in ('N','n'):
        print ('ChemShell input files for QM/MM optimisations will not be generated.')
        gen_chemshell_opt_input_flag = False
        with open(log_file,'a') as f:
          f.write('ChemShell input files for QM/MM optimisations will not be generated.\n')
      else:
        print (user_input + ' is not a valid option.')
        user_input = None
        print (dash_line)
  elif output_format == 'X':
    with open(log_file,'a') as f:
      f.write('User requested .xyz co-ordinates files to be output.\n')
  print (dash_line) 

  user_input = None
  while not user_input:
    user_input = raw_input('Please give output directory: ').strip()
    with open(log_file,'a') as f:
      f.write('The output directory is \'' + user_input + '\'.\n')
    if not op.isdir(user_input):
      #print('The directory \'' + user_input + '\' does not exist.')
      print('Creating directory \'' + user_input + '\'')
      # VC: restrict alphanumeric and underscore characters?
      output_dir = user_input
      os.makedirs(output_dir)
      with open(log_file,'a') as f:
        f.write('Output directory \'' + user_input + '\' created.\n')
    else:
      output_dir = user_input
  
  print (dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('Please give output basename (no spaces and only alphanumeric and underscore characters allowed):',80) + ' ').strip()
    if len(user_input.split()) > 1:
      print ('Please enter basename without any spaces')
      user_input = None
      print (dash_line)
    else:
      if not basename_p.match(user_input):
        print ('Please only use alphanumeric and underscore characters')
        user_input = None
        print (dash_line)
      else:
        output_basename = user_input
  with open(log_file,'a') as f:
    f.write('The basename of the folders for each cluster is \'' + output_basename + '\'.\n')
  shutil.copy(bulk_system_info_file,output_dir)
  with open(log_file,'a') as f:
    f.write('Info file of original bulk system is copied to \'' + output_dir + '\'.\n')
  os.remove(bulk_system_info_file)
   
  print (dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('Do you want the bulk and cluster data objects to be generated for debug purposes [Y/N]? ',80) + ' ').strip().upper()
    if user_input == 'Y':
      output_data_obj = True
      with open(log_file,'a') as f:
        f.write('Bulk and cluster data objects will be generated for debug purposes.\n')
    elif user_input == 'N':
      output_data_obj = False
    else:
      user_input = None
      print ('Please enter \'Y\' or \'N\'')
      print (dash_line)
  
  #####################
  # Generate clusters #
  #####################
  for n,s in enumerate(system_list):
    cluster_dir = op.join(output_dir,'cluster_' + str(n+1))
    cluster_name = 'CLUSTER ' + str(n+1)
    with open(log_file,'a') as f:
      f.write(dash_line + '\n')
      f.write('Cluster ' + str(n+1) + ' of ' + str(len(system_list)) + ':\n')
    if not op.isdir(cluster_dir):
      print (dash_line)
      print('Creating directory \'' + cluster_dir + '\'')
      print (dash_line)
      os.makedirs(cluster_dir)
      with open(log_file,'a') as f:
        f.write('-->Directory \'' + cluster_dir + '\' created.\n')
    if output_data_obj == True:
      RI.saveSystemObj(s,op.join(cluster_dir,'bulk.system'))
      with open(log_file,'a') as f:
        f.write('Bulk data object written to ' + op.join(cluster_dir,'bulk.system') + '.\n')
    if cluster_construction_type == 'R':
      cluster_obj = s.cut_cluster_of_defined_radius(radius,origin_type_var,origin_type,radius_type,neutral)
    #elif cluster_construction_type == 'T':
    #  pass
    elif cluster_construction_type == 'M':
      cluster_obj = s.cut_cluster_of_specified_n_of_each_molecule(cluster_mol_numbers_list,origin_type_var,origin_type)
    cluster_obj.name = cluster_name
    print ('Cluster ' + str(n+1) + ' out of ' + str(len(system_list)) + ' generated...') 
    #with open(log_file,'a') as f:
    #  f.write('-->Cluster ' + str(n+1) + ' out of ' + str(len(system_list)) + ' generated.\n')
    if output_format == 'X':
      output_file = op.join(cluster_dir,output_basename + '_' + str(n+1) + '.xyz')
      cluster_obj.write_xyz(output_file,'Cluster generated by MolCluster')
      with open(log_file,'a') as f:
        f.write('-->.xyz co-ordinates file \'' + output_file + '\' generated.\n')
    elif output_format == 'C':
      output_file = op.join(cluster_dir,output_basename + '_' + str(n+1) + '.chm')
      cluster_obj.write_chemshell_fragment(output_file,'cluster.pun')
      with open(log_file,'a') as f:
        f.write('-->ChemShell fragment input file \'' + output_file + '\' generated.\n')
      xyz_output_file = op.join(cluster_dir,output_basename + '_' + str(n+1) + '.xyz')
      with open(log_file,'a') as f:
        f.write('-->.xyz co-ordinates file \'' + xyz_output_file + '\' generated.\n')
      cluster_obj.write_xyz(xyz_output_file,'Cluster generated by MolCluster')
      ff_file = op.join(cluster_dir,'ff.dat')
      print ('ChemShell potential file ' +  ff_file + '...')
      cluster_obj.gen_chemshell_ff_file(ff_file)
      with open(log_file,'a') as f:
        f.write('-->ChemShell force field file \'' + ff_file + '\' generated.\n')
      print ('...generated')
      print (dash_line)

      ####################################################
      # Generate ChemShell QM/MM optimisation input file #
      ####################################################
      if gen_chemshell_opt_input_flag == True:

        ######################
        # Generate QM region #
        ######################
        if qm_region_type == 'S':
          cluster_obj.det_atom_list_for_given_radius(qm_radius,'QM')
          print ('There are ' + str(len(cluster_obj.qm_mol_list)) + ' molecules (' + str(len(cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
        elif qm_region_type == 'M':
          if n_qm_mol > len(cluster_obj.mol_list):
            print (textwrap.fill('Requested number of molecules to be included in QM region is greater than the total number of molecules in the cluster. All molecules will be in the QM region.',80))
            cluster_obj.det_atom_list_for_given_n_mol(len(cluster_obj.mol_list),'QM') 
          else:
            cluster_obj.det_atom_list_for_given_n_mol(n_qm_mol,'QM') 
          print ('There are ' + str(len(cluster_obj.qm_mol_list)) + ' molecules (' + str(len(cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
        elif qm_region_type == 'A':
          #qm_atom_list,qm_mol_list = cluster_obj.det_atom_list_for_given_n_atoms(n_qm_atoms) 
          if n_qm_atoms > len(cluster_obj.atom_list):
            print (textwrap.fill('Requested number of atoms to be included in QM region is greater than the total number of atoms in the cluster. All molecules will be in the QM region.',80))
            cluster_obj.det_atom_list_for_given_n_atoms(len(cluster_obj.atom_list),'QM')
          else:
            cluster_obj.det_atom_list_for_given_n_atoms(n_qm_atoms,'QM')
          print ('There are ' + str(len(cluster_obj.qm_mol_list)) + ' molecules (' + str(len(cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
        elif qm_region_type == 'N':
          if n == 0: 
            qm_mol_numbers_list = CF.define_region_of_number_of_each_molecule_type('QM region',cluster_obj,'CLUSTER SYSTEM INFO',False,log_file)
          cluster_obj.closest_n_molecules_of_each_type(qm_mol_numbers_list,True,True,'QM')
        elif qm_region_type == 'C':
          if active_radius == 'A':
            cluster_obj.choose_qm_molecules(cluster_name,None,'QM')
          elif active_radius == None:
            cluster_obj.choose_qm_molecules(cluster_name,None,'QM')
          else:
            cluster_obj.choose_qm_molecules(cluster_name,active_radius,'QM')
        elif qm_region_type == 'O': 
          qm_regions = CF.choose_qm_regions(cluster_obj,log_file)
          #print (qm_regions)
          cluster_obj.det_atom_list_for_set_of_qm_regions(qm_regions,'QM')
          print ('There are ' + str(len(cluster_obj.qm_mol_list)) + ' molecules (' + str(len(cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
        print ('The total charge of the QM region is ' + str(cluster_obj.qm_charge) + '.')
        print (dash_line)
        with open(log_file,'a') as f:
          f.write('-->There are ' + str(len(cluster_obj.qm_mol_list)) + ' molecules (' + str(len(cluster_obj.qm_atom_list)) + ' atoms) in the QM region:\n')
          for m in cluster_obj.qm_mol_list:
            f.write('---->Molecule ' + str(m) + ' (type: ' + cluster_obj.mol_list[m-1].label + '; charge: ' + str(cluster_obj.mol_list[m-1].charge) + '; dist. to origin: ' + str(cluster_obj.mol_list[m-1].dist_to_origin) + ' ' + ANGSTROM_STR + ')\n')
          f.write('-->The total charge of the QM region is ' + str(cluster_obj.qm_charge) + '.\n')

        ##########################
        # Generate active region #
        ##########################
        if active_region_type == 'S':
          if active_radius == 'A':
            pass
            #cluster_obj.gen_tiered_label('1','2','3',True)
          else:
            cluster_obj.det_atom_list_for_given_radius(active_radius,'Active')
        elif active_region_type == 'M':
          if n_active_mol > len(cluster_obj.mol_list):
            print (textwrap.fill('Requested number of molecules to be included in active region is greater than the total number of molecules in the cluster. All molecules will be in the QM region.',80))
            cluster_obj.det_atom_list_for_given_n_mol(len(cluster_obj.mol_list),'Active') 
          else:
            active_list,active_mol_list = cluster_obj.det_atom_list_for_given_n_mol(n_active_mol) 
          print ('There are ' + str(len(cluster_obj.active_mol_list)) + ' molecules (' + str(len(cluster_obj.active_atom_list)) + ' atoms) in the active region.')
        elif active_region_type == 'A':
          if n_qm_atoms > len(cluster_obj.atom_list):
            print (textwrap.fill('Requested number of atoms to be included in active region is greater than the total number of atoms in the cluster. All molecules will be in the QM region.',80))
            cluster_obj.det_atom_list_for_given_n_atoms(len(cluster_obj.atom_list),'Active')
          else:
            cluster_obj.det_atom_list_for_given_n_atoms(n_active_atoms,'Active')
          print ('There are ' + str(len(cluster_obj.active_mol_list)) + ' molecules (' + str(len(cluster_obj.active_atom_list)) + ' atoms) in the active region.')
        elif active_region_type == 'N':
          if n == 0: 
            active_mol_numbers_list = CF.define_region_of_number_of_each_molecule_type('active region',cluster_obj,'CLUSTER SYSTEM INFO',False,log_file)
          cluster_obj.closest_n_molecules_of_each_type(active_mol_numbers_list,True,True,'Active')
        with open(log_file,'a') as f:
          if len(cluster_obj.active_mol_list) > 0:
            f.write('-->There are ' + str(len(cluster_obj.active_mol_list)) + ' molecules (' + str(len(cluster_obj.active_atom_list)) + ' atoms) in the active region.\n')
          else:
            f.write('-->All molecules are in the active region.\n')
        
        #####################################
        # Generate tiered cluster .xyz file #
        #####################################
        cluster_obj.gen_tiered_label('1','2','3')
        tiered_cluster_file = op.join(cluster_dir,output_basename + '_' + str(n+1) + '_tiered.xyz')
        cluster_obj.write_xyz(tiered_cluster_file,'Cluster generated by MolCluster',True)
        with open(log_file,'a') as f:
          f.write('-->Tiered cluster .xyz file written to ' + tiered_cluster_file + '.\n')

        ####################################################
        # Generate ChemShell QM/MM optimisation input file #
        ####################################################
        CSIF.qmmm_opt_input(cluster_obj,op.join(cluster_dir,'opt.chm'),'cluster.pun','cluster_opt.pun',cluster_obj.qm_charge,1,True,log_file,coord_system)
        with open(log_file,'a') as f:
          f.write('-->ChemShell QM/MM optimisation input written to \'' + op.join(cluster_dir,'opt.chm') + '\'.\n')
        if output_data_obj == True:
          RI.saveSystemObj(cluster_obj,op.join(cluster_dir,'cluster.system'))
          with open(log_file,'a') as f:
            f.write('-->Cluster data object file written to \'' + op.join(cluster_dir,'cluster.system') + '\'.\n')
    
    system_info_output = op.join(cluster_dir,'system_info.txt')
    cluster_obj.print_system_info('CLUSTER ' + str(n+1) + ' INFORMATION',system_info_output)
    cluster_obj.update_atom_and_mol_region_attribute()
    cluster_obj.gen_atom_by_atom_info(op.join(cluster_dir,'atom_by_atom.info'))
    cluster_obj.gen_mol_by_mol_info(op.join(cluster_dir,'mol_by_mol.info'))
    with open(log_file,'a') as f:
      f.write('Cluster system info written to ' + system_info_output + '.\n')
      f.write('Atom by atom info file written to ' + op.join(cluster_dir,'atom_by_atom.info') + '.\n')
      f.write('Molecule by molecule info file written to ' + op.join(cluster_dir,'mol_by_mol.info') + '.\n')

  print (dash_line)
  print ('')
  print ('I'*80)
  print ('III' + ('  All clusters generated.  ').center(74) + 'III') 
  print ('I'*80)
  print ('')
  print (dash_line) 
  print (textwrap.fill(bugs_line,80))
  with open(log_file,'a') as f:
    f.write(dash_line + '\n')
    f.write('All clusters generated.\n')
    f.write(exit_line)
  shutil.copy(log_file,output_dir)
  print (dash_line)
  print ('O'*80)
  log_info = textwrap.wrap('Log file copied to \'' + output_dir + '\'.',70)
  for l in log_info:
    print ('OOO  ' + l.center(70) + '  OOO')
  #print ('OOO' + ('Log file copied to \'' + output_dir + '\'.').center(74) + 'OOO')
  print ('O'*80)
  print (dash_line)
  print (exit_line)
  return

if __name__ == "__main__":
  sys.exit(main())

#####################
# DEFUNCT FUNCTIONS #
#####################

#-------------------------------------------------------------------------------
# Interface to choose a custom set of one or more QM regions.
#-------------------------------------------------------------------------------
#def choose_qm_regions(cluster_obj,log_file=None):
#  n_regions = None
#  while not n_regions:
#    n_regions = raw_input('Enter the number of QM regions you would like to define. --> ').strip()
#    if n_regions.isdigit() == False:
#      print (n_regions + ' is not a number')
#      n_regions = None
#      print (dash_line)
#    elif n_regions < 1:
#      print ('Please enter an integer greater than 0')
#      n_regions = None
#      print (dash_line)
#    else:
#      n_regions = int(n_regions)
#  if log_file:
#    with open(log_file,'a') as f: 
#      f.write('The QM region will consist of ' + str(n_regions) + ' spherical regions in total.\n')
#  print (dash_line)
#  regions_def = []
#  while len(regions_def) < n_regions:
#    origin_type = None
#    while not origin_type:
#      print ('QM region no. ' + str(len(regions_def)+1) + '. Possible origin types:')
#      print ('-->[A]tom')
#      print ('-->[M]olecule')
#      print ('-->[C]luster origin')
#      origin_type = raw_input('Please choose origin type [A/M/C] --> ').strip()
#      if origin_type == 'A':
#        #ordered_atom_list = cluster_obj.ret_ordered_atom_list()
#        origin_atom_no = cluster_obj.choose_atom_origin(len(regions_def)+1)
#        origin_coords = cluster_obj.atom_list[origin_atom_no-1].coords
#      elif origin_type == 'M':
#        cluster_obj.print_mol_list(cluster_obj.name)
#        origin_mol_no = None
#        while not origin_mol_no:
#          origin_mol_no = raw_input(textwrap.fill('Please choose the molecule to be used as the origin of QM region no. ' + str(len(regions_def)+1) + '. Note that the centre of mass of the molecule will be used. --> ',80) + ' ').strip()
#          if not origin_mol_no.isdigit():
#            print (origin_mol_no + ' is not valid. Please enter an integer between 1 and ' + str(len(cluster_obj.mol_list)) + '.')
#            print (dash_line)
#            origin_mol_no = None
#          elif int(origin_mol_no) < 1 or int(origin_mol_no) > len(cluster_obj.mol_list):
#            print (origin_mol_no + ' is not valid. Please enter an integer between 1 and ' + str(len(cluster_obj.mol_list)) + '.')
#            print (dash_line)
#            origin_mol_no = None
#          else:
#            origin_mol_no = int(origin_mol_no)
#            origin_coords = cluster_obj.mol_list[origin_mol_no-1].centre_of_mass
#      elif origin_type == 'C':
#        origin_coords = (0.0,0.0,0.0)
#      else:
#        print (origin_type + ' is not a valid option')
#        print (dash_line)
#        origin_type = None
#        continue 
#      qm_radius = CF.choose_qm_radius(cluster_obj.radius,cluster_obj)
#      if log_file:
#        with open(log_file,'a') as f:
#          if origin_type == 'A':
#            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about atom ' + str(origin_atom_no) + ' [type: ' + cluster_obj.atom_list[origin_atom_no-1].label + '; coords: ' + coords_to_string(origin_coords) + '] with a radius of ' + str(qm_radius) + ' ' + u'\u212B'.encode('utf-8') + '.\n')
#          elif origin_type =='M':
#            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about the centre of mass of molecule ' + str(origin_mol_no) + ' [type: ' + cluster_obj.mol_list[origin_mol_no-1].label + '; coords: ' + coords_to_string(origin_coords) + '] with a radius of ' + str(qm_radius) + ' ' + u'\u212B'.encode('utf-8') + '.\n')
#          elif origin_type == 'C':
#            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about the cluster origin with a radius of ' + str(qm_radius) + u'\u212B'.encode('utf-8') + '.\n')
#    regions_def.append((origin_coords,qm_radius))
#  return (regions_def)


#def coords_to_string(coords):
#  string = '('
#  for i in coords:
#    string += str(i) + ','
#  string = string[:-1] + ')'
#  return (string)


#-------------------------------------------------------------------------------
# 1. Generate bulk system fragment inputs - Done
# 2. Embedded cluster inputs - TODO
# 3. Generate submit scripts - TODO
# NOTE: This function is obsolete now.. remove soon
#-------------------------------------------------------------------------------
#def user_interface_for_chemshell(system_list,mol_info,vdw_info):
#  user_input = None
#  while not user_input:
#    user_input = raw_input('Please enter base directory for output: ').strip()
#    if len(user_input.split()) > 1:
#      print ('Please enter directory name without any spaces')
#      user_input = None
#      print (dash_line)
#    elif not op.isdir(user_input):
#      print ('The directory ' + user_input + ' does not exist, creating directory.')
#      print (dash_line)
#      base_dir = user_input
#      os.makedirs(base_dir)
#    else:
#      base_dir = user_input
#  print (dash_line)
#  ff_file = base_dir + os.sep + 'ff.dat'
#  print ('ChemShell potential file ' +  ff_file + '...')
#  ReadInput.gen_chemshell_ff_file(mol_info,vdw_info,ff_file)
#  user_input = None
#  while not user_input:
#    user_input = raw_input('Please enter directory prefix for each fragment input (no spaces and only alphanumeric and underscore characters allowed): ').strip()
#    if len(user_input.split()) > 1:
#      print ('Please enter directory prefix without any spaces')
#      user_input = None
#      print (dash_line)
#    else:
#      if not basename_p.match(user_input):
#        print ('Please only use alphanumeric and underscore characters')
#        user_input = None
#        print (dash_line)
#      else:
#        output_dir_prefix = user_input
#   
#  for i,s in enumerate(system_list):
#    output_dir = op.join(base_dir,output_dir_prefix + str(i+1))
#    if not op.isdir(output_dir):
#      os.makedirs(output_dir)
#      shutil.copyfile(ff_file,op.join(output_dir,'ff.dat'))
#    fragment_file = op.join(output_dir,'bulk_fragment.chm')
#    print ('Creating fragment input ' + fragment_file + ' ...') 
#    s.write_chemshell_fragment(fragment_file)
#    print ('...complete.')
#  return

#  ######################################################################
#  # Possible merge of embedded cluster point charges with bare cluster #
#  ######################################################################
#  print (dash_line) 
#  user_input = None
#  while not user_input:
#    user_input = raw_input(textwrap.fill('Would you like to include point charges from a ChemShell embedded cluster fragment?')) 
