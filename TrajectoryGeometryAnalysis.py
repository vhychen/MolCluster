#!/usr/bin/python

from decimal import Decimal
import CommonFunctions as CF
import GlobalVariables as GV
import sys
import ReadInput as RI
import os
import MolClass
import textwrap
import copy
import os.path as op

#-------------------------------------------------------------------------------
# Function main
#
# Operation: Main interface to analyse a DL_POLY trajectory
#-------------------------------------------------------------------------------
def main(workdir):
  RI.genTempDirs()
  field_file = workdir + os.sep + 'FIELD'
  history_file = workdir + os.sep + 'HISTORY'

  mol_info,vdw_info,n_atoms = RI.field(field_file)
  n_configs,timestep_size,ns_between_configs = RI.configurationsInHistory(history_file,n_atoms)
  print (str(n_configs) + ' configurations found in HISTORY')
  user_input = None
  while not user_input:
    print ('Choose configurations in HISTORY for analysis:')
    print ('[1] All configurations')
    print ('[2] User defined range')
    print ('[3] Custom list')
    user_input = raw_input('Please choose \'1\' or \'2\' or \'3\': ').strip()
    if user_input == '1':
      configs_list = range(1,n_configs+1)
    elif user_input == '2':
      print (GV.dash_line)
      user_input2 = None
      while not user_input2:
        user_input2 = raw_input('Please enter starting configuration number (1-' + str(n_configs) + '): ').strip()
        if user_input2.isdigit():
          if int(user_input2) < 1 or int(user_input2) > n_configs:
            print ('Please enter an integer between 1 and ' + str(n_configs))
            user_input2 = None
            print (GV.dash_line)
          else:
            start_config_no = int(user_input2)
        else:
          print (str(user_input2) + ' is not an integer.')
          user_input2 = None
          print (GV.dash_line)

      user_input2 = None
      while not user_input2:
        user_input2 = raw_input('Please enter final configuration number (' + str(start_config_no) + '-' + str(n_configs) + '): ').strip()
        if user_input2.isdigit():
          if int(user_input2) < start_config_no or int(user_input2) > n_configs:
            print ('Please enter an integer between ' + str(start_config_no) +  ' and ' + str(n_configs))
            user_input2 = None
            print (GV.dash_line)
          else:
            final_config_no = int(user_input2)
        else:
          print (str(user_input2) + ' is not an integer.')
          user_input2 = None
          print (GV.dash_line)
      configs_list = range(start_config_no,final_config_no+1)
    elif user_input == '3':
      configs_list = sorted(CF.chooseCustomConfigurations(n_configs))
    else:
      user_input = None
      print ('Please enter \'1\' or \'2\' or \'3\'.')
      print (GV.dash_line)
  coordinates_list,timestep_list,unit_cell_length_list = RI.genCoordinatesListFromHistory(history_file,mol_info,copy.deepcopy(configs_list),n_atoms,timestep_size)
  all_info = {'coords_list':coordinates_list,\
              'timestep_list':timestep_list,\
              'unit_cell_length_list':unit_cell_length_list,\
              'configs_no_list':configs_list,\
              'mol_info':mol_info,\
              'vdw_info':vdw_info,\
              'n_atoms':n_atoms,\
              'timestep_size':timestep_size\
             }

  ####################################################
  ### Different analysis options available to user ###
  ####################################################
  user_input = None
  while not user_input:
    print ('Available actions:')
    print ('-->[1] Distance of target species within a certain distance threshold to a chosen origin')
    print ('-->[2] Distance between two chosen atoms')
    print ('-->[3] Average distance over time of target species to a chosen origin')
    print ('-->[4] Min distance of target species to a chosen origin over whole trajectory')
    print ('-->[5] Distance between one or more particles to a chosen origin')
    print ('-->[6] Distance from particles to origin for a set of configurations')
    print ('-->[7] Convert configuration(s) from HISTORY into .xyz format')
    user_input = raw_input('Please choose option: ').strip()
    if user_input.isdigit():
      user_input = int(user_input)
      if user_input == 1:
        analyse_target_species_within_distance_threshold(all_info,workdir)
      elif user_input == 2:
        output_distance_between_two_specified_particles(all_info,workdir)
      elif user_input == 3:
        average_distance_to_chosen_origin(all_info,workdir)
      elif user_input == 4:
        min_distance_to_chosen_origin(all_info,workdir)
      elif user_input == 5:
        output_distance_from_particles_to_chosen_origin(all_info,workdir)
      elif user_input == 6:
        output_ordered_distance_to_chosen_origin_for_specified_configurations(all_info,workdir)
      elif user_input == 7:
        convert_configs_into_xyz(all_info,workdir)
      else:
        ('Please enter \'1\' or \'2\' or \'3\' or \'4\' or \'5\' or \'6\' or \'7\'.')
        user_input = None
    else:
      print (user_input + ' is not an integer.')
      user_input = None
  return


#-------------------------------------------------------------------------------
# OPTION 1
# 
# Function analyse_target_species_within_distance_threshold
#-------------------------------------------------------------------------------
def analyse_target_species_within_distance_threshold(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],all_info['vdw_info'])
  
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')
  origin_atom_no = CF.chooseOriginAtomNo(all_info['n_atoms'])

  #####################################################################
  # Choose lower and upper bounds of distance from origin to consider #
  #####################################################################
  user_input = None
  while not user_input:
    user_input = raw_input('Please enter lower bound of distance from origin: ').strip()
    try: 
      lower_bound = float(user_input)
      if lower_bound < 0.0:
        print ('You entered a negative number!')
        user_input = None
        print (GV.dash_line)
    except:
      print (user_input + ' is not a number')
      print (GV.dash_line)
      user_input = None

  user_input = None
  while not user_input:
    user_input = raw_input('Please enter upper bound of distance from origin: ').strip()
    try: 
      upper_bound = float(user_input)
      if upper_bound < lower_bound:
        print ('The upper bound must be larger than the lower bound!')
        print (GV.dash_line)
        user_input = None
    except:
      print (user_input + ' is not a number')
      user_input = None
      print (GV.dash_line)
  
  ##################################################
  # Choose target atom(s)/molecule(s) for analysis #
  ##################################################
  user_input = None
  while not user_input:
    print ('Do you wish to target:')
    print ('-->[A]toms')
    print ('-->[M]olecules')
    user_input = raw_input('Please enter \'A\' or \'M\': ').strip()
    if user_input.upper() == 'A':
      target_type = 'A'
      target_type_string = 'atom'
      target_label = CF.chooseAtomType(first_system)

    elif user_input.upper() == 'M':
      target_type = 'M'
      target_type_string = 'molecule'
      target_label = CF.chooseMolType(first_system)
    else:
      print ('Please enter \'A\' or \'M\'.')
      user_input = None
      print (GV.dash_line)

  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('Do you want to consider ' + target_type_string + 's that started outside the distance threshold but entered the distance threshold at sometime during the trajectory [Y/N]? -->',80) + ' ').strip().upper()
    if user_input == 'Y':
      scan_all_configs = True
    elif user_input == 'N':
      scan_all_configs = False
    else:
      user_input = None
      print ('Please enter \'Y\' or \'N\'.')
      print (GV.dash_line)
  
  target_no_list = []
  if target_type == 'M':
    #output = op.join(workdir,'dist_from_mol_within_threshold_to_origin_vs_t.txt')
    if scan_all_configs == True:
      output = op.join(workdir,'dist_vs_t_from_mol_within_' + str(lower_bound) + '-' + str(upper_bound) + '_to_atom_' + str(origin_atom_no) + '_at_any_time.txt')
      for i,c in enumerate(all_info['configs_no_list']):
        #print ('Analysing configuration ' + str(i+1) + ' of ' + str(len(configuration_list)) + '.')
        print ('Analysing configuration No. ' + str(c) + ' (' + str(i+1) + ' of ' + str(len(all_info['configs_no_list'])) + ').')
        s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],all_info['vdw_info'],all_info['timestep_list'][i],all_info['timestep_list'][i]*all_info['timestep_size'])
        #print (s.atom_list[origin_atom_no-1].coords)
        s.shift_origin_to_coords(s.atom_list[origin_atom_no-1].coords)
        #print (s.atom_list[origin_atom_no-1].coords)
        for j,m in enumerate(s.mol_list):
          if m.label != target_label:
            continue
          elif j+1 in target_no_list:
            continue
          elif m.dist_to_origin > lower_bound and m.dist_to_origin < upper_bound:
            #print m.centre_of_mass
            #print m.dist_to_origin
            target_no_list.append(j+1)
    elif scan_all_config == False:
      output = op.join(workdir,'dist_vs_t_from_mol_within_' + str(lower_bound) + '-' + str(upper_bound) + '_to_atom_' + str(origin_atom_no) + '_at_the_start.txt')
      first_system.shift_origin_to_coords(first_system.atom_list[origin_atom_no-1].coords)
      for j,m in enumerate(first_system.mol_list):
        if m.label != target_label:
          continue
        elif m.dist_to_origin > lower_bound and m.dist_to_origin < upper_bound:
          target_no_list.append(j+1)

    with open(output,'w') as f:
      f.write('t')
      for m in target_no_list:
        f.write(' ' + str(m))
      f.write('\n')
      for i,c in enumerate(all_info['configs_no_list']):
        print ('Writing distances to origin for configuration No.' + str(c) + ' (' +  str(i+1) + ' of ' + str(len(all_info['configs_no_list'])) + ').')
        s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None,all_info['timestep_list'][i],all_info['timestep_list'][i]*all_info['timestep_size'])
        s.shift_origin_to_coords(s.atom_list[origin_atom_no-1].coords)
        f.write(str(i))
        for m in target_no_list:
          f.write(' ' + str(Decimal(s.mol_list[m-1].dist_to_origin).quantize(dp_obj)))
        f.write('\n')
  elif target_type == 'A':
    pass
  print (target_no_list)
  return

#-------------------------------------------------------------------------------
# OPTION 2
#-------------------------------------------------------------------------------
def output_distance_between_two_specified_particles(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],all_info['vdw_info'])
  n_atoms = all_info['n_atoms']
  n_mols = len(first_system.mol_list)
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')
  
  user_input = None
  while not user_input:
    print ('Choose particle 1 type:')
    print ('-->[A] Atom')
    print ('-->[M] Molecule')
    user_input = raw_input('Please enter \'A\' or \'M\' --> ').strip().upper()
    if user_input == 'A':
      particle_1_type = 'A' 
      particle_1_no = CF.chooseAtomNo(n_atoms)
    elif user_input == 'M':
      particle_1_type = 'M' 
      particle_1_no = CF.chooseMolNo(n_mols)
    else:
      print ('Please enter \'A\' or \'M\'.')
      print (GV.dash_line)
      user_input = None

  print (GV.dash_line) 
  user_input = None
  while not user_input:
    print ('Choose particle 2 type:')
    print ('-->[A] Atom')
    print ('-->[M] Molecule')
    user_input = raw_input('Please enter \'A\' or \'M\' --> ').strip().upper()
    if user_input == 'A':
      particle_2_type = 'A' 
      particle_2_no = CF.chooseAtomNo(n_atoms)
    elif user_input == 'M':
      particle_2_type = 'M' 
      particle_2_no = CF.chooseMolNo(n_mols)
    else:
      print ('Please enter \'A\' or \'M\'.')
      print (GV.dash_line)
      user_input = None
  dist_list = [] 
  for i,c in enumerate(all_info['configs_no_list']):
    output = 'dist_vs_t_between_'
    print ('Analysing configurtaion No. ' + str(c) + ' ('+ str(i+1) + ' of ' + str(len(all_info['configs_no_list'])) + ').')
    s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None,all_info['timestep_list'][i],all_info['timestep_list'][i]*all_info['timestep_size'])
    #s = MolClass.BulkSystem(c,mol_info,unit_cell_length,vdw_info)
    if particle_1_type == 'A':
      particle_1 = s.atom_list[particle_1_no-1] 
      output += 'atom_' + str(particle_1_no) + '_and_'
    elif particle_1_type == 'M':
      particle_1 = s.mol_list[particle_1_no-1] 
      output += 'mol_' + str(particle_1_no) + '_and_'
    if particle_2_type == 'A':
      particle_2 = s.atom_list[particle_2_no-1] 
      output += 'atom_' + str(particle_2_no) + '.txt'
    elif particle_2_type == 'M':
      particle_2 = s.mol_list[particle_2_no-1] 
      output += 'mol_' + str(particle_2_no) + '.txt'
    #print (s.atom_list[atom_1_no-1].coords)
    #print (s.atom_list[atom_2_no-1].coords)
    dist_list.append(particle_1.dist_to_particle(particle_2))
  output = op.join(workdir,output)
  with open(output,'w') as f:
    for d in dist_list:
      f.write(str(d) + '\n')
  return

#-------------------------------------------------------------------------------
# OPTION 3
#-------------------------------------------------------------------------------
def average_distance_to_chosen_origin(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],None)
  origin_atom_no = CF.chooseOriginAtomNo(len(first_system.atom_list))
  dist_avg_list = []
  n_systems = len(all_info['configs_no_list'])
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')
  user_input = None
  while not user_input:
    print ('Do you wish to target:')
    print ('-->[A]toms')
    print ('-->[M]olecules')
    user_input = raw_input('Please enter \'A\' or \'M\': ').strip().upper()
    if user_input == 'A':
      target_type = 'A'
      target_type_str = 'atoms'
    elif user_input == 'M':
      target_type = 'M'
      target_type_str = 'mols'
    else:
      print ('Please enter \'A\' or \'M\'.')
      user_input = None
      print (GV.dash_line)
  for i,c in enumerate(all_info['configs_no_list']):
    print ('Analysing configurtaion No. ' + str(c) + ' (' +  str(i+1) + ' of ' + str(n_systems) + ').')
    s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None,all_info['timestep_list'][i],all_info['timestep_list'][i]*all_info['timestep_size'])
    origin_coords = s.atom_list[origin_atom_no-1].coords
    if target_type == 'A':
      if i == 0:
        dist_avg_list = [0]*len(s.atom_list)
      for j,a in enumerate(s.atom_list):
        dist_avg_list[j] += a.dist_to_coords(origin_coords)
    elif target_type == 'M':
      if i == 0:
        dist_avg_list = [0]*len(s.mol_list)
      for j,m in enumerate(s.mol_list):
        dist_avg_list[j] += m.dist_to_coords(origin_coords)
  for p,d in enumerate(dist_avg_list):
    dist_avg_list[p] = Decimal(str(d/n_systems)).quantize(dp_obj)
  dist_avg_list = zip(range(1,len(dist_avg_list)+1),dist_avg_list)
  unordered_output = op.join(workdir,'avg_unordered_dist_of_all_' + target_type_str + '_to_atom_' + str(origin_atom_no) +'.txt')
  with open(unordered_output,'w') as f:
    for i in dist_avg_list:
      f.write(str(i[0]) + '   ' + str(i[1]) + '\n')
  dist_avg_list = sorted(dist_avg_list,key=lambda x: x[1])
  ordered_output = op.join(workdir,'avg_ordered_dist_of_all_' + target_type_str + '_to_atom_' + str(origin_atom_no) +'.txt')
  with open(ordered_output,'w') as f:
    for i in dist_avg_list:
      f.write(str(i[0]) + '   ' + str(i[1]) + '\n')
  return

#-------------------------------------------------------------------------------
# OPTION 4
#-------------------------------------------------------------------------------
def min_distance_to_chosen_origin(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],None)
  origin_atom_no = CF.chooseOriginAtomNo(len(first_system.atom_list))
  dist_min_list = []
  n_systems = len(all_info['configs_no_list'])
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')

  user_input = None
  while not user_input:
    print ('Do you wish to target:')
    print ('-->[A]toms')
    print ('-->[M]olecules')
    user_input = raw_input('Please enter \'A\' or \'M\': ').strip().upper()
    if user_input == 'A':
      target_type = 'A'
      unordered_output = op.join(workdir,'min_unordered_dist_from_all_atoms_to_atom_' + str(origin_atom_no) + '.txt')
      ordered_output = op.join(workdir,'min_ordered_dist_from_all_atoms_to_atom_' + str(origin_atom_no) + '.txt')
    elif user_input == 'M':
      target_type = 'M'
      unordered_output = op.join(workdir,'min_unordered_dist_from_all_mols_to_atom_' + str(origin_atom_no) + '.txt')
      ordered_output = op.join(workdir,'min_ordered_dist_from_all_mols_to_atom_' + str(origin_atom_no) + '.txt')
    else:
      print ('Please enter \'A\' or \'M\'.')
      user_input = None
      print (GV.dash_line)
  for i,c in enumerate(all_info['configs_no_list']):
    print ('Analysing configurtaion ' + str(i+1) + ' of ' + str(n_systems) + '.')
    s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None)
    origin_coords = s.atom_list[origin_atom_no-1].coords
    if target_type == 'A':
      if i == 0:
        dist_min_list = [s.unit_cell_length*10]*len(s.atom_list)
      for j,a in enumerate(s.atom_list):
        dist = a.dist_to_coords(origin_coords)
        if dist < dist_min_list[j]:
          dist_min_list = dist
    elif target_type == 'M':
      if i == 0:
        dist_min_list = [s.unit_cell_length*10]*len(s.mol_list)
      for j,m in enumerate(s.mol_list):
        dist = m.dist_to_coords(origin_coords)
        if dist < dist_min_list[j]:
          dist_min_list[j] = dist
  for i,d in enumerate(dist_min_list):
    dist_min_list[i] = Decimal(str(d)).quantize(dp_obj)
  dist_min_list = zip(range(1,len(dist_min_list)+1),dist_min_list)
  
  with open(unordered_output,'w') as f:
    for i in dist_min_list:
      f.write(str(i[0]) + '   ' + str(i[1]) + '\n')
  dist_min_list = sorted(dist_min_list,key=lambda x: x[1])
  with open(ordered_output,'w') as f:
    for i in dist_min_list:
      f.write(str(i[0]) + '   ' + str(i[1]) + '\n')
  return


#-------------------------------------------------------------------------------
# OPTION 5
#-------------------------------------------------------------------------------
def output_distance_from_particles_to_chosen_origin(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],None)
  origin_atom_no = CF.chooseOriginAtomNo(len(first_system.atom_list))
  dist_min_list = []
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')
  n_systems = len(all_info['configs_no_list'])

  user_input = None
  while not user_input:
    print ('Do you wish to target:')
    print ('-->[A]toms')
    print ('-->[M]olecules')
    user_input = raw_input('Please enter \'A\' or \'M\': ').strip().upper()
    if user_input == 'A':
      target_type = 'A'
      atom_no_set = CF.chooseAtomNoSet(first_system)
    elif user_input == 'M':
      target_type = 'M'
      mol_no_set = CF.chooseMolNoSet(first_system)
    else:
      print ('Please enter \'A\' or \'M\'.')
      user_input = None
  if target_type == 'A':
    output_file = op.join(workdir,'dist_vs_t_from_atom_set_to_atom_' + str(origin_atom_no) + '.txt')
    with open(output_file,'w') as f:
      f.write('t')
      for a in atom_no_set:
        f.write(' ' + str(a))
      f.write('\n')
  elif target_type == 'M':
    output_file = op.join(workdir,'dist_vs_t_from_mol_set_to_atom_' + str(origin_atom_no) + '.txt')
    with open(output_file,'w') as f:
      f.write('t')
      for m in mol_no_set:
        f.write(' ' + str(m))
      f.write('\n')
  with open(output_file,'a') as f:
    for i,c in enumerate(all_info['configs_no_list']):
      print ('Analysing configuration No. ' + str(c) + ' (' + str(i+1) + ' of ' + str(n_systems) + ').')
      f.write(str(i+1))
      s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None)
      origin_coords = s.atom_list[origin_atom_no-1].coords
      print ('origin_coords: ' + str(origin_coords))
      if target_type == 'A':
        for a in atom_no_set:
          d = s.atom_list[a-1].dist_to_coords(origin_coords)
          f.write(' ' + str(Decimal(str(d)).quantize(dp_obj))) 
        f.write('\n')
      elif target_type == 'M':
        for m in mol_no_set:
          d = s.mol_list[m-1].dist_to_coords(origin_coords)
          f.write(' ' + str(Decimal(str(d)).quantize(dp_obj))) 
        f.write('\n')
  return


#-------------------------------------------------------------------------------
# OPTION 6
#-------------------------------------------------------------------------------
def output_ordered_distance_to_chosen_origin_for_specified_configurations(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],None)
  origin_atom_no = CF.chooseOriginAtomNo(len(first_system.atom_list))
  dist_min_list = []
  dp_obj = Decimal('0.' + (decimal_places-1)*'0' + '1')
  n_systems = len(all_info['configs_no_list'])

  user_input = None
  while not user_input:
    print ('Do you wish to target:')
    print ('-->[A]toms')
    print ('-->[M]olecules')
    user_input = raw_input('Please enter \'A\' or \'M\': ').strip().upper()
    if user_input == 'A':
      target_type = 'A'
      target_type_str = 'atoms'
      max_particles = len(first_system.atom_list)
    elif user_input == 'M':
      target_type = 'M'
      target_type_str = 'molecules'
      max_particles = len(first_system.mol_list)
    else:
      print ('Please enter \'A\' or \'M\'.')
      user_input = None
      print (GV.dash_line)
  
  print (GV.dash_line)
  user_input = None
  while not user_input:
    print ('Do you wish to include:')
    print ('-->[1] A defined number of particles')
    print ('-->[2] All particles within a distance threshold')
    user_input = raw_input('Please enter \'1\' or \'2\': ').strip().upper()
    if user_input == '1':
      analysis_type = 'N'
    elif user_input == '2':
      analysis_type = 'D'
    else:
      print ('Please enter \'1\' or \'2\'.')
      user_input = None
      print (GV.dash_line)

  if analysis_type == 'N':
    print (GV.dash_line)
    user_input = None
    while not user_input:
      user_input = raw_input('Include how many of the closest ' + target_type_str + '(<' + str(max_particles-1) + ')? ')
      if user_input.isdigit():
        user_input = int(user_input)
        if user_input < 1:
          print ('Please enter an integer between 1 and ' + str(max_particles) + '.')
          user_input = None
          print (GV.dash_line)
        elif user_input > max_particles:
          print ('Please enter an integer between 1 and ' + str(max_particles) + '.')
          user_input = None
          print (GV.dash_line)
        else:
          n_closest_particles = user_input
      else:
        print (user_input + ' is not an integer.')
        user_input = None
        print (GV.dash_line)

  elif analysis_type == 'D':
    print (GV.dash_line)
    user_input = None
    while not user_input:
      user_input = raw_input('Include all partcles within what threshold? ')
      try: 
        user_input = float(user_input)
        if user_input < 0.0:
          print ('Please enter a number greater than 0.')
          user_input = None
          print (GV.dash_line)
        else:  
          dist_thres = user_input
      except:
        print (user_input + ' is not a number.')
        user_input = None
        print (GV.dash_line)

        
  if target_type == 'A':
    if analysis_type == 'N':
      output_file = op.join(workdir,'closest_' + str(n_closest_particles) + '_atoms_to_atom_' + str(origin_atom_no) + '_for_each_config.txt')
      heading_n = n_closest_particles
    elif analysis_type == 'D':
      output_file = op.join(workdir,'atoms_within_' + str(dist_thres) + 'A_to_atom_' + str(origin_atom_no) + '_for_each_config.txt')
      heading_n = len(first_system.atom_list)
  elif target_type == 'M':
    if analysis_type == 'N':
      output_file = op.join(workdir,'closest_' + str(n_closest_particles) + '_molecules_to_atom_' + str(origin_atom_no) + '_for_each_config.txt')
      heading_n = n_closest_particles
    elif analysis_type == 'D':
      output_file = op.join(workdir,'molecules_within_' + str(dist_thres) + 'A_to_atom_' + str(origin_atom_no) + '_for_each_config.txt')
      heading_n = len(first_system.mol_list)
  with open(output_file,'w') as f:
    f.write('t')
    if analysis_type == 'D':
      f.write(' n')
    for i in range(heading_n):
      f.write(' ' + target_type + str(i+1) + ' dist2origin')
    f.write('\n')

  with open(output_file,'a') as f:
    for i,c in enumerate(all_info['configs_no_list']):
      dist_list = []
      print ('Analysing configuration No. ' + str(c) + ' (' + str(i+1) + ' of ' + str(n_systems) + ').')
      f.write(str(i+1))
      s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None)
      origin_coords = s.atom_list[origin_atom_no-1].coords
      print ('origin_coords: ' + str(origin_coords))
      if target_type == 'A':
        for a,a_obj in enumerate(s.atom_list):
          if a == origin_atom_no - 1:
            continue
          d = a_obj.dist_to_coords(origin_coords)
          if analysis_type == 'D' and d > dist_thres:
            continue
          if len(dist_list) == 0:
            dist_list.append((a+1,d))
          else: 
            for x,y in enumerate(dist_list):
              if d < y[1]:
                dist_list.insert(x,(a+1,d))
                break
              elif x == len(dist_list)-1:
                dist_list.append((a+1,d))
                break
          if analysis_type == 'N' and len(dist_list) > n_closest_particles:
            del dist_list[-1]
      elif target_type == 'M':
        for m,m_obj in enumerate(s.mol_list):
          d = m_obj.dist_to_coords(origin_coords)
          if analysis_type == 'D' and d > dist_thres:
            continue
          if len(dist_list) == 0:
            dist_list.append((m+1,d))
          else: 
            for x,y in enumerate(dist_list):
              if d < y[1]:
                dist_list.insert(x,(m+1,d))
                break
              elif x == len(dist_list)-1:
                dist_list.append((m+1,d))
                break
          if analysis_type == 'N' and len(dist_list) > n_closest_particles:
            del dist_list[-1]
      if analysis_type == 'D':
        f.write(' ' + str(len(dist_list)))
      for info in dist_list:
        f.write(' ' + str(info[0]) + ' ' + str(Decimal(str(info[1])).quantize(dp_obj)))
      f.write('\n')

  return

#-------------------------------------------------------------------------------
# OPTION 7
#-------------------------------------------------------------------------------
def convert_configs_into_xyz(all_info,workdir,decimal_places=4):
  first_system = MolClass.BulkSystem(all_info['coords_list'][0],all_info['mol_info'],all_info['unit_cell_length_list'][0],all_info['vdw_info'])
  n_mols = len(first_system.mol_list)
  n_systems = len(all_info['configs_no_list'])
  user_input = None 
  while not user_input:
    print ('Please choose new origin:')
    print ('-->[1] Use old origin')
    print ('-->[2] Atom')
    print ('-->[3] Molecule')
    user_input = raw_input('Please choose \'1\', \'2\' or \'3\': ').strip()
    if user_input == '1':
      origin_type = None
    elif user_input == '2':
      origin_type = 'A' 
      origin_no = CF.chooseAtomNo(all_info['n_atoms'])
    elif user_input == '3':
      origin_type = 'M' 
      origin_no = CF.chooseMolNo(n_mols)
    else:
      user_input = None
      print ('Please enter \'1\' or \'2\' or \'3\'.')
      print (GV.dash_line)
  print (GV.dash_line)
  
  user_input = None
  while not user_input:
    print ('Output options:')
    print ('-->[1] Single file containing all configurations')
    print ('-->[2] Individual files for each configuration')
    user_input = raw_input('Please choose \'1\', \'2\': ').strip()
    if user_input == '1':
      single_file = True 
    elif user_input == '2':
      single_file = False
    else:
      user_input = None
      print ('Please enter \'1\' or \'2\'.')
      print (GV.dash_line)
  

  for i,c in enumerate(all_info['configs_no_list']):
    print ('Generating output for configuration No. ' + str(c) + ' (' + str(i+1) + ' of ' + str(n_systems) + ').')
    s = MolClass.BulkSystem(all_info['coords_list'][i],all_info['mol_info'],all_info['unit_cell_length_list'][i],None,all_info['timestep_list'][i],(all_info['timestep_list'][i]-1)*all_info['timestep_size'])
    if origin_type == 'M':
      s.shift_origin_to_coords(s.mol_list[origin_no-1].coords)
    elif origin_type == 'A':
      s.shift_origin_to_coords(s.atom_list[origin_no-1].coords)

    title = 'Configuration ' + str(c) + '; time=' + str(s.time)
    if single_file == True:
      output = op.join(workdir,'configuration_set.xyz')
      if i == 0:
        append = False
      else:
        append = True
    elif single_file == False:
      output = op.join(workdir,'configuration_' + str(c) + '.xyz')
      append = False
    s.write_xyz(output,title,False,append) 
  return

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
if __name__ == "__main__":
  #print (sys.argv)
  workdir = sys.argv[1]
  sys.exit(main(workdir))

