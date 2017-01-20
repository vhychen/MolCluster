import textwrap
import GlobalVariables as GV
import MolClass as MC
from decimal import Decimal
import math

#-------------------------------------------------------------------------------
# Function printBugsLine
# 
# Operation: Prints out the 'bugs' line
#-------------------------------------------------------------------------------
def printBugsLine(temp_log_file):
  bugs_line = 'Please report any bugs to Vincent Chen (vc0489@gmail.com). The log file can be found in \'' + str(temp_log_file) + '\'.'
  print (textwrap.fill(bugs_line,80))
  return


#-------------------------------------------------------------------------------
# Function coords_to_string
#
# Operation: Converts a set of coordinates into string format
#-------------------------------------------------------------------------------
def coords_to_string(coords):
  string = '('
  for i in coords:
    string += str(i) + ','
  string = string[:-1] + ')'
  return (string)
  
#-------------------------------------------------------------------------------
# Function choose_QM_regions
#
# Operation: Interface to choose a custom set of one or more QM regions.
#-------------------------------------------------------------------------------
def choose_qm_regions(cluster_obj,log_file=None):
  n_regions = None
  while not n_regions:
    n_regions = raw_input('Enter the number of QM regions you would like to define. --> ').strip()
    if n_regions.isdigit() == False:
      print (n_regions + ' is not a number')
      n_regions = None
      print (GV.dash_line)
    elif n_regions < 1:
      print ('Please enter an integer greater than 0')
      n_regions = None
      print (GV.dash_line)
    else:
      n_regions = int(n_regions)
  if log_file:
    with open(log_file,'a') as f: 
      f.write('The QM region will consist of ' + str(n_regions) + ' spherical regions in total.\n')
  print (GV.dash_line)
  regions_def = []
  while len(regions_def) < n_regions:
    origin_type = None
    while not origin_type:
      print ('QM region no. ' + str(len(regions_def)+1) + '. Possible origin types:')
      print ('-->[A]tom')
      print ('-->[M]olecule')
      print ('-->[C]luster origin')
      origin_type = raw_input('Please choose origin type [A/M/C] --> ').strip()
      if origin_type == 'A':
        #ordered_atom_list = cluster_obj.ret_ordered_atom_list()
        origin_atom_no = cluster_obj.choose_atom_origin(len(regions_def)+1)
        origin_coords = cluster_obj.atom_list[origin_atom_no-1].coords
      elif origin_type == 'M':
        cluster_obj.print_mol_list(cluster_obj.name)
        origin_mol_no = None
        while not origin_mol_no:
          origin_mol_no = raw_input(textwrap.fill('Please choose the molecule to be used as the origin of QM region no. ' + str(len(regions_def)+1) + '. Note that the centre of mass of the molecule will be used. --> ',80) + ' ').strip()
          if not origin_mol_no.isdigit():
            print (origin_mol_no + ' is not valid. Please enter an integer between 1 and ' + str(len(cluster_obj.mol_list)) + '.')
            print (GV.dash_line)
            origin_mol_no = None
          elif int(origin_mol_no) < 1 or int(origin_mol_no) > len(cluster_obj.mol_list):
            print (origin_mol_no + ' is not valid. Please enter an integer between 1 and ' + str(len(cluster_obj.mol_list)) + '.')
            print (GV.dash_line)
            origin_mol_no = None
          else:
            origin_mol_no = int(origin_mol_no)
            origin_coords = cluster_obj.mol_list[origin_mol_no-1].centre_of_mass
      elif origin_type == 'C':
        origin_coords = (0.0,0.0,0.0)
      else:
        print (origin_type + ' is not a valid option')
        print (GV.dash_line)
        origin_type = None
        continue 
      qm_radius = choose_qm_radius(cluster_obj.radius,cluster_obj)
      if log_file:
        with open(log_file,'a') as f:
          if origin_type == 'A':
            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about atom ' + str(origin_atom_no) + ' [type: ' + cluster_obj.atom_list[origin_atom_no-1].label + '; coords: ' + coords_to_string(origin_coords) + '] with a radius of ' + str(qm_radius) + ' ' + u'\u212B'.encode('utf-8') + '.\n')
          elif origin_type =='M':
            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about the centre of mass of molecule ' + str(origin_mol_no) + ' [type: ' + cluster_obj.mol_list[origin_mol_no-1].label + '; coords: ' + coords_to_string(origin_coords) + '] with a radius of ' + str(qm_radius) + ' ' + u'\u212B'.encode('utf-8') + '.\n')
          elif origin_type == 'C':
            f.write('-->QM region no. ' + str(len(regions_def) + 1) + ' is centred about the cluster origin with a radius of ' + str(qm_radius) + u'\u212B'.encode('utf-8') + '.\n')
    regions_def.append((origin_coords,qm_radius))
  return (regions_def)

#-------------------------------------------------------------------------------
# Function define_region_of_number_of_each_molecule_type
#
# Operation: Allows the user to define the number of each type of molecule in a
#            specific region.
#-------------------------------------------------------------------------------
def define_region_of_number_of_each_molecule_type(region_string,system,system_type_title_str='',calc_est_radius=False,log_file=None):
  region_mol_numbers_list = []
  region_tot_charge = 0.0
  #region_tot_molecules = 0
  region_tot_atoms = 0
  while not region_mol_numbers_list:
    print ('')
    system.print_system_info(system_type_title_str)
    for m in system.mol_info:
      user_input = None
      while not user_input: 
        user_input = raw_input('Please enter the number of ' + m['molname'] + ' molecules you want in each ' + str(region_string) + ': ').strip()
        if user_input.isdigit() == False or int(user_input) < 0:
          print (user_input + ' is not an integer >= 0.')
          print (GV.dash_line)
          user_input = None
        #elif int(user_input) > m['nmols']
        #  user_input = raw_input('WARNING: There are only ' + str(m['nmols']) + ' in the periodic unit cell. A supercell will be created. Do you wish to continue [Y/N]?')
        else:
          region_mol_numbers_list.append(int(user_input))
          region_tot_charge += int(user_input)*sum(m['chargelist'])
          region_tot_atoms += int(user_input)*m['natoms']
      user_input = None
    if calc_est_radius == True:
      est_radius = math.pow(3.0/(4.0*math.pi)*region_tot_atoms/system.density,1.0/3.0)
    print (GV.dash_line)
    print ('Each ' + region_string + ' will be formed of:')
    for i,m in enumerate(system.mol_info):
      print ('--> ' + str(region_mol_numbers_list[i]) + ' X ' + m['molname'])
    print ('')
    print ('I'*80)
    print ('III' + ('  There will be ' + str(sum(region_mol_numbers_list)) + ' molecules (' + str(region_tot_atoms) + ' atoms) in each ' + region_string + '.  ').center(74) + 'III')
    print ('III' + ('  The total charge of each ' + region_string +  ' will be ' + str(Decimal(region_tot_charge).quantize(Decimal('.01'))) + ' e.  ').center(74) + 'III')
    if calc_est_radius == True:
      print ('III' + ('  The estimated radius of each ' + region_string + ' is ~' + str(Decimal(est_radius).quantize(Decimal('.1'))) + ' ' + GV.ANGSTROM_STR + '.  ').center(76) + 'III') 
    print ('I'*80)
    print ('')

    second_user_input = None
    while not second_user_input:
      second_user_input = raw_input('Do you wish to continue [Y/N]? ').strip()
      if second_user_input.upper() == 'Y':
        pass
      elif second_user_input.upper() == 'N':
        region_mol_numbers_list = []
        region_tot_atoms = 0
      else:
        second_user_input = None
        print ('Please enter \'Y\' or \'N\'')
        print (GV.dash_line)
  if log_file:
    with open(log_file,'a') as f:
      f.write('Each ' + region_string + ' to be formed of:\n')
      for i,m in enumerate(system.mol_info):
        f.write('--> ' + str(region_mol_numbers_list[i]) + ' X ' + m['molname'] + '\n')
      f.write('Total of ' + str(sum(region_mol_numbers_list)) + ' molecules (' + str(region_tot_atoms) + ' atoms) in each ' + region_string + '.\n')
      f.write('Total charge of each ' + region_string + ' is ' + str(Decimal(region_tot_charge).quantize(Decimal('.01'))) + ' e.\n')
      if calc_est_radius == True:
        f.write('The estimated radius of each ' + region_string + ' is ~' + str(Decimal(est_radius).quantize(Decimal('.1'))) + ' ' + GV.ANGSTROM_STR + '.\n')
  if calc_est_radius == True:
    return (region_mol_numbers_list,est_radius)
  else:
    return (region_mol_numbers_list)

#-------------------------------------------------------------------------------
# Function choose_region_radius
# 
# Operation: Prompts user to choose the radius of a region within a cluster.
#-------------------------------------------------------------------------------
def choose_region_radius(cluster_radius,system,region_type_str,all_option = False):
  print ('')
  region_radius = None
  while not region_radius:
    if cluster_radius:
      system.print_estimate_n_atoms_in_cluster(3,cluster_radius,1.0)
    else:
      system.print_estimate_n_atoms_in_cluster(3,10.0,0.5)
    if all_option == True:
      region_radius = raw_input(textwrap.fill('Please enter radius of ' + region_type_str + ' region in Angstroms. Note that a molecule will be included in region if its centre of mass is within the ' + region_type_str + ' boundary. If you want all molecules to be in the ' + region_type_str + ' region enter \'ALL\' -->',80) + ' ').strip()
    else:
      region_radius = raw_input(textwrap.fill('Please enter radius of ' + region_type_str + ' region in Angstroms. Note that a molecule will be included in region if its centre of mass is within the ' + region_type_str + ' boundary. -->',80) + ' ').strip()
    try:
      region_radius = float(region_radius)
      if cluster_radius and region_radius > cluster_radius:
        #print ('You requested a radius greater than the cluster radius (' + str(cluster_radius) + ').')
        print ('')
        print (textwrap.fill('***WARNING: You requested a ' + region_type_str + ' radius which may be larger than the cluster radius.***',80))
        print ('')
        print (GV.dash_line)
        #region_radius = None
      elif region_radius < 0:
        print ('You requested a negative radius!')
        print (GV.dash_line)
        region_radius = None
    except:
      if all_option == True and region_radius.upper() == 'ALL':
        region_radius = 'A'
      else:
        print (str(region_radius) + ' is not a number.')
        print (GV.dash_line)
        region_radius = None
  return (region_radius)



#-------------------------------------------------------------------------------
# Function choose_qm_radius
#
# Operation: Interface to choose the radius of a QM region.
#-------------------------------------------------------------------------------
def choose_qm_radius(cluster_radius,system,qm_region_no=None):
  qm_radius = choose_region_radius(cluster_radius,system,'QM')
  return (qm_radius)

#-------------------------------------------------------------------------------
# Function: chooseActiveRegionDef
#
# Operation: Interface to define the active region.
#-------------------------------------------------------------------------------
def chooseActiveRegionDef(cluster_radius,system,log_file):
  active_region_type = None
  active_region_var = None
  while not active_region_type:
    print ('Available active region definitions:')
    print ('-->[S] Spherical region centred about the origin')
    print ('-->[M] By number of total molcules closest to the origin')
    print ('-->[A] By an upper limit of total atoms closest to the origin')
    print ('-->[N] Specific numbers of each molecule type')
    active_region_type = raw_input('Please choose active region definition type [S/N] --> ').strip()
    if active_region_type.upper() == 'S':
      active_region_type = 'S'
      print (GV.dash_line)
      active_radius = choose_region_radius(cluster_radius,system,'active',True)
      active_region_var = active_radius
      if log_file:
        with open(log_file,'a') as f:
          if active_radius == 'A':
            f.write('User requested all atoms to be included in the active region.\n')
          else:
            f.write('User requested the active region to be defined as a spherical region with radius ' + str(active_radius) + ' ' + GV.ANGSTROM_STR + ' centred about the origin.\n')
    elif active_region_type.upper() == 'M':
      active_region_type = 'M'
      n_active_mol = None
      print (GV.dash_line)
      while not n_active_mol:
        n_active_mol = raw_input(textwrap.fill('Please enter the number of molecules to be included in the active region. Note that the molecules closest to the origin will be included. -->',80) + ' ').strip()
        if n_active_mol.isdigit() == False:
          print (n_active_mol + ' is not a positive integer.')
          print (GV.dash_line)
          n_active_mol = None
        elif int(n_active_mol) < 1:
          print ('Please enter an integer greater than 0.') 
          print (GV.dash_line)
          n_active_mol = None
        else:
          n_active_mol = int(n_active_mol)
          #print (GV.dash_line)
          active_region_var = n_active_mol
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the active region to be defined as the ' + str(n_active_mol) + ' molecules closest to the origin.\n')
    elif active_region_type.upper() == 'A':
      active_region_type = 'A'
      n_active_atoms = None
      print (GV.dash_line)
      while not n_active_atoms:
        n_active_atoms = raw_input(textwrap.fill('Please enter the upper limit for the number of atoms to be included in the active region. Note that only whole molecules will be considered and those closest to the origin will be included. -->',80) + ' ').strip()
        if n_active_atoms.isdigit() == False:
          print (n_active_atoms + ' is not a positive integer.')
          print (GV.dash_line)
          n_active_atoms = None
        elif int(n_active_atoms) < 1:
          print ('Please enter an integer greater than 0.')
          print (GV.dash_line)
          n_active_atoms = None
        else:
          n_active_atoms = int(n_active_atoms)
          active_region_var = n_active_atoms
          print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the active region to be defined as all the molecules closest to the origin without exceeding the upper limit of ' + str(n_active_atoms) + ' atoms.\n')
    elif active_region_type.upper() == 'N':
      active_region_type = 'N'
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the active region to be defined as specific number of molecules closest to the origin for each molecule type.\n')
      #print (GV.dash_line)
    else:
      print (active_region_type + ' is not a valid option.')
      active_region_type = None
      print (GV.dash_line)
  return (active_region_type,active_region_var)

#-------------------------------------------------------------------------------
# Function chooseQMRegionDef
#
# Operation: Interface to define QM region.
#-------------------------------------------------------------------------------
def chooseQMRegionDef(cluster_radius,system,log_file):
  qm_region_type = None
  qm_region_var = None
  print (GV.dash_line)
  while not qm_region_type:
    print ('Available QM region definitions:')
    print ('-->[S] Spherical region centred about the origin')
    print ('-->[M] By number of total molcules closest to the origin')
    print ('-->[A] By an upper limit of total atoms closest to the origin')
    print ('-->[N] Specific number of molecules closest to the origin of each molecule type')
    print ('-->[C] Custom region')
    print ('-->[O] Spherical volumes centred about one or more user-defined origins')
    qm_region_type = raw_input('Please choose QM region definition type [S/M/A/N/C/O] --> ').strip()
    if qm_region_type.upper() == 'S':
      qm_region_type = 'S'
      qm_radius = None
      print (GV.dash_line)
      qm_radius = choose_qm_radius(cluster_radius,system)
      qm_region_var = qm_radius
      print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be defined as a spherical region with radius ' + str(qm_radius) + ' ' + GV.ANGSTROM_STR + ' centred about the origin.\n')
    elif qm_region_type.upper() == 'M':
      qm_region_type = 'M'
      n_qm_mol = None
      print (GV.dash_line)
      while not n_qm_mol:
        n_qm_mol = raw_input(textwrap.fill('Please enter the number of molecules to be included in the QM region. Note that the molecules closest to the origin will be included. -->',80) + ' ').strip()
        if n_qm_mol.isdigit() == False:
          print (n_qm_mol + ' is not a positive integer.')
          print (GV.dash_line)
          n_qm_mol = None
        elif int(n_qm_mol) < 1:
          print ('Please enter an integer greater than 0.') 
          print (GV.dash_line)
          n_qm_mol = None
        else:
          n_qm_mol = int(n_qm_mol)
          qm_region_var = n_qm_mol
          print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be defined as the ' + str(n_qm_mol) + ' molecules closest to the origin.\n')
    elif qm_region_type.upper() == 'A':
      qm_region_type = 'A'
      n_qm_atoms = None
      print (GV.dash_line)
      while not n_qm_atoms:
        n_qm_atoms = raw_input(textwrap.fill('Please enter the upper limit for the number of atoms to be included in the QM region. Note that only whole molecules will be considered and those closest to the origin will be included. -->',80) + ' ').strip()
        if n_qm_atoms.isdigit() == False:
          print (n_qm_atoms + ' is not a positive integer.')
          print (GV.dash_line)
          n_qm_atoms = None
        elif int(n_qm_atoms) < 1:
          print ('Please enter an integer greater than 0.')
          print (GV.dash_line)
          n_qm_atoms = None
        else:
          n_qm_atoms = int(n_qm_atoms)
          qm_region_var = n_qm_atoms
          print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be defined as all molecules closest to the origin without exceeding the upper limit of ' + str(n_qm_atoms) + ' atoms.\n')
    elif qm_region_type.upper() == 'N':
      qm_region_type = 'N'
      print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be defined as specific number of molecules closest to the origin for each molecule type.\n')
    elif qm_region_type.upper() == 'C':
      qm_region_type = 'C'
      print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be customly defined for each cluster.\n') 
    elif qm_region_type in ('O','o'):
      qm_region_type = 'O'
      print (GV.dash_line)
      if log_file:
        with open(log_file,'a') as f:
          f.write('User requested the QM region to be one of more spherical volumes centred about user-defined origin(s).\n') 
    else:
      print (qm_region_type + ' is not a valid option.')
      qm_region_type = None
      print (GV.dash_line)
  return (qm_region_type,qm_region_var)

#-------------------------------------------------------------------------------
# Function chooseOrigin
#
# Operation: Interface to select the origin of the system.
#-------------------------------------------------------------------------------
def chooseOrigin(options,log_file,system):
  options_string_dict = {'AN':'Atom number',\
                         'AT':'Atom type',\
                         'MN':'Molecule number',\
                         'MT':'Molecule type',\
                         'C':'Custom co-ordinates'}
  if not options:
    options = ['AN','AT','MN','MT','C']
  user_input = None 
  n_atoms = len(system.atom_list)
  n_mols = len(system.mol_list)
  while not user_input:
    print ('Available origin types:')
    for o in options:
      print ('-->[' + o + '] ' + options_string_dict[o])
    user_input = raw_input('Please choose origin type --> ').strip().upper()
    if user_input not in options:
      print (user_input + ' is not a valid option')
      user_input = None
      print (GV.dash_line)
    elif user_input == 'AN':
      second_user_input = None
      print (GV.dash_line)
      while not second_user_input:
        second_user_input = raw_input(textwrap.fill('Please enter the number of the atom which will be the origin [integer between 1-' + str(n_atoms) + ']:',80) + ' ').strip()
        try:
          second_user_input = int(second_user_input)
          if second_user_input < 1 or second_user_input > n_atoms:
            print ('You did not enter an interger between 1 and ' + str(n_atoms) + '.') 
            print (GV.dash_line)
            second_user_input = None
          else:
            print ('I'*80)
            print ('III' + ('  Atom ' + str(second_user_input) + ' (' + system.atom_list[second_user_input-1].label + ') will be the origin.  ').center(74) + 'III')
            print ('I'*80)
            origin_type = user_input
            origin_type_var = second_user_input
            with open(log_file,'a') as f:
              f.write('User chose atom no. ' + str(origin_type_var) + ' to be the origin.\n')
              f.write('Atom ' + str(second_user_input) + ' (' + system.atom_list[second_user_input-1].label + ') will be the origin.\n')
        except:
          print (str(second_user_input) + ' is not an integer.')
          print (GV.dash_line)
          second_user_input = None

    elif user_input == 'AT':
      print (GV.dash_line)
      prompt_message = 'Please choose the type of atom for the origin (the first atom of that type will be the origin)'
      origin_atom_type = chooseAtomType(system,prompt_message)
      origin_type = user_input
      origin_type_var = origin_atom_type
      for i,a in enumerate(system.atom_list):
        if a.label == origin_type_var:
          print ('')
          print ('I'*80)
          print ('III' + ('  Atom ' + str(i+1) + ' (' + origin_type_var + ') will be the origin.  ').center(74) + 'III') 
          print ('I'*80)
          print ('')
          with open(log_file,'a') as f:
            f.write('User chose origin to be first atom of type ' + origin_type_var + '.\n')
            f.write('Atom ' + str(i+1) + ' (' + origin_type_var + ') will be the origin.\n')
          break
    elif user_input == 'MN':
      second_user_input = None
      print (GV.dash_line)
      while not second_user_input:
        print ('Molecule numbers:')
        c_mols = 0
        for t in mol_info:
          if t['nmols'] == 1:
            print ('[' + str(c_mols+1) + '] ' + t['molname'])
          else:
            print ('[' + str(c_mols+1) + '-' + str(c_mols+t['nmols']) + '] ' + t['molname'])
          c_mols += t['nmols']
        second_user_input = raw_input(textwrap.fill('Please enter the number of the molecule whose centre of mass will be the origin [integer between 1-' + str(c_mols) + ']:',80) + ' ').strip()
        try:
          second_user_input = int(second_user_input)
          if second_user_input < 1 or second_user_input > c_mols:
            print ('You did not enter an interger between 1 and ' + str(c_mols) + '.') 
            print (GV.dash_line)
            second_user_input = None
          else:
            print ('')
            print ('I'*80)
            print ('III' + ('  The centre of mass of molcule ' + str(second_user_input) + ' (' + system.mol_list[second_user_input-1].label + ') will be the origin.  ').center(74) + 'III')
            print ('I'*80)
            print ('')
            origin_type = user_input
            origin_type_var = second_user_input
            with open(log_file,'a') as f:
              f.write('User chose origin to be the centre of mass of molecule no. ' + str(origin_type_var) + '.\n')
              f.write('The centre of mass of molecule ' + str(second_user_input) + ' (' + system.mol_list[origin_type_var-1].label + ') will be the origin.\n')
        except:
          print (str(second_user_input) + ' is not an integer.')
          print (GV.dash_line)
          second_user_input = None
    elif user_input == 'MT':
      print (GV.dash_line)
      origin_type = user_input
      prompt_message = 'Please choose the type of molecule for the origin (the centre of mass of the first molecule of that type will be the origin)'
      origin_mol_type = chooseMolType(system,prompt_message) 
      origin_type_var = origin_mol_type
      for i,m in enumerate(system.mol_list):
        if m.label == origin_type_var:
          print ('')
          print ('I'*80)
          print ('III' + ('  The centre of mass of molecule ' + str(i+1) + ' (' + origin_type_var + ') will be the origin.  ').center(74) + 'III') 
          print ('I'*80)
          print ('')
          break
      with open(log_file,'a') as f:
        f.write('User chose origin to be the centre of mass of the first molecule of type ' + origin_type_var + '.\n')
        f.write('The centre of mass of molecule ' + str(i+1) + ' (' + origin_type_var + ') will be the origin.\n')
    elif user_input == 'C':
      second_user_input = None
      while not second_user_input:
        second_user_input = raw_input('Please enter the x,y and z co-ordinates of the origin separated by space: ').strip()
        coords = second_user_input.split()
        axis_labels = ['x','y','z']
        if len(coords) == 3:
          origin_type_var = []
          adjust_flag = False
          orig_coords = copy.deepcopy(coords)
          try:
            for i,c in enumerate(coords):
              coords[i] = float(c)
              orig_coords[i] = float(c)
              while (coords[i] > unit_cell_length/2.0):
                coords[i] -= unit_cell_length
                adjust_flag = True
              while (coords[i] < -unit_cell_length/2.0):
                coords[i] += unit_cell_length
                adjust_flag = True
              origin_type_var.append(coords[i])
            origin_type = user_input
            with open(log_file,'a') as f:
              f.write('User chose origin to be the co-ordinates ' + coords_to_string(orig_coords) + '.\n')
            if adjust_flag == True:
              print ('The input co-ordinates ' + coords_to_string(orig_coords) + ' is out of the unit cell. Origin adjusted to ' + coords_to_string(coords) + '.')
              with open(log_file,'a') as f:
                f.write('The input co-ordinates ' + coords_to_string(orig_coords) + ' is out of the unit cell. Origin adjusted to ' + coords_to_string(coords) + '.\n')
            else:
              print ('')
              print ('I'*80)
              print ('III' + ('  The origin is ' + coords_to_string(coords) + '  ').center(74) + 'III')
              print ('I'*80)
              print ('')
          except:
            print ('The ' +  axis_labels[i] + '-coordinate ' + str(c) + ' is invalid.')
            second_user_input = None
            print (GV.dash_line)
        else: 
          print ('You did not enter 3 co-ordinates.')
          print (GV.dash_line)
          second_user_input = None
  return (origin_type,origin_type_var)

#-------------------------------------------------------------------------------
# Function chooseConfigurations
#
# Operation: User interface to choose configurations from the trajectory to
#            create system objects out of which cluster are cut from.
#-------------------------------------------------------------------------------
def chooseConfigurations(n_configs,prompt_message=None,log_file=None):
  config_no_list = []
  if not prompt_message:
    prompt_message = 'Would you like to use:' 
  user_input = None
  while not user_input:
    print(textwrap.fill('A total of ' + str(n_configs) + ' configuration(s) found in HISTORY. ' + prompt_message,80))
    print('-->Every [n] configurations (starting with the 1st).')
    print('-->A [c]ustom set of configurations.') 
    user_input = raw_input('Please choose option [n/c]: ').strip()
    if user_input in ('C','c'):
      #regular_configs_flag = False
      user_input2 = None
      while not user_input2:
        user_input2 = raw_input(textwrap.fill('Please select configurations by entering a set of integers (between 1 and ' + str(n_configs) + ') separated by space (e.g. "1 2 4 15"):',80) + ' ').strip()
        user_config_nos = user_input2.split()
        config_no_list = []
        for i,n in enumerate(user_config_nos):
          if n.isdigit() == False:
            print (n + ' is not an integer')
            print (GV.dash_line)
            user_input2 = None
            break
          elif int(n) < 1 or int(n) > n_configs:
            print (n + ' is not between 1 and ' + str(n_configs))
            print (GV.dash_line)
            user_input2 = None
          else:
            if int(n) in config_no_list:
              print ('Config. no. ' + str(n) + ' already chosen, removing duplicate.')
            else:
              config_no_list.append(int(n))
              
        if user_input2:
          config_no_list = sorted(config_no_list)
          print (GV.dash_line)
          print ('I'*80)
          print ('III  ' + ('You have selected ' + str(len(config_no_list)) + ' configuration(s).').center(70) + '  III')
          print ('I'*80)
          with open (log_file,'a') as f:
            config_nos_str = ''
            for c in config_no_list:
              config_nos_str += str(c) + ' '
            f.write('Configuration(s) chosen manually. (Configuration no(s) ' + config_nos_str + 'included - ' + str(len(config_no_list)) + ' in total.)\n') 
    elif user_input in ('N','n'):
      user_input2 = None
      #regular_configs_flag = True
      while not user_input2:
        user_input2 = raw_input('Please enter n (an integer): ').strip()
        try:
          read_every_n_configs = int(user_input2.strip())
          if read_every_n_configs > n_configs:
            final_config_nos = [1]
            print ('')
            print ('I'*80)
            output_lines = textwrap.wrap(str(read_every_n_configs) + ' is larger than the total number of configurations. Only the first configuration in HISTORY will be used.',70)
            for l in output_lines:
              print ('III  ' + l.center(70) + '  III')
            print ('I'*80)
            print ('')
            if log_file:
              with open(log_file,'a') as f:
                f.write('Only the first configuration will be used.\n')
              read_every_n_configs = 0
          else:
            c = 1
            while c <= n_configs:
              config_no_list.append(c)
              c += read_every_n_configs
            #if n_configs % read_every_n_configs != 0:
            #  configs_to_extract = n_configs/read_every_n_configs + 1
            #else:
            #  configs_to_extract = n_configs/read_every_n_configs
            print ('')
            print ('I'*80)
            print ('III  ' + (str(len(config_no_list)) + ' configuration(s) will be extracted.').center(70) + '  III')
            print ('I'*80)
            print ('')
            if log_file:
              with open(log_file,'a') as f:
                f.write('Configurations every ' + str(read_every_n_configs) + ' snapshots will be used, giving ' + str(len(config_no_list)) + ' configurations in total.\n')
          break
        except:
          print (str(user_input2) + ' is not an integer.') 
          print (GV.dash_line)
          user_input2 = None
    else:
      user_input = None
      print ('Please enter \'n\' or \'c\'.')
      print (GV.dash_line)
  print (GV.dash_line)
  #print ('Loading configuration(s)...')
  #with open(log_file,'a') as f:
  #  f.write('Loading configuration(s)...\n')
  #configuration_list,unit_cell_length = ReadInput.readHistory(history_file,mol_info,n_configs,-1,config_nos)
  return (config_no_list)

#-------------------------------------------------------------------------------
# Function chooseCustomConfigurations
#
# Operation: Prompts user to select a set of configurations by configuration 
#            number within the trajectory file (to generate systems for cluster 
#            cutting)
#-------------------------------------------------------------------------------
def chooseCustomConfigurations(n_configs):
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('Please select configurations by entering a set of integers (between 1 and ' + str(n_configs) + ') separated by space (e.g. "1 2 4 15"):',80) + ' ').strip()
    user_config_nos = user_input.split()
    final_config_nos = [] 
    for i,n in enumerate(user_config_nos):
      if n.isdigit() == False:
        print (n + ' is not an integer')
        print (GV.dash_line)
        user_input = None
        break
      elif int(n) < 1 or int(n) > n_configs:
        print (n + ' is not between 1 and ' + str(n_configs))
        print (GV.dash_line)
        user_input = None
      else:
        if int(n) in final_config_nos:
          print ('Config. no. ' + str(n) + ' already chosen, removing duplicate.')
        else:
          final_config_nos.append(int(n))
    if user_input:
      print (GV.dash_line)
      print ('I'*80)
      print ('III  ' + ('You have selected ' + str(len(final_config_nos)) + ' configuration(s).').center(70) + '  III')
      print ('I'*80)
  return (final_config_nos)


#-------------------------------------------------------------------------------
# Function chooseParticleNo
#
# Operation: Prompts user to choose a particle by its particle number
#-------------------------------------------------------------------------------
def chooseParticleNo(n_particles,particle_name):
  user_input = None
  while not user_input:
    user_input = raw_input('Please enter ' + particle_name + ' number [1-' + str(n_particles) + ']: ').strip()
    if user_input.isdigit() == False:
      print (user_input + ' is not an integer.')
      user_input = None
    else:
      if int(user_input) < 1 or int(user_input) > n_particles:
        print (user_input + ' is not between 1 and ' + str(n_particles) + '.')
        user_input = None
      else:
        particle_no = int(user_input)
  return (particle_no)

#-------------------------------------------------------------------------------
# Function chooseMolNo
#-------------------------------------------------------------------------------
def chooseMolNo(n_mols,mol_name=None):
  #user_input = None
  #while not user_input:
  #  user_input = raw_input('Please enter molecule number [1-' + str(n_mols) + ']: ').strip()
  #  if user_input.isdigit() == False:
  #    print (user_input + ' is not an integer.')
  #    user_input = None
  #  else:
  #    if int(user_input) < 1 or int(user_input) > n_mols:
  #      print (user_input + ' is not between 1 and ' + str(n_mols) + '.')
  #      user_input = None
  #    else:
  #      mol_no = int(user_input)
  mol_no = chooseParticleNo(n_mols,'molecule')
  return (mol_no)

#-------------------------------------------------------------------------------
# Function chooseAtomNo
#-------------------------------------------------------------------------------
def chooseAtomNo(n_atoms,atom_name=None):
  #user_input = None
  #while not user_input:
  #  user_input = raw_input('Please enter atom number [1-' + str(n_atoms) + ']: ').strip()
  #  if user_input.isdigit() == False:
  #    print (user_input + ' is not an integer.')
  #    user_input = None
  #  else:
  #    if int(user_input) < 1 or int(user_input) > n_atoms:
  #      print (user_input + ' is not between 1 and ' + str(n_atoms) + '.')
  #      user_input = None
  #    else:
  #      atom_no = int(user_input)
  atom_no = chooseParticleNo(n_atoms,'atom')
  return (atom_no)

#-------------------------------------------------------------------------------
# Function chooseMolType
#
# Operation: Prompts user to choose a molecule type
#-------------------------------------------------------------------------------
def chooseMolType(system,prompt_message=None):
  if not prompt_message:
    prompt_message = 'Please choose molecule type'
  n_mol_types = len(system.mol_info)
  user_input = None
  while not user_input:
    print ('Molecule types in the system')
    for i,m in enumerate(system.mol_info):
      print ('[' + str(i+1) + '] ' + m['molname'])
    user_input = raw_input(textwrap.fill(prompt_message + ' [1-' + str(n_mol_types) + ']:',80) + ' ').strip()
    if user_input.isdigit():
      user_input = int(user_input)
      if user_input < 1 or user_input > len(system.mol_info):
        print ('You did not enter an integer between 1 and ' + str(n_mol_types) + '.')
        print (GV.dash_line)
        user_input = None
      else:
        mol_type = system.mol_info[user_input-1]['molname']
    else:
      print (user_input + ' is not an integer')
      user_input = None
      print (GV.dash_line)
  return (mol_type)

#-------------------------------------------------------------------------------
# Function chooseAtomType
#
# Operation: Prompts user to choose an atom type.
#-------------------------------------------------------------------------------
def chooseAtomType(system,prompt_message=None):
  if not prompt_message:
    prompt_message = 'Please choose atom type'
  n_atom_types = len(system.atom_info.keys()) 
  user_input = None
  #atoms_names_and_mols = []
  atom_names = sorted(system.atom_info.keys())
  while not user_input:
    print ('Atom types in the system:')
    for i,a in enumerate(atom_names):
      print ('-->[' + str(i+1) + '] ' + a)
    user_input = raw_input(textwrap.fill(prompt_message + ' [1-' + str(n_atom_types) + ']:',80) + ' ').strip()
    if user_input.isdigit():
      user_input = int(user_input)
      if user_input < 1 or user_input > len(system.atom_info.keys()):
        print ('You did not enter an integer between 1 and ' + str(n_atom_types) + '.')
        print (GV.dash_line)
        user_input = None
      else:
        atom_type = atom_names[user_input-1]
    else:
      print (str(user_input) + ' is not an integer.')
      print (GV.dash_line)
      user_input = None
  return (atom_type)


#-------------------------------------------------------------------------------
# Function chooseSetOfAtomTypes
#
# Operation: Prompts user to select a set of atom types.
#            Used in the calculation of the bond correlation function.
#-------------------------------------------------------------------------------
def chooseSetOfAtomTypes(system,prompt_message=None):
  atom_types = []
  guide_message = 'To choose more than one type enter the label numbers separated by space (e.g. 1 3 4):'
  if not prompt_message: 
    prompt_message = 'Please choose atom type(s).'
  prompt_message = prompt_message + ' ' + guide_message
  n_atom_types = len(system.atom_info.keys()) 
  user_input = None
  while not user_input:
    print ('Atom types in the system:')
    for i,a in enumerate(sorted(system.atom_info.keys())):
      print ('-->[' + str(i+1) + ']' + a)
    user_input = raw_input(textwrap.fill(prompt_message + ' [1-' + str(n_atom_types) + ']:',80) + ' ').strip()
    user_input_list = user_input.split()
    for n in user_input_list:
      if n.isdigit():
        n = int(n)
        if sorted(system.atom_info.keys())[n-1] in atom_types:
          print ('More than one instance of ' + str(n) + ', removing duplicate.')
        elif n < 1 or n > len(system.atom_info.keys()):
          print (str(n) + ' is not an integer between 1 and ' + str(n_atom_types) + '.')
          print (GV.dash_line)
          user_input = None
          break
        else:
          atom_type = sorted(system.atom_info.keys())[n-1]
          atom_types.append(atom_type)
      else:
        print (str(n) + ' is not an integer.')
        print (GV.dash_line)
        user_input = None
        break
  return (atom_types)


#-------------------------------------------------------------------------------
# Function chooseAtomNoSet
#
# Operation: Prompts user to choose a set of atoms by atom number.
#-------------------------------------------------------------------------------
def chooseAtomNoSet(system):
  user_input = None
  n_atoms = len(system.atom_list)
  while not user_input:
    user_input = raw_input(textwrap.fill('Please select set of atom numbers (integers between 1 and ' + str(n_atoms) + ' separated by space (e.g. "1 5 6 13"):',80) + ' ').strip()
    atom_nos = user_input.split()
    final_atom_nos = []
    for i,n in enumerate(atom_nos):
      if n.isdigit() == False:
        print (n + ' is not an integer. Please try again.')
        user_input = None
        #break
      elif int(n) < 1 or int(n) > n_atoms:
        print (n + ' is not between 1 and ' + str(n_atoms) + '. Please try again.')
        user_input = None
      else:
        if int(n) in final_atom_nos:
          print ('Atom no. ' + str(n) + ' already chosen, removing dupicate.')
        else:
          final_atom_nos.append(int(n))

    if user_input:
      print ('You have selected the following atoms:')
      for a in final_atom_nos:
        print ('-->No. ' + str(a) + ' (' + system.atom_list[a-1].label + ')')
      user_confirmation = None
      while not user_confirmation: 
        user_confirmation = raw_input('Do you wish to continue [Y/N]? ').strip().upper()
        if user_confirmation == 'Y':
          pass
        elif user_confirmation == 'N':
          user_input = None 
        else:
          print ('Please enter \'Y\' or \'N\'.')
          user_confirmation = None
  return (final_atom_nos)

#-------------------------------------------------------------------------------
# Function chooseMolNoSet
# 
# Operation: Prompts user to choose a set of molecules by molcule number
#-------------------------------------------------------------------------------
def chooseMolNoSet(system):
  user_input = None
  n_mol = len(system.mol_list)
  while not user_input:
    user_input = raw_input(textwrap.fill('Please select set of molecule numbers (integers between 1 and ' + str(n_mol) + ' separated by space (e.g. "1 5 6 13"):',80) + ' ').strip()
    mol_nos = user_input.split()
    final_mol_nos = []
    for i,n in enumerate(mol_nos):
      if n.isdigit() == False:
        print (n + ' is not an integer. Please try again.')
        user_input = None
        #break
      elif int(n) < 1 or int(n) > n_mol:
        print (n + ' is not between 1 and ' + str(n_mol) + '. Please try again.')
        user_input = None
      else:
        if int(n) in final_mol_nos:
          print ('Atom no. ' + str(n) + ' already chosen, removing dupicate.')
        else:
          final_mol_nos.append(int(n))

    if user_input:
      print ('You have selected the following molecules:')
      for a in final_mol_nos:
        print ('-->No. ' + str(a) + ' (' + system.mol_list[a-1].label + ')')
      user_confirmation = None
      while not user_confirmation: 
        user_confirmation = raw_input('Do you wish to continue [Y/N]? ').strip().upper()
        if user_confirmation == 'Y':
          pass
        elif user_confirmation == 'N':
          user_input = None 
        else:
          print ('Please enter \'Y\' or \'N\'.')
          user_confirmation = None
  return (final_mol_nos)

#-------------------------------------------------------------------------------
# Function chooseOriginAtomNo
#
# Operation: Prompts user to choose an origin by atom number.
#-------------------------------------------------------------------------------
def chooseOriginAtomNo(n_atoms):
  ######################
  # Origin atom number #
  ######################
  user_input = None
  while not user_input:
    user_input = raw_input('Please choose origin atom number [1-' + str(n_atoms) + ']: ').strip()
    if user_input.isdigit() == False:
      print (user_input + ' is not an integer.')
      user_input = None
    else:
      if int(user_input) < 1 or int(user_input) > n_atoms:
        print (user_input + ' is not between 1 and ' + str(n_atoms) + '.')
        user_input = None
      else:
        origin_atom_no = int(user_input)
  return (origin_atom_no)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
#def retDecimalString(number,decimal_places):
#  return
