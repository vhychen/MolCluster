# read parameter input file

from decimal import Decimal
import MolClass as MC
import math
import os
import os.path as op
import pickle
import shutil

kJ_to_kcal = 0.2390057361377
eV_to_kcal = 23.06054194533
K_to_kcal = 0.001987191686486
internal_to_kcal = 0.002390057361377


TEMP_DIR = 'temp'
TEMP_SYSTEM_DIR=op.join('temp','SYS_OBJECTS')
TEMP_SYSTEM_BASE=op.join('temp','SYS_OBJECTS','config')
TEMP_SPLIT_CONFIG_DIR=op.join('temp','HISTORY_PARTS')
TEMP_SPLIT_CONFIG_BASE=op.join('temp','HISTORY_PARTS','config')

#-------------------------------------------------------------------------------
# Function tail
#-------------------------------------------------------------------------------
def tail(filename,count,offset=1024):
  """
  A more efficent way of getting the last few lines of a file.
  Depending on the length of your lines, you will want to modify offset
  to get better performance.
  """
  f_size = os.stat(filename).st_size
  if f_size == 0:
    return []
  with open(filename, 'r') as f:
    if f_size <= offset:
      offset = int(f_size / 2)
    while True:
      seek_to = min(f_size - offset, 0)
      f.seek(seek_to)
      lines = f.readlines()
      # Empty file
      if seek_to <= 0 and len(lines) == 0:
        return []
      # count is larger than lines in file
      if seek_to == 0 and len(lines) < count:
        return lines
      # Standard case
      if len(lines) >= (count + 1):
        return lines[count*-1:]

#-------------------------------------------------------------------------------
# Function genTempDirs
#
# Operation: Generates directories to store temporary files
#-------------------------------------------------------------------------------
def genTempDirs():
  if not op.isdir(TEMP_DIR):
    os.makedirs(TEMP_DIR)
  if not op.isdir(TEMP_SYSTEM_DIR):
    os.makedirs(TEMP_SYSTEM_DIR)
  return

#-------------------------------------------------------------------------------
# Function genTempDirs
#
# Operation: Deletes temporary directories and their contents
#-------------------------------------------------------------------------------
def clearTempDirs():
  shutil.rmtree(TEMP_SYSTEM_DIR)
  shutil.rmtree(TEMP_SPLIT_CONFIG_DIR)
  return

#-------------------------------------------------------------------------------
# Function unitsConverter
#
# IMPORTANT: Assume lengths in units of Angstroms
#
#-------------------------------------------------------------------------------
def unitsConverter(orig_units,quantity):
  if orig_units.lower() == 'kj':
    return (quantity*kJ_to_kcal)
  elif orig_units.lower() == 'ev':
    return (quantity*eV_to_kcal)
  elif orig_units.lower() == 'k':
    return (quantity*K_to_kcal)
  elif orig_units.lower() == 'internal':
    return (quantity*internal_to_kcal)
  elif orig_units.lower() == 'kcal':
    return (quantity)
  else:
    print ('Defaulting to DL_POLY internal units (10 J/mol)')
    return (quantity*internal_to_kcal)

#-------------------------------------------------------------------------------
# Function saveSystemObj
#
# Operation: Save the system object file so that it can be loaded another time
#-------------------------------------------------------------------------------
def saveSystemObj(sys_obj,filename):
  with open(filename,'wb') as f:
    pickle.dump(sys_obj,f,pickle.HIGHEST_PROTOCOL)
  return

#-------------------------------------------------------------------------------
# Function loadSystemObj
#
# Operation: Load a saved system object file
#-------------------------------------------------------------------------------
def loadSystemObj(sys_obj_file):
  with open(sys_obj_file) as f:
    sys_obj = pickle.load(f)
  return (sys_obj)
  
#-------------------------------------------------------------------------------
# Function field
#
# Operation: Read a DL_POLY Field file
#-------------------------------------------------------------------------------
def field(fieldfile='example_inputs/FIELD'):
  molecule_array = []
  forcefield_array = []
  with open(fieldfile) as f:
    firstline = f.readline()
    unitsline = f.readline()
    units = unitsline.split()[-1]
    typesline = f.readline()
    types_of_molecules = int(typesline.split()[-1])

    # Read the composition and force field parameters for each molecule type
    for n in range(types_of_molecules):
      molecule_array.append({})
      molecule_name = f.readline().strip()
      molecule_array[-1]["molname"] = molecule_name
      molecule_array[-1]["nmols"] = int(f.readline().split()[-1])
      molecule_array[-1]["natoms"] = int(f.readline().split()[-1])
      atomlist = []
      masslist = []
      chargelist = []
      atominfo = {}
      for a in range(molecule_array[-1]["natoms"]):
        atomline = f.readline().split()
        atom_label = atomline[0]
        mass = float(atomline[1])
        charge = float(atomline[2])
        atomlist.append(atom_label) 
        if atom_label not in atominfo:
          atominfo[atom_label] = {'natoms':1,\
                                   'mass':mass,\
                                   'charge':charge}
        else:
          atominfo[atom_label]['natoms'] += 1
        masslist.append(mass)
        chargelist.append(charge)
      molecule_array[-1]["atomlist"] = atomlist
      molecule_array[-1]["masslist"] = masslist
      molecule_array[-1]["chargelist"] = chargelist
      molecule_array[-1]["atominfo"] = atominfo
      
      nextline = f.readline().split()  
      molecule_array[-1]['potential'] = {}
      molecule_array[-1]['connectivity'] = []
      while nextline[0] != "finish":
        n_entries = int(nextline[1])
        if nextline[0] == 'bonds':
          molecule_array[-1]['potential']['bonds'] = []
          for i in range(n_entries):
            entryline = f.readline().split()
            if entryline[0] == 'harm':
              atom1 = int(entryline[1])
              atom2 = int(entryline[2])
              k = unitsConverter(units,float(entryline[3]))
              r0 = float(entryline[4])
              molecule_array[-1]['potential']['bonds'].append({'type':entryline[0],'atom1':atom1,'atom2':atom2,'k':k,'r0':r0})
              molecule_array[-1]['connectivity'].append((atom1,atom2))
            else:
              print ('Bond type not supported')
        elif nextline[0] == 'angles':
          molecule_array[-1]['potential']['angles'] = []
          for i in range(n_entries):
            entryline = f.readline().split()
            if entryline[0] == 'harm':
              atom1 = int(entryline[1])
              atom2 = int(entryline[2])
              atom3 = int(entryline[3])
              k = unitsConverter(units,float(entryline[4]))
              theta0 = float(entryline[5])
              molecule_array[-1]['potential']['angles'].append({'type':entryline[0],'atom1':atom1,'atom2':atom2,'atom3':atom3,'k':k,'theta0':theta0})
            else:
              print ('Angle type not supported')
        elif nextline[0] == 'dihedrals':
          molecule_array[-1]['potential']['dihedrals'] = []
          for i in range(n_entries):
            entryline = f.readline().split()
            if entryline[0] == 'cos3':
              atom1 = int(entryline[1])
              atom2 = int(entryline[2])
              atom3 = int(entryline[3])
              atom4 = int(entryline[4])
              A1 = unitsConverter(units,float(entryline[5]))
              A2 = unitsConverter(units,float(entryline[6]))
              A3 = unitsConverter(units,float(entryline[7]))
              molecule_array[-1]['potential']['dihedrals'].append({'type':entryline[0],'atom1':atom1,'atom2':atom2,'atom3':atom3,'atom4':atom4,'A1':A1,'A2':A2,'A3':A3})
            elif entryline[0] == 'cos':
              atom1 = int(entryline[1])
              atom2 = int(entryline[2])
              atom3 = int(entryline[3])
              atom4 = int(entryline[4])
              A = unitsConverter(units,float(entryline[5]))
              delta = float(entryline[6])
              m = float(entryline[7])
              molecule_array[-1]['potential']['dihedrals'].append({'type':entryline[0],'atom1':atom1,'atom2':atom2,'atom3':atom3,'atom4':atom4,'A':A,'delta':delta,'m':m})
            else:
              print ('Dihedral type not supported')
        elif nextline[0] == 'rigid':
          n_entries = int(nextline[1])
          for i in range(n_entries):
            entryline = f.readline()
            #print (i)
            #print (entryline)
            entryline = entryline.split()
            n_atoms = int(entryline[0])
            rigid_list = entryline[1:]
            # VC: Insert user-defined connectivity  
            molecule_array[-1]['rigid'] = rigid_list
            if 'constraints' not in molecule_array[-1]['potential']:
              molecule_array[-1]['potential']['constraints'] = []
            #print ('Molcule ' + molecule_name + ' has a rigid component')
            for i in range(len(rigid_list)):
              for j in range(i+1,len(rigid_list)):
                molecule_array[-1]['potential']['constraints'].append({'atom_1':rigid_list[i],'atom_2':rigid_list[j],'dist':'NA'})
                #print ('--> constraint defined between ' + rigid_list[i] + ' and ' + rigid_list[j])
        elif nextline[0] == 'constraints':
          molecule_array[-1]['potential']['constraints'] = []
          n_entries = int(nextline[1])
          for i in range(n_entries):
            entryline = f.readline().split()
            atom_1,atom_2,dist = entryline
            molecule_array[-1]['potential']['constraints'].append({'atom_1':atom_1,'atom_2':atom_2,'dist':dist})
            molecule_array[-1]['connectivity'].append((int(atom_1),int(atom_2)))
        else:
          print ('Interaction type ' + entryline[0] + ' not supported')
          print (nextline)
          for i in range(n_entries):
            entryline = f.readline()

        nextline = f.readline().split()
        #print (nextline)
    nextline = f.readline().split()
    vdw_array = [] 
    while nextline[0] != 'close':
      n_entries = int(nextline[1])
      if nextline[0].lower() == 'vdw':
        for i in range(n_entries):
          entryline = f.readline().split()
          if entryline[2] == 'lj':
            atom_type1 = entryline[0] 
            atom_type2 = entryline[1]
            epsilon = unitsConverter(units,float(entryline[3]))
            sigma = float(entryline[4])
            c6 = 4*epsilon*math.pow(sigma,6)
            c12 = 4*epsilon*math.pow(sigma,12)
            vdw_array.append({'atom_type1':atom_type1,'atom_type2':atom_type2,'c6':c6,'c12':c12})
          elif entryline[2] == '12-6':
            atom_type1 = entryline[0] 
            atom_type2 = entryline[1]
            A = unitsConverter(units,float(entryline[3]))
            B = unitsConverter(units,float(entryline[4]))
            c6 = B
            c12 = A
            vdw_array.append({'atom_type1':atom_type1,'atom_type2':atom_type2,'c6':c6,'c12':c12})
          else:
            print ('van der Waals term not supported')
      else:
        print ('Non-bonded interaction type ' + nextline[0] + ' not supported')
      nextline = f.readline().split()
  n_atoms = 0
  for t in molecule_array:
    n_atoms += t['nmols']*t['natoms']
  return(molecule_array,vdw_array,n_atoms)     

def gen_chemshell_ff_file(molecule_array,vdw_array,ff_file='ff.dat'):
  atom_list = []
  with open(ff_file,'w') as f:
    for m in molecule_array:
      for a in m['atominfo']:
        if a not in atom_list:
          atom_list.append(a)
          f.write('declare ' + a + '\n')    
    for vdw_term in vdw_array:
      f.write('vdw ' + vdw_term['atom_type1'] + ' ' + vdw_term['atom_type2'] + ' ' + str(vdw_term['c6']) + ' ' + str(vdw_term['c12']) + '\n')
    for m in molecule_array:
      #print (m['potential'].keys())
      if 'bonds' in m['potential']:
        for bond_term in m['potential']['bonds']:
          f.write('bond ' + m['atomlist'][bond_term['atom1']-1] + ' ' + m['atomlist'][bond_term['atom2']-1] + ' ' + str(bond_term['k']) + ' ' + str(bond_term['r0']) + '\n')
      if 'angles' in m['potential']:
        #print ('angles found')
        for angle_term in m['potential']['angles']:
          f.write('angle ' + m['atomlist'][angle_term['atom1']-1] + ' ' + m['atomlist'][angle_term['atom2']-1] + ' ' + m['atomlist'][angle_term['atom3']-1] + ' ' + str(angle_term['k']) + ' ' + str(angle_term['theta0']) + '\n')
      if 'dihedrals' in m['potential']:
        for dihedral_term in m['potential']['dihedrals']:
          if dihedral_term['type'] == 'cos3':
            f.write('mm2ptor ' + m['atomlist'][dihedral_term['atom1']-1] + ' ' + m['atomlist'][dihedral_term['atom2']-1] + ' ' +  m['atomlist'][dihedral_term['atom3']-1] + ' ' + m['atomlist'][dihedral_term['atom4']-1] + ' ' + str(dihedral_term['A1']) + ' ' + str(dihedral_term['A2']) + ' ' + str(dihedral_term['A3']) + '\n')
          elif dihedral_term['type'] == 'cos':
            f.write('ptor ' + m['atomlist'][dihedral_term['atom1']-1] + ' ' + m['atomlist'][dihedral_term['atom2']-1] + ' ' +  m['atomlist'][dihedral_term['atom3']-1] + ' ' + m['atomlist'][dihedral_term['atom4']-1] + ' ' + str(dihedral_term['k']) + ' ' + str(dihedral_term['delta']) + ' ' + str(dihedral_term['m']) + 'i-j-k-l' +  '\n')
          else:
            print ('dihedral_term not recognised')
  return

#-------------------------------------------------------------------------------
# Function configurationsInHistory
#
#
#-------------------------------------------------------------------------------
def configurationsInHistory(history_file,n_atoms):
  single_config_lines = n_atoms*2+4
  with open(history_file) as f:
    title_line = f.readline()
    history_n_atoms = int(f.readline().split()[-1])
    if history_n_atoms != n_atoms:
      return (-1)
    timestep_line = f.readline().split()
    first_timestep = int(timestep_line[1])
    print ('first_timestep = ' + str(first_timestep))
    #timestep_size_dp = len(timestep_line[-1].split('.')[1])
    timestep_size = Decimal(timestep_line[-1])
    for i in range(3):
      f.readline()
    for i in range(n_atoms):
      f.readline()
      f.readline()
  
    second_timestep_line = f.readline().split()
    second_timestep = int(second_timestep_line[1])
    print ('second_timestep = ' + str(second_timestep))
  final_timestep_line = tail(history_file,single_config_lines)[0].split()
  if final_timestep_line[0] != 'timestep':
    return (-1)
  final_timestep = int(final_timestep_line[1])
  print ('final_timestep = ' + str(final_timestep))
  ns_between_configs = (second_timestep - first_timestep)*timestep_size/1000
  n_configs = (final_timestep - first_timestep)/(second_timestep - first_timestep) + 1
  return (n_configs,timestep_size,ns_between_configs)

#-------------------------------------------------------------------------------
# Function readSingleSnapshot
#
#-------------------------------------------------------------------------------
def readSingleSnapshot(f,mol_array,read_header=True,return_timestep_size=False):
  coordinates_array = []
  if read_header == True:
    timestep_line = f.readline()
    timestep = int(timestep_line.split()[1])
    if return_timestep_size == True:
      timestep_size = float(timestep_line.split()[-1])
    unit_cell_length = float(f.readline().split()[0])
    f.readline()
    f.readline()
    #f.readline()
  else:
    timestep = 1
    unit_cell_length = None
  for m in mol_array:
    for i in range(m['nmols']):
      coordinates_array.append([])
      for j in range (m["natoms"]):
        f.readline()
        coordinates = f.readline().split()
        for k in range(3):
          coordinates[k] = float(coordinates[k])
        coordinates_array[-1].append(coordinates)
  print ('Loaded configuration from timestep: ' + str(timestep))
  if return_timestep_size == True:
    return (coordinates_array,timestep,unit_cell_length,timestep_size)
  else:
    return (coordinates_array,timestep,unit_cell_length)

#def splitHistory(history_file,config_file_base=op.join('temp','HISTORY_PARTS','config')):
#  return

#-------------------------------------------------------------------------------
# Function genSystemsFromHistory
#
# Operation:Read configuration(s) and generate System objects from HISTORY file
#-------------------------------------------------------------------------------
def genSystemsFromHistory(history_file,mol_info,vdw_info,config_nos,n_atoms,timestep_size,temp_system_base=TEMP_SYSTEM_BASE):
  config_block_lines = 3*n_atoms+4
  with open(history_file) as f:
    firstline = f.readline()
    secondline = f.readline()
    # Read custom list of configs
    config = 1  
    while len(config_nos) > 0:
      if config in config_nos:
        output_file = temp_system_base + str(config) + '.system'
        coordinates,timestep,unit_cell_length = readSingleSnapshot(f,mol_info)
        print ('Generating BulkSystem object for configuration ' + str(config) + '; timestep = ' + str(timestep) + '; time = ' + str(timestep*timestep_size-timestep_size))
        system_obj = MC.BulkSystem(coordinates,mol_info,unit_cell_length,vdw_info,timestep,timestep*timestep_size-timestep_size)
        saveSystemObj(system_obj,output_file)
        #snapshots_array.append(snapshot)
        #timestep_array.append(timestep)
        config_nos.remove(config)
      else:
        for j in range(config_block_lines):
          f.readline()
      config += 1
  return 

def genConfigurationInfoForListOfConfigs(history_file,mol_info,config_no_lit,n_atoms,timestep_size,log_file):
  print ('Loading configuration(s)...')
  if log_file:
    with open(log_file,'a') as f:
      f.write('Loading configuration(s)...\n')
  configuration_list,timestep_list,unit_cell_length_list = ReadInput.genCoordinatesListFromHistory(history_file,mol_info,config_no_list,n_atoms,timestep_size)
  print ('...complete.')
  return (configuration_list,timestep_list,unit_cell_length_list)

#-------------------------------------------------------------------------------
# Function genCoordinatesListFromHistory
#
# Operation: Read configuration(s) and generate coordinates list from HISTORY 
#            file
#-------------------------------------------------------------------------------
def genCoordinatesListFromHistory(history_file,mol_info,config_nos,n_atoms,timestep_size,temp_system_base=TEMP_SYSTEM_BASE):
  config_block_lines = 2*n_atoms+4
  coordinates_list = []
  timestep_list = []
  unit_cell_length_list = []
  with open(history_file) as f:
    firstline = f.readline()
    secondline = f.readline()
    # Read custom list of configs
    config = 1  
    while len(config_nos) > 0:
      if config in config_nos:
        #output_file = temp_system_base + str(config) + '.system'
        coordinates,timestep,unit_cell_length = readSingleSnapshot(f,mol_info)
        print ('Generated coordinates list for configuration ' + str(config) + '; timestep = ' + str(timestep) + '; time = ' + str((timestep-1)*timestep_size))
        coordinates_list.append(coordinates)
        timestep_list.append(timestep)
        unit_cell_length_list.append(unit_cell_length)
        #system_obj = MC.BulkSystem(coordinates,mol_info,unit_cell_length,vdw_info,timestep,timestep*timestep_size-timestep_size)
        #saveSystemObj(system_obj,output_file)
        #snapshots_array.append(snapshot)
        #timestep_array.append(timestep)
        config_nos.remove(config)
      else:
        for j in range(config_block_lines):
          f.readline()
      config += 1
  return (coordinates_list,timestep_list,unit_cell_length_list)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
def readxyz(xyzfile, molecule_array):
  #molecule_array = field(fieldfile)
  with open(xyzfile) as f:
    coordinates_array = []
    while 1:
      atomnumberline = f.readline()
      if not atomnumberline: break
      else:
        number_of_atoms = int(atomnumberline.split()[-1])
      commentline = f.readline()
      for molecule in molecule_array: 
        for i in range (molecule["nmols"]):
          coordinates_array.append([])
          for j in range (molecule["natoms"]):
            coordinates = f.readline().split()[1:4]
            coordinates_array[-1].append(coordinates)
#  print coordinates_array     
  return (coordinates_array)

#########################
### DEFUNCT FUNCTIONS ###
#########################

#def buildSystemsFromHISTORY(field_file='example_inputs/FIELD',history_file='example_inputs/HISTORY',read_every_n_configs=0):
#
#  # Read inputs, starting with FIELD file getting info about the molecules and potentials...
#  mol_info,vdw_info = ReadInput.field(field_file)
#
#  # ... then HISTORY file getting the configurational info.
#  snapshots_array,unit_cell_length,timestep_array = ReadInput.readHistory(history_file,mol_info,read_every_n_configs)
#
#  # Create a array of BulkSystem objects.
#  systems_array = []
#  for s in snapshots_array:
#    systems_array.append(MolClass.BulkSystem(s,mol_info,unit_cell_length,vdw_info))
#  return (systems_array)
#
#def outputClustersFromHISTORY(field_file='example_inputs/FIELD',history_file='example_inputs/HISTORY',read_every_n_configs=0,output_directory='clusters/',basename='cluster'):
#  bulk_systems = buildSystemsFromHISTORY(field_file,history_file,read_every_n_configs)
#  return

#-------------------------------------------------------------------------------
# Function readHistory
#
# Read configuration(s) from whole HISTORY file
# *** DEFUNCT - please use genCoordinatesListFromHistory ***
#-------------------------------------------------------------------------------
#def readHistory(historyfile,mol_info,tot_configs,read_every_n_configs=0,config_nos=None):
#  #molecule_array = field(fieldfile)
#  with open(historyfile) as f:
#    snapshots_array = []
#    unit_cell_length = None
#    firstline = f.readline()
#    secondline = f.readline()
#    number_of_atoms = int(secondline.split()[-1])
#    coordinates_array = []
#    timestepline = f.readline()
#    timestep = int(timestepline.split()[1])
#    timestep_size = float(timestepline.split()[-1])
#    unit_cell_length = float(f.readline().split()[0])
#    f.readline()
#    f.readline()
#    timestep_array = []
#    # Read first snapshot
#    if read_every_n_configs == 0:
#      snapshot,timestep,unit_cell_length_check = readSingleSnapshot(f,mol_info,False)
#      #print ('Loaded configuration from timestep: 1')
#      snapshots_array.append(snapshot)
#      timestep_array.append(timestep)
#
#    # Read every 'read_every_n_configs' configs
#    elif read_every_n_configs >= 1:
#      snapshot,timestep,unit_cell_length_check = readSingleSnapshot(f,mol_info,False)
#      #print ('Loaded configuration from timestep: 1')
#      snapshots_array.append(snapshot)
#      timestep_array.append(int(timestep))
#      for i in range((tot_configs-1)/read_every_n_configs):
#        for j in range(read_every_n_configs - 1):
#          for k in range(4 + number_of_atoms*2):
#            f.readline()
#        snapshot,timestep,unit_cell_length_check = readSingleSnapshot(f,mol_info)
#        snapshots_array.append(snapshot)
#        timestep_array.append(timestep)
#
#    # Read custom list of configs
#    elif read_every_n_configs == -1:
#      count = 1  
#      while len(config_nos) > 0:
#        if count in config_nos:
#          if count == 1: 
#            snapshot,timestep,unit_cell_length_check = readSingleSnapshot(f,mol_info,False)
#            #print ('Loaded configuration from timestep: 1')
#            snapshots_array.append(snapshot)
#            timestep_array.append(int(timestep))
#          else:
#            snapshot,timestep,unit_cell_length_check = readSingleSnapshot(f,mol_info)
#            snapshots_array.append(snapshot)
#            timestep_array.append(timestep)
#          config_nos.remove(count)
#        else:
#          if count > 1:
#            for j in range(4):
#              f.readline()
#          for j in range(number_of_atoms*2):
#            f.readline()
#        count += 1
#  return (snapshots_array,unit_cell_length,timestep_array)
