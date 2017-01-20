#!/usr/bin/python

#import math
import MolClass as MC
import ReadInput as RI
import CommonFunctions as CF
import GlobalVariables as GV
import sys
import os
import os.path as op
import textwrap
import ChemShellInputFunctions as CSIF
import shutil

bohr_in_ang = 1.889725989

#-------------------------------------------------------------------------------
# Function genFullEmbeddedClusters
#
# Operation: 1) Load bare clusters objects generated previously
#            2) Read point charge positions and charge magnitudes
#            3) Combine information to create embeded cluster objects
#            4) Generate ChemShell QM/MM optimisation file
#-------------------------------------------------------------------------------
def genFullEmbeddedClusters(workdir,embedded_cluster_specs_obj=None,log_file=None):
  if not embedded_cluster_specs_obj: 
    embedded_cluster_specs_obj = RI.loadSystemObj(op.join(workdir,'embedded_cluster_specs.obj'))

  n_clusters = embedded_cluster_specs_obj.n_clusters

  ##########################
  # Load first bulk system #
  ##########################
  first_system = RI.loadSystemObj(op.join(workdir,'embedded_cluster_1','bulk.obj'))
 
  #########################
  # Initial cluster specs #
  #########################
  cluster_construction_type = 'R'
  active_region_type = 'S'
  cluster_radius = embedded_cluster_specs_obj.r_cluster_in_ang
  active_radius = embedded_cluster_specs_obj.r_active_in_ang
  
  user_input = None
  while not user_input:
    #embedded_cluster_specs_obj.printSpecs()  
    print ('===Cluster specifications===')
    if cluster_construction_type == 'R':
      print ('-->Cluster: sphere with radius = ' + str(cluster_radius) + ' Angstroms')
    elif cluster_construction_type == 'M':
      print ('-->Cluster: fixed number of each molecule type')
      for i,m in enumerate(cluster_mol_numbers_list):
        print ('---->' + str(m) + ' X ' + first_system.mol_info[i]['molname'])
    if active_region_type == 'S':
      print ('-->Active region: sphere with radius of ' + str(active_radius) + ' Angstroms')
    elif active_region_type == 'M':
      print ('-->Active region: fixed number of total molecules = ' + str(active_region_n_mol))
    elif active_region_type == 'A':
      print ('-->Active region: fixed upper limit of total atoms = ' + str(active_region_n_atoms))
    elif active_region_type == 'N':
      print ('-->Active region: fixed number of each molecule type')
      for i,m in enumerate(active_mol_number_list):
        print ('---->' + str(m) + ' X ' + first_system.mol_info[i]['molname'])

    user_input = raw_input('Use these specs [Y/N]? ').strip().upper()
    if user_input == 'Y':
      with open(log_file,'a') as f:
        if cluster_construction_type == 'R':
          f.write('Cluster(s) of radius ' + str(cluster_radius) + ' ' + GV.ANGSTROM_STR + ' to be cut.\n')
        elif cluster_construction_type == 'M':
          f.write('Cluster(s) with specified numbers of each molecule type to be cut. Each cluster to be formed of:\n')
          for i,m in enumerate(cluster_mol_numbers_list):
            f.write('-->' + str(m) + ' X ' + first_system.mol_info[i]['molname'] + '\n')
           
        if active_region_type == 'S':
          if active_radius == 'A':
            f.write('User requested all atoms to be included in the active region.\n')
          else:
            f.write('User requested the active region to be defined as a spherical region with radius ' + str(active_radius) + ' ' + GV.ANGSTROM_STR + ' centred about the origin.\n')
        elif active_region_type == 'M':
          f.write('User requested the active region to be defined as the ' + str(active_region_n_mol) + ' molecules closest to the origin.\n')
        elif active_region_type == 'A':
          f.write('User requested the active region to be defined as all molecules closest to the origin without exceeding the upper limit of ' + str(active_region_n_atoms) + ' atoms.\n')
        elif active_region_type == 'N':
          f.write('User requested the active region to be defined as specific number of molecules closest to the origin for each molecule type.\n')
    elif user_input == 'N':
      user_input = None
      print (GV.dash_line) 
      user_input_2 = None
      while user_input_2 == None:
        print ('[C]luster definition')
        print ('[A]ctive region definition')
        user_input_2 = raw_input('What would you like to change [C/A]? ').strip().upper()
        if user_input_2 == 'C':
          user_input_3 = None
          while not user_input_3:
            print ('Available cluster construction definition:')
            print ('-->[R] Fixed radius')
            print ('-->[M] Specify number of each molecule type') 
            user_input_3 = raw_input('Please choose cluster size definition [R/M] --> ').strip().upper()
            if user_input_3 == 'R':
              # User defines cluster size (radius/number of atoms)
              cluster_construction_type = 'R'
              user_input_4 = None
              first_system.print_estimate_n_atoms_in_cluster()
              while not user_input_4:
                user_input_4 = raw_input('Enter the radius (in Angstroms) of the cluster(s) to be cut: ').strip()
                try:
                  radius = float(user_input_4)
                  if radius < 0:
                    user_input_4 = None
                    print ('Please enter a positive number.')
                    print (GV.dash_line)
                  else:
                    cluster_radius = radius
                    pass
                except:
                  print (str(user_input_4) + ' is not a valid radius length.')
                  print (GV.dash_line)
                  user_input_4 = None
            elif user_input_3 == 'M':
              cluster_mol_numbers_list,est_cluster_radius = CF.define_region_of_number_of_each_molecule_type('cluster',first_system,'BULK SYSTEM INFO',True,None)
              cluster_construction_type = 'M'
            else:
              print (user_input_3 + ' is not a valid option')
              user_input_3 = None
              print (GV.dash_line)
        elif user_input_2 == 'A':
          active_region_type,active_region_var = CF.chooseActiveRegionDef(cluster_radius,first_system,None)
          active_radius = None
          if active_region_type == 'S':
            active_radius = active_region_var
          elif active_region_type == 'M':
            active_region_n_mol = active_region_var
          elif active_region_type == 'A':
            active_region_n_atoms = active_region_var
          elif active_region_type == 'N':
            active_mol_number_list = active_region_var
        else:
          print ('Please enter \'C\' or \'A\'.') 
          user_input_2 = None
    else:
      user_input = None
      print ('Please enter \'Y\' or \'N\'.')
      print (GV.dash_line)

  ####################
  # Define QM region #
  ####################
  qm_region_type,qm_region_var = CF.chooseQMRegionDef(cluster_radius,first_system,log_file) 
  
  ######################
  # Loop over clusters # 
  ######################
  for i in range(n_clusters):
    cluster_dir = op.join(workdir,'embedded_cluster_' + str(i+1))
     
    ####################
    # Read bulk object #
    ####################
    bulk_obj = RI.loadSystemObj(op.join(cluster_dir,'bulk.obj'))
     
    #################### 
    # Cut bare cluster #
    #################### 
    if cluster_construction_type == 'R': 
      bare_cluster_obj = bulk_obj.cut_cluster_of_defined_radius(cluster_radius)
    elif cluster_construction_type == 'M':
      bare_cluster_obj = bulk_obj.cut_cluster_of_specified_n_of_each_molecule(cluster_mol_numbers_list,(0.0,0.0,0.0),'C')

    ff_file = op.join(cluster_dir,'ff.dat')
    bare_cluster_obj.gen_chemshell_ff_file(ff_file)
    with open(log_file,'a') as f:
      f.write('-->ChemShell force field file \'' + ff_file + '\' generated.\n')
    print ('...generated')

    ######################
    # Generate QM region #
    ######################
    if qm_region_type == 'S':
      bare_cluster_obj.det_atom_list_for_given_radius(qm_region_var,'QM')
    elif qm_region_type == 'M':
      if qm_region_var > len(bare_cluster_obj.mol_list):
        print (textwrap.fill('Requested number of molecules to be included in QM region is greater than the total number of molecules in the cluster. All molecules will be in the QM region.',80))
        bare_cluster_obj.det_atom_list_for_given_n_mol(len(bare_cluster_obj.mol_list),'QM') 
      else:
        bare_cluster_obj.det_atom_list_for_given_n_mol(qm_region_var,'QM') 
      print ('There are ' + str(len(bare_cluster_obj.qm_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
    elif qm_region_type == 'A':
      if qm_region_var > len(bare_cluster_obj.atom_list):
        print (textwrap.fill('Requested number of atoms to be included in QM region is greater than the total number of atoms in the cluster. All molecules will be in the QM region.',80))
        bare_cluster_obj.det_atom_list_for_given_n_atoms(len(bare_cluster_obj.atom_list),'QM')
      else:
        bare_cluster_obj.det_atom_list_for_given_n_atoms(qm_region_var,'QM')
      print ('There are ' + str(len(bare_cluster_obj.qm_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
    elif qm_region_type == 'N':
      if n == 0: 
        qm_mol_numbers_list = CF.define_region_of_number_of_each_molecule_type('QM region',bare_cluster_obj,'CLUSTER SYSTEM INFO',False,log_file)
      bare_cluster_obj.closest_n_molecules_of_each_type(qm_mol_numbers_list,True,True,'QM')
    elif qm_region_type == 'C':
      if active_radius == 'A':
        bare_cluster_obj.choose_qm_molecules(cluster_name,None,'QM')
      elif active_radius == None:
        bare_cluster_obj.choose_qm_molecules(cluster_name,None,'QM')
      else:
        bare_cluster_obj.choose_qm_molecules(cluster_name,active_radius,'QM')
    elif qm_region_type == 'O': 
      qm_regions = CF.choose_qm_regions(bare_cluster_obj,log_file)
      #print (qm_regions)
      bare_cluster_obj.det_atom_list_for_set_of_qm_regions(qm_regions,'QM')
      print ('There are ' + str(len(bare_cluster_obj.qm_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.qm_atom_list)) + ' atoms) in the QM region.')
    print ('The total charge of the QM region is ' + str(bare_cluster_obj.qm_charge) + '.')
    print (GV.dash_line)
    with open(log_file,'a') as f:
      f.write('-->There are ' + str(len(bare_cluster_obj.qm_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.qm_atom_list)) + ' atoms) in the QM region:\n')
      for m in bare_cluster_obj.qm_mol_list:
        f.write('---->Molecule ' + str(m) + ' (type: ' + bare_cluster_obj.mol_list[m-1].label + '; charge: ' + str(bare_cluster_obj.mol_list[m-1].charge) + '; dist. to origin: ' + str(bare_cluster_obj.mol_list[m-1].dist_to_origin) + ' ' + u'\u212B'.encode('utf-8') + ')\n')
      f.write('-->The total charge of the QM region is ' + str(bare_cluster_obj.qm_charge) + '.\n')
    
    ##########################
    # Generate active region #
    ##########################
    if active_region_type == 'S':
      if active_radius == 'A':
        pass
        #cluster_obj.gen_tiered_label('1','2','3',True)
      else:
        bare_cluster_obj.det_atom_list_for_given_radius(active_radius,'Active')
        with open(log_file,'a') as f:
          f.write('-->There are ' + str(len(bare_cluster_obj.active_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.active_atom_list)) + ' atoms) in the active region.\n')
    elif active_region_type == 'M':
      if active_region_n_mol > len(bare_cluster_obj.mol_list):
        print (textwrap.fill('Requested number of molecules to be included in active region is greater than the total number of molecules in the cluster. All molecules will be in the QM region.',80))
        bare_cluster_obj.det_atom_list_for_given_n_mol(len(bare_cluster_obj.mol_list),'Active') 
      else:
        active_list,active_mol_list = bare_cluster_obj.det_atom_list_for_given_n_mol(active_region_n_mol) 
      print ('There are ' + str(len(bare_cluster_obj.active_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.active_atom_list)) + ' atoms) in the active region.')
    elif active_region_type == 'A':
      if active_region_n_atoms > len(bare_cluster_obj.atom_list):
        print (textwrap.fill('Requested number of atoms to be included in active region is greater than the total number of atoms in the cluster. All molecules will be in the QM region.',80))
        bare_cluster_obj.det_atom_list_for_given_n_atoms(len(bare_cluster_obj.atom_list),'Active')
      else:
        bare_cluster_obj.det_atom_list_for_given_n_atoms(active_region_n_atoms,'Active')
      print ('There are ' + str(len(bare_cluster_obj.active_mol_list)) + ' molecules (' + str(len(bare_cluster_obj.active_atom_list)) + ' atoms) in the active region.')
    elif active_region_type == 'N':
      if n == 0: 
        active_mol_numbers_list = CF.define_region_of_number_of_each_molecule_type('active region',bare_cluster_obj,'CLUSTER SYSTEM INFO',False,log_file)
      bare_cluster_obj.closest_n_molecules_of_each_type(active_mol_numbers_list,True,True,'Active')
    
    ##########################################
    # Read ChemShell embedded cluster output #
    # --> Obtain point charge information    #
    ##########################################
    atom_mapping_list = []
    point_charge_list = []
    with open(op.join(cluster_dir,'embedded_cluster.pun')) as f:
      while 1:
        line = f.readline()
        #print (line)
        if len(line.split()) > 2 and line.split()[2] == 'coordinates':
          n_particles = int(line.split()[-1])
          break
     
      n_atoms = None
      for j in range(n_particles):
        line = f.readline()
        if line.split()[0] == 'F':
          if not n_atoms:
            n_atoms = j
            n_pt_charges = n_particles - n_atoms
          coordinates = line.split()[1:]
          for k,c in enumerate(coordinates):
            coordinates[k] = float(c)/bohr_in_ang
          point_charge_list.append(coordinates)
      
      while 1:
        line = f.readline()
        #print (line)
        if len(line.split()) > 2 and line.split()[2] == 'atom_charges':
          break
      for j in range(n_atoms):
        f.readline()
      for j in range(n_pt_charges):
        charge = f.readline().strip()
        point_charge_list[j].append(charge)

    ###############################################################
    # Combine with bare cluster to obtain embedded cluster object #
    ###############################################################
    embedded_cluster_obj = MC.EmbeddedClusterSystem(bare_cluster_obj,cluster_radius,point_charge_list,embedded_cluster_specs_obj.point_charge_symbol)
    RI.saveSystemObj(embedded_cluster_obj,op.join(cluster_dir,'embedded_cluster.obj'))
    output_file = op.join(cluster_dir,'embedded_cluster_' + str(i+1) + '.chm')
    embedded_cluster_obj.write_chemshell_fragment(output_file,'embedded_cluster_' + str(i+1) + '.pun')
    with open(log_file,'a') as f:
      f.write('-->ChemShell fragment input file \'' + output_file + '\' generated.\n')
    xyz_output_file = op.join(cluster_dir,'embedded_cluster_' + str(i+1) + '.xyz')
    with open(log_file,'a') as f:
      f.write('-->.xyz co-ordinates file \'' + xyz_output_file + '\' generated.\n')
    embedded_cluster_obj.write_xyz(xyz_output_file,'Cluster generated by MolCluster')

    #####################################
    # Generate tiered cluster .xyz file #
    #####################################
    embedded_cluster_obj.gen_tiered_label('1','2','3','4')
    embedded_cluster_obj.write_xyz(op.join(cluster_dir,'embedded_cluster_' + str(i+1) + '_tiered.xyz'),'Cluster generated by MolCluster',True)
    with open(log_file,'a') as f:
      f.write('-->Tiered cluster .xyz file written to \'' + op.join(cluster_dir,'embedded_cluster_tiered.xyz') + '\'.\n')

    #####################################
    # Generate QM/MM optimisation input #
    #####################################
    CSIF.qmmm_opt_input(embedded_cluster_obj,op.join(cluster_dir,'opt.chm'),'embedded_cluster.pun','embedded_cluster_opt.pun',0,1,True,log_file,'C')
    with open(log_file,'a') as f:
      f.write('-->ChemShell QM/MM optimisation input written to \'' + op.join(cluster_dir,'opt.chm') + '\'.\n')
    
    #############################################
    # Generate EmbeddedCluster data object file #
    #############################################
    RI.saveSystemObj(embedded_cluster_obj,op.join(cluster_dir,'embedded_cluster.system'))
    with open(log_file,'a') as f:
      f.write('-->EmbeddedCluster data object file written to \'' + op.join(cluster_dir,'embedded_cluster.system') + '\'.\n')
  
  print (GV.dash_line)
  print ('')
  print ('I'*80)
  print ('III' + ('  All embedded clusters generated.  ').center(74) + 'III')
  print ('I'*80)
  print ('')
  print (GV.dash_line)
  with open(log_file,'a') as f:
    f.write(GV.dash_line + '\n')
    f.write('All embedded clusters generated.\n')
    f.write(GV.exit_line)
  shutil.copy(log_file,workdir)
  CF.printBugsLine(op.join(workdir,'molcluster.log'))
  print (GV.dash_line)
  print ('O'*80)
  log_info = textwrap.wrap('Log file copied to \'' + workdir + '\'.',70)
  for l in log_info:
    print ('OOO  ' + l.center(70) + '  OOO')
  print ('O'*80)
  print (GV.dash_line)
  print (GV.exit_line)
  return

#-------------------------------------------------------------------------------
# Function main
#
# Operation: 
#-------------------------------------------------------------------------------
def main(workdir,history_file,first_system,log_file,n_configs,mol_info,timestep_size,vdw_info):
  user_input = None 
  n_atoms = len(first_system.atom_list)
  n_mols = len(first_system.mol_list)
  
  
  #############################################
  # Choose configurations to extract clusters #
  #############################################
  prompt_message = 'Would you like to generate embedded clusters from:'
  config_no_list = CF.chooseConfigurations(n_configs,prompt_message,log_file)
  print ('config_no_list: ' + str(config_no_list))
  print (GV.dash_line)

  ####################################
  # Load configurations from HISTORY #
  ####################################
  try:
    configuration_list,timestep_list,unit_cell_length_list = RI.genCoordinatesListFromHistory(history_file,mol_info,config_no_list,n_atoms,timestep_size)
  except:
    print (textwrap.fill('Error when extracting configuration(s) from HISTORY. ' + GV.bugs_line,80))
    print (GV.exit_line)
    with open(log_file,'a') as f:
      f.write('Error when extracting configuration(s) from HISTORY.\n')
      f.write(GV.exit_line)
  print (GV.dash_line)
  

  ############################
  # Build BulkSystem objects #
  ############################
  print ('Creating BulkSystem objects from configuration(s)...')
  system_list = []
  try:
    for i,c in enumerate(configuration_list):
      s = MC.BulkSystem(c,mol_info,unit_cell_length_list[i],vdw_info,timestep_list[i],(timestep_list[i]-1)*timestep_size)
      system_list.append(s)
      with open(log_file,'a') as f:
        f.write('-->Configuration extracted from timestep ' + str(timestep_list[i]) + ' (time = ' + str(s.time/1000) + ' ns).\n')
  except:
    print ('Error when building BulkSystem object(s)')
    CF.printBugsLine(log_file)
    print (GV.exit_line)
    with open(log_file,'a') as f:
      f.write('Error when building BulkSystem object(s).\n')
      f.write(GV.exit_line)
    return
  print ('...complete.')
  print (GV.dash_line)
  with open(log_file,'a') as f:
    f.write('BulkSystem object created for all configuration(s) extracted.\n')
   
  ##############################
  # Define cluster origin type #
  ##############################
  origin_type,origin_type_var = CF.chooseOrigin(None,log_file,first_system)
  print (GV.dash_line)

  #########################
  # Choose cluster radius #
  #########################
  user_input = None
  while not user_input:
    user_input = raw_input('Enter the radius (in Angstroms) of the embedded cluster(s) to be cut: ').strip()
    try:
      radius = float(user_input)
      if radius < 0:
        user_input = None
        print ('Please enter a positive number.')
        print (GV.dash_line)
      else:
        pass
    except:
      print (user_input + ' is not a number.')
      user_input = None
  
  ##################################
  # Choose projected active radius #
  ##################################
  print (GV.dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input('Enter the projected radius (in Angstroms) of the active regions: ').strip()
    try:
      active_radius = float(user_input)
      if active_radius < 0:
        user_input = None
        print ('Please enter a positive number.')
        print (GV.dash_line)
      elif active_radius > radius:
        print ('Active radius larger than the cluster radius. Active radius set equal to the cluster radius')
        active_radius = radius
      else:
        pass
    except:
      print (user_input + ' is not a number.')
      user_input = None
      print (GV.dash_line)
  
  ##########################################
  # Choose charge margin (dist to cluster) #
  ##########################################
  print (GV.dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('What would like the charge margin - distance (in Angstroms) from the cluster boundary to the outer point charges - to be?',80) + ' ').strip()
    try:
      charge_margin = float(user_input)
      if charge_margin < 0:
        user_input = None
        print ('Please enter a positive number.')
        print (GV.dash_line)
    except:
      print (user_input + ' is not a number.')
      user_input = None
      print (GV.dash_line)
  
  ###################################
  # Choose density of point charges #
  ###################################
  print (GV.dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('What would you like the charge density to be? The number of added point charges = (4 * density * density + 2). Please enter an integer > 0:',80) + ' ').strip()
    if user_input.isdigit():
      if user_input < 1:
        print ('Please enter an integer greater than 0.')
        print (GV.dash_line)
      else:
        bq_density = int(user_input)
    else:
      print (user_input + ' is not an integer.')
      user_input = None
      print (GV.dash_line)

  ###################################
  # Choose symbol for point charges #
  ###################################
  print (GV.dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input(textwrap.fill('Choose a symbol (e.g. \'X\') to represent the point charges (only alphanumerical characters allowed):',80) + ' ').strip()
    if len(user_input.split()) > 1:
      print ('Please enter a symbol without any spaces.')
      print (GV.dash_line)
      user_input = None
    else:
      if not GV.alphanumeric_p.match(user_input):
        print ('Please only use alphanumeric characters.')
        user_input = None
        print (GV.dash_line)
      else:
        charge_symbol = user_input

  ###########################
  # Choose output directory #
  ###########################
  print (GV.dash_line)
  user_input = None
  while not user_input:
    user_input = raw_input('Please enter output directory (no spaces): ').strip()
    if len(user_input.split()) > 1:
      print ('Please enter directory name without any spaces.')
      user_input = None
      print (GV.dash_line)
    else:
      output_dir = user_input
      if not op.isdir(user_input):
        print ('Creating directory \'' + user_input + '\'')
        os.makedirs(output_dir)
      with open(log_file,'a') as f:
        f.write('The output directory is \'' + user_input + '\'')
  
  print (GV.dash_line)
  print ('origin_type: ' + origin_type)
  print ('origin_var: ' + str(origin_type_var))
  print ('radius: ' + str(radius))
  print ('active_radius: ' + str(active_radius))
  print ('charge_margin: ' + str(charge_margin))
  print ('bq_density: ' + str(bq_density))
  print ('charge_symbol: ' + str(charge_symbol))
   
  embedded_cluster_specs_obj = MC.EmbeddedClusterSpecs(output_dir,configuration_list,active_radius,radius,origin_type,origin_type_var,charge_margin,bq_density,charge_symbol)
  embedded_cluster_prefix = 'embedded_cluster'

  ############################################
  # Generate input file for embedded cluster #
  ############################################
  #embedded_input = op.join(workdir,'')
  #with open():
  for n,s in enumerate(system_list):
    cluster_dir = op.join(output_dir,embedded_cluster_prefix + '_' + str(n+1))
    if not op.isdir(cluster_dir):
      print (GV.dash_line)
      print ('Creating directory \'' + cluster_dir + '\'')
      print (GV.dash_line)
      os.makedirs(cluster_dir)
      with open(log_file,'a') as f:
        f.write('-->Directory \'' + cluster_dir + '\' created.\n')
    origin_coords = s.det_origin_coords(origin_type_var,origin_type)
    s.shift_origin_to_coords(origin_coords)
    RI.saveSystemObj(s,op.join(cluster_dir,'bulk.obj'))
    chemshell_fragment_input = op.join(cluster_dir,'bulk_fragment.chm')
    bulk_fragment = 'bulk.pun'
    s.write_chemshell_fragment(chemshell_fragment_input,bulk_fragment)
    writeEmbeddedClusterInput(radius,active_radius,charge_margin,bq_density,charge_symbol,op.join(cluster_dir,'gen_embedded_cluster.chm'),bulk_fragment)
  embedded_cluster_specs_obj.embedded_cluster_prefix = embedded_cluster_prefix
  embedded_cluster_specs_obj.n_clusters = len(system_list) 
  RI.saveSystemObj(embedded_cluster_specs_obj,op.join(output_dir,'embedded_cluster_specs.obj'))
  return

#-------------------------------------------------------------------------------
# Function writeEmbeddedClusterInput
#
# Operation: Generate ChemShell input file to generate a embedded cluster
#            (correction point charges surrounding a bare cluster) from a bulk
#            periodic system.
#-------------------------------------------------------------------------------
def writeEmbeddedClusterInput(radius,active_radius,charge_margin,bq_density,charge_symbol,input_file,bulk_fragment,cluster_fragment='embedded_cluster.pun'):
  bohr_in_ang = 1.889725989
  radius = radius*bohr_in_ang
  active_radius = active_radius*bohr_in_ang
  charge_margin = charge_margin*bohr_in_ang
  with open(input_file,'w') as f:
    f.write('construct_cluster coords=' + bulk_fragment + ' cluster=' + cluster_fragment+ ' crystal_type=molecular \\\n')
    f.write('  origin_point= {0.0 0.0 0.0} radius_cluster=' + str(radius) + ' radius_active=' + str(active_radius) + '\\\n')
    f.write('  bq_margin=' + str(charge_margin) + ' bq_density=' + str(bq_density))
  return

if __name__ == "__main__":
  workdir,first_system,log_file = sys.argv
  sys.exit(main(workdir,first_system,log_file))
