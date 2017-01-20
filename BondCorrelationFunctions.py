from decimal import Decimal
import MolClass as MC
import ReadInput as RI
import CommonFunctions as CF
import GlobalVariables as GV
import textwrap
import os.path as op

def userInterface(history_file,mol_info,n_atoms,first_system,n_configs,timestep_size,ns_between_configs,workdir):
  prompt_message = 'Please choose first set of atom types'
  atom_1_type = CF.chooseSetOfAtomTypes(first_system,prompt_message)
  print ('atom_types_1: ' + str(atom_1_type))
  print (GV.dash_line)
  prompt_message = 'Please choose second set of atom types'
  atom_2_type = CF.chooseSetOfAtomTypes(first_system,prompt_message)
  print ('atom_types_2: ' + str(atom_2_type))
  print (GV.dash_line)
  
  user_input = None
  while not user_input:
    user_input = raw_input('Enter minimum distance threshold in Angstroms: ').strip()
    try:
      user_input = float(user_input)
      if user_input < 0.0:
        print ('Please enter a positive number.')
        print (GV.dash_line)
        user_input = None
      else:
        min_dist_thres = user_input
    except:
      print (user_input + ' is not a number.')
      print (GV.dash_line)
      user_input = None
    
  user_input = None
  while not user_input:
    user_input = raw_input('Enter maximum distance threshold in Angstroms (>' + str(min_dist_thres) + '): ').strip()
    try:
      user_input = float(user_input)
      if user_input < 0.0:
        print ('Please enter a positive number.')
        print (GV.dash_line)
        user_input = None
      elif user_input < min_dist_thres:
        print ('Please enter a number greater than the minimum distance threshold.')
        print (GV.dash_line)
        user_input = None
      else:
        max_dist_thres = user_input
    except:
      print (user_input + ' is not a number.')
      print (GV.dash_line)
      user_input = None

  print (GV.dash_line)
  traj_length = Decimal(ns_between_configs*(n_configs-1)).quantize(Decimal('.1'))
  user_input = None
  #print ('timestep_size = ' + str(timestep_size))
  #print ('n_configs = ' + str(n_configs))
  #print ('traj_length = ' + str(traj_length))
  while not user_input:
    user_input = raw_input('Please enter window length (in ns) for analysis (Suggested ~ ' + str(Decimal(str(traj_length/2)).quantize(Decimal('.1'))) + ' ns): ').strip()
    try:
      user_input = float(user_input)
      if user_input > traj_length:
        print('You entered a window length greater than the trajectory length (' + str(traj_length) + ')!')
        user_input = None
        print (GV.dash_line)
      elif user_input < 0.0:
        print ('You entered a window length less than 0!')
        user_input = None
        print (GV.dash_line)
      else:
        window_length_in_configs = int(user_input/float(ns_between_configs) + 1)
        if window_length_in_configs == 1:
          window_length_in_configs = 2
          print(textwrap.fill('Window length less than time between configurations in trajectory, increasing window length to ' + str(ns_between_configs) + ' ns.',80))
          if window_length_in_configs < 0.2*n_configs:
            user_input_2 = None
            while not user_input_2:
              user_input_2 = raw_input(textwrap.fill('WARNING: Window length is less than 20% of the trajectory length, do you wish to continue [Y/N]?',80) + ' ').strip().upper()
              if user_input_2 == 'Y':
                pass
              elif user_input_2 == 'N':
                user_input = None
              else:
                print ('Please enter \'Y\' or \'N\'.')
                user_input_2 = None
                print (GV.dash_line)
          elif window_length_in_configs > 0.8*n_configs:
            user_input_2 = None
            while not user_input_2:
              user_input_2 = raw_input(textwrap.fill('WARNING: Window length is more than 80% of the trajectory length (you may not have enough samples), do you wish to continue [Y/N]?',80) + ' ').strip().upper()
              if user_input_2 == 'Y':
                pass
              elif user_input_2 == 'N':
                user_input = None
              else:
                print ('Please enter \'Y\' or \'N\'.')
                user_input_2 = None
                print (GV.dash_line)
    except:
      print (user_input + ' is not a number.')
      print (GV.dash_line)
      user_input = None
  n_windows = n_configs - window_length_in_configs + 1
  print ('n_configs = ' + str(n_configs) + '; ' + str(traj_length) + ' ns.')
  print ('window_length = ' + str(window_length_in_configs) + '; ' + str((window_length_in_configs-1)*ns_between_configs) + ' ns.')
  print ('n_windows = ' + str(n_windows))
  coordinates_list,timestep_list,unit_cell_length_list = RI.genCoordinatesListFromHistory(history_file,mol_info,range(1,n_configs+1),n_atoms,timestep_size)
  config_info = {'coords_list':coordinates_list,\
                 'mol_info':mol_info,\
                 'unit_cell_length_list':unit_cell_length_list,\
                }
  calcBondCorrelation(first_system,atom_1_type,atom_2_type,min_dist_thres,max_dist_thres,window_length_in_configs,n_windows,config_info,history_file,mol_info,ns_between_configs,workdir)
  return


def calcBondCorrelation(first_system,set_1_labels,set_2_labels,bond_min,bond_max,window_length,windows,config_info,history_file,mol_info,ns_between_configs,workdir):
  n_atoms = len(first_system.atom_list)
  system_list = ret_systems_list(mol_info,history_file,n_atoms)
  #system_list[0].__dict__.keys()
  #first_system = MolClass.BulkSystem(config_info['coords_list'][0],config_info['mol_info'],config_info['unit_cell_length_list'][0],None)

  #oxygen_list = [a.atom_no for a in system_list[0].atom_list if a.label == 'OW']
  #hydrogen_list = [a.atom_no for a in system_list[0].atom_list if a.label[0] == 'H']
  set_1_list = [a.atom_no for a in first_system.atom_list if a.label in set_1_labels]
  set_2_list = [a.atom_no for a in first_system.atom_list if a.label in set_2_labels]
  
  #hbondmin = 1.65
  hbondmin = bond_min

  #hbondmax = 2.45
  hbondmax = bond_max

  # frames = length of window
  #frames = 50
  frames = window_length
  #frames = 50

  # samples = number of windows = n_configs - frames + 1
  #samples = 50
  samples = windows
  #samples = 50

  count = 0

  normalise = [0.0 for i in range(samples)]
  C = [[0.0 for i in range(frames)] for j in range(samples)]
  C_Av = [0.0 for i in range(frames)] 

  #for a1 in set_1_list:
  #  for hydrogen in hydrogen_list:
  #    count += 1

  count = len(set_1_list)*len(set_2_list)

  bonded_pair_list = []
  H = [[0 for i in range(count)] for l in range(frames)] 
  zero_array = [[0 for i in range(count)] for l in range(frames)] 
  
  items_to_analyse = samples*frames
  items_analysed = 0
  for s in range(samples):
    #system_1 = MC.BulkSystem(config_info['coords_list'][s],config_info['mol_info'],config_info['unit_cell_length_list'][s],None)
    system_1 = system_list[s]
    if s == 0:
      RI.saveSystemObj(system_1,'wtf.system')
    H = [[0 for i in range(count)] for l in range(frames)] 
#    H = zero_array
    bonded_pair_list = []
    for i,a1 in enumerate(set_1_list):
      for j,a2 in enumerate(set_2_list):
        #OHdistance = system_list[s].atom_list[a1-1].dist_to_atom(system_list[s].atom_list[a2-1])
        OHdistance = system_1.atom_list[a1-1].dist_to_atom(system_1.atom_list[a2-1])
        if OHdistance > hbondmin and OHdistance < hbondmax: 
          p = i*len(set_2_list) + j
          bonded_pair_list.append(p)
          H[0][p] = 1
#          bonded_pair_list.append(o*len(hydrogen_list)+h)

#    for p in bonded_pair_list:
#      H[0][p] = 1
    print('s=' + str(s) + '; len(bonded_pair_list)=' + str(len(bonded_pair_list)))
    normalise[s] = len(bonded_pair_list)/(1.0*count)
    for f in range(frames):
      if f == 0:
        system_2 = system_1
      else:
        #system_2 = MC.BulkSystem(config_info['coords_list'][s+f],config_info['mol_info'],config_info['unit_cell_length_list'][s+f],None)
        system_2 = system_list[s+f]
      for p in bonded_pair_list:
        j = p % len(set_2_list)
        i = (p - j) / len(set_2_list)
        OHdistance = system_2.atom_list[set_1_list[i]-1].dist_to_atom(system_2.atom_list[set_2_list[j]-1])
        if OHdistance > hbondmin and OHdistance < hbondmax: 
          H[f][p] = 1
        #elif f == 0:
        #  print ('VC: BUG!!!!!!!!')
        #  with open('BUG.txt','a') as f:
        #    f.write('Something wrong found at item ' + str(items_analysed))
        C_Av[f] += ((H[0][p]*H[f][p]))/(count*samples*normalise[s])
        #if f == 0:
        #  print ('s=' + str(s) + '; H[0][p]=' + str(H[0][p]) + '; H[f][p]=' + str(H[f][p]))
      items_analysed += 1
      #print ('Analysed ' + str(items_analysed) + ' out of ' + str(items_to_analyse) + '.')
  with open(op.join(workdir,'bond_correlation_function.txt'),'w') as out_f:
    out_f.write('Frame Time/ns C_Av\n')
    for f in range(0,frames):
      print  f, C_Av[f]
#      out_f.write(str(f) + ' ' + str(t*f) + ' ' + str(C_Av[f]) + '\n')
      out_f.write(str(f) + ' ' + str(f*ns_between_configs) + ' ' + str(C_Av[f]) + '\n')
  return

#-------------------------------------------------------------------------------
# Function ret_systems_list
#-------------------------------------------------------------------------------
def ret_systems_list(mol_info,history_file,n_atoms):
  #mol_info,vdw_info = RI.field('example_inputs4/FIELD')
  #history_file = 'example_inputs4/216short'
  n_configs,timestep_size,ns_between_configs = RI.configurationsInHistory(history_file,n_atoms)
  #config_list,timestep_list,unit_cell_len = RI.readHistory(history_file,mol_info,n_configs,1)
  coordinates_list,timestep_list,unit_cell_len_list = RI.genCoordinatesListFromHistory(history_file,mol_info,range(1,n_configs+1),n_atoms,timestep_size)
  systems_list = [MC.BulkSystem(c,mol_info,unit_cell_len_list[i],None) for i,c in enumerate(coordinates_list)]
  return (systems_list)

