import os
import textwrap
import math
import random
import copy
from GlobalVariables import *
from decimal import Decimal

#===============================================================================
# Class EmbeddedClusterSpecs
#
#       An object contains information about the various parameters used to
#       generate a embedded cluster.
#===============================================================================
class EmbeddedClusterSpecs:
  def __init__(self,workdir,config_list,r_active_in_ang,r_cluster_in_ang,origin_type,origin_var,charge_margin_in_ang,bq_density,point_charge_symbol):
    self.workdir = workdir
    self.config_list = config_list
    self.r_active_in_ang = r_active_in_ang
    self.r_cluster_in_ang = r_cluster_in_ang
    self.origin_type = origin_type
    self.origin_var = origin_var
    self.charge_margin_in_ang = charge_margin_in_ang
    self.bq_density = bq_density
    self.point_charge_symbol = point_charge_symbol
    self.n_clusters = len(config_list)
    
    return
  
  def printSpecs(self):
    print ('Specs')
    print ('Radius of clusters: ' + str(self.r_cluster_in_ang) + ' Angstroms')
    print ('Active')
    return

#===============================================================================
# Class Particle
#     
#       A particle object represents a general particle with coordinates and
#       various atomic properties. 
#       This is the parent class to the PointCharge, Atom and Mol classes.
#===============================================================================
class Particle:
  def __init__(self,label,charge,coords,unit_cell_length=None):
    self.label = label
    #self.mass = mass
    self.charge = charge
    self.coords = coords
    self.unit_cell_length = unit_cell_length
    #self.update_dist_to_origin()
    self.region = None
    return

  #-----------------------------------------------------------------------------
  # Function dist_to_coords
  #
  # Operation: Returns the distance from the particle to a given set of 
  #            coordinates.
  #-----------------------------------------------------------------------------
  def dist_to_coords(self,coords,print_coords=False):
    if print_coords == True:
      print (self.coords)
      print (coords)
      #print ('=='*40)
    r_sq = 0.0
    for i in range(3):
      i_dist = math.fabs(self.coords[i] - coords[i])
      if self.unit_cell_length and i_dist > self.unit_cell_length/2.0:
        i_dist = self.unit_cell_length - i_dist 
      r_sq += i_dist*i_dist
      if print_coords == True:
        print (i_dist)
    if print_coords == True:
      print (r_sq)
      print ('=='*40)
    return (math.sqrt(r_sq))
  
  #-----------------------------------------------------------------------------
  # Function dist_to_particle
  #
  # Operation: Returns the distance from the particle to another particle
  #-----------------------------------------------------------------------------
  def dist_to_particle(self,particle,print_coords=False):
    return (self.dist_to_coords(particle.coords,print_coords))

  #-----------------------------------------------------------------------------
  # Function dist_to_origin
  #
  # Operation: Returns the distance from the particle to the origin
  #-----------------------------------------------------------------------------
  def dist_to_origin(self):
    return (self.dist_to_coords([0.0,0.0,0.0]))
  
  #-----------------------------------------------------------------------------
  # Function coords_to_string
  #
  # Operation: Converts the coordinates of the particle to a string with 8 or
  #            user specified number of decimal places
  #-----------------------------------------------------------------------------
  def coords_to_string(self,decimal_places=8,tiered=False):
    if tiered == False:
      string = self.label.ljust(6)
    else:
      string = self.tiered_label.ljust(6)
    for i in self.coords:
      string +=  str(Decimal(str(i)).quantize(Decimal('0.' + (decimal_places-1)*'0' + '1'))).rjust(decimal_places+8)
    return (string)

  
  #-----------------------------------------------------------------------------
  # Function update_dist_to_origin
  #
  # Operation: Initialises/updates the attribute dist_to_origin 
  #-----------------------------------------------------------------------------
  def update_dist_to_origin(self):
    self.dist_to_origin = self.dist_to_coords([0.0,0.0,0.0])
    return

#===============================================================================
# Class PointCharge
#
#       A PointCharge object is a Particle with a attribute labelling the 
#       point charge number
#===============================================================================
class PointCharge(Particle):
  def __init__(self,label,charge,coords,pt_charge_no,unit_cell_length=None):
    self.pt_charge_no = pt_charge_no
    Particle.__init__(self,label,charge,coords,unit_cell_length)
    return

#===============================================================================
# Class Atom
#
#       An Atom object is a Particle with attributes describing the atom.
#===============================================================================
class Atom(Particle):
  def __init__(self,label,mass,charge,element,coords,atom_no,mol_label,unit_cell_length=None):
    #print (coords)
    #self.label = label
    self.mass = mass
    self.atom_no = atom_no
    self.orig_atom_no = atom_no
    self.element = element
    #self.coords = coords
    #self.charge = charge
    #self.unit_cell_length = unit_cell_length
    self.mol_label = mol_label
    Particle.__init__(self,label,charge,coords,unit_cell_length)
    self.update_dist_to_origin()
    return
  
  #-----------------------------------------------------------------------------
  # Function dist_to_atom
  #
  # Operation: Returns the distance between the atom to another atom.
  #-----------------------------------------------------------------------------
  def dist_to_atom(self,atom_obj):
    return (self.dist_to_particle(atom_obj))
  
#===============================================================================
# Class Mol
#      
#       A Mol object represents a molecule consisting of several Atom objects
#===============================================================================
class Mol(Particle):
  # The order of atoms in atom_list is assumed to be the same as the order that appears in the HISTORY/xyz file.
  def __init__(self,label,atom_list,mol_no,unit_cell_length=None,potential=None):
    self.unit_cell_length = unit_cell_length
    self.label = label
    self.mol_no = mol_no
    self.orig_mol_no = mol_no
    #self.tot_mass = 0
    self.mass = 0
    charge = 0
    self.charge = 0.0
    for a in atom_list:
      #self.tot_mass += float(a.mass)
      self.mass += float(a.mass)
      #self.charge += a.charge
      charge += a.charge
      a.parent_mol_obj = self
    self.atom_list = atom_list
    self.no_atoms = len(atom_list)
    if self.no_atoms > 1:
      self.uncrossUnitCell()
    self.centre_of_mass = self.retCentreOfMass(self.mass)
    Particle.__init__(self,label,charge,self.centre_of_mass,unit_cell_length)
    self.potential = potential
    #self.coords = self.centre_of_mass
    self.gen_atominfo()
    self.update_dist_to_origin()
    return

  #-----------------------------------------------------------------------------
  # Function gen_atominfo
  #
  # Operation: Fill in the atominfo attribute with information of each atom
  #-----------------------------------------------------------------------------
  def gen_atominfo(self):
    self.atominfo = {}
    for a in self.atom_list:
      if a.label not in self.atominfo:
        self.atominfo[a.label] = {'natoms':1,\
                             'mass':a.mass,\
                             'charge':a.charge,\
                             }
        self.atominfo[a.label]['natoms'] += 1
    return 
  
  #-----------------------------------------------------------------------------
  # Function uncrossUnitCell
  #
  # Operation: Ensures that atoms of molecule are at the same side of the unit 
  #            cell
  #-----------------------------------------------------------------------------
  def uncrossUnitCell(self):
    for a in self.atom_list[1:]:
      for i in range(3):
        if a.coords[i] - self.atom_list[0].coords[i] > self.unit_cell_length / 2.0:
          a.coords[i] -= self.unit_cell_length
        elif a.coords[i] - self.atom_list[0].coords[i] < - self.unit_cell_length / 2.0:
          a.coords[i] += self.unit_cell_length
    return
  
  #-----------------------------------------------------------------------------
  # Function retCentreOfMass
  #
  # Operation: Returns the centre of the mass of the Mol
  #-----------------------------------------------------------------------------
  def retCentreOfMass(self,mass=None):
    if not mass:
      mass = self.mass
    mass_coord_product = [0.0,0.0,0.0]
    centre_of_mass = []
    for a in self.atom_list:
      for i in range(3):
        mass_coord_product[i] += a.coords[i]*a.mass
    for i in range(3):
      centre_of_mass.append(mass_coord_product[i]/mass)
    out_of_bounds_flag = False
    for i in range(3):
      if centre_of_mass[i] > self.unit_cell_length/2.0:
        out_of_bounds_flag = True
        for a in self.atom_list:
          #print (a.atom_no)
          a.coords[i] -= self.unit_cell_length
      elif centre_of_mass[i] < -self.unit_cell_length/2.0:
        out_of_bounds_flag = True
        for a in self.atom_list:
          a.coords[i] += self.unit_cell_length
    if out_of_bounds_flag == True:
      #print ('Old centre of mass: ' + str(centre_of_mass))
      centre_of_mass = self.retCentreOfMass(mass)
      #print ('New centre of mass: ' + str(centre_of_mass))
    return (centre_of_mass)
  
  #-----------------------------------------------------------------------------
  # Function update_dist_to_origin
  #
  # Operation: Initialises/updates the different distance to origin attributes
  #-----------------------------------------------------------------------------
  def update_dist_to_origin(self):
    self.dist_to_origin = self.dist_to_coords([0.0,0.0,0.0])
    self.min_dist_to_origin = self.min_dist_to_coords([0.0,0.0,0.0])
    self.max_dist_to_origin = self.max_dist_to_coords([0.0,0.0,0.0])
    for a in self.atom_list:
      a.update_dist_to_origin()
    return

  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def dist_to_molecule(self):
  #  return
  
  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def dist_to_atom(self):
  #  return
  
  #-----------------------------------------------------------------------------
  # Function min_dist_to_coords
  #
  # Operation: Returns the minimum distance from the Mol to a specified set of
  #            coordinates
  #-----------------------------------------------------------------------------
  def min_dist_to_coords(self,coords):
    min_distance = 99999999
    for a in self.atom_list: 
      r = a.dist_to_coords(coords)
      if r < min_distance:
        min_distance = r
    return (min_distance) 

  #-----------------------------------------------------------------------------
  # Function max_dist_to_coords
  #
  # Operation: Returns the maximum distance from the Mol to a specified set of
  #            coordinates
  #-----------------------------------------------------------------------------
  def max_dist_to_coords(self,coords):
    max_distance = 0.0
    for a in self.atom_list: 
      r = a.dist_to_coords(coords)
      if r > max_distance:
        max_distance = r
    return (max_distance) 
  
#===============================================================================
# Class GenSystem
#
#       A GenSystem object represents a system (unit cell) with a number of
#       Mol objects
#===============================================================================
class GenSystem():
  def __init__(self,mol_info,mol_list,atom_info,atom_list,unit_cell_length=None,density=None,vdw_info=None,time_step=None,time=None):
    self.density = density
    self.mol_info = mol_info
    self.mol_list = mol_list 
    self.atom_info = atom_info
    self.atom_list = atom_list
    self.unit_cell_length = unit_cell_length
    self.vdw_info = vdw_info
    self.time_step = time_step
    self.time = time
    #self.tot_atoms = len(self.atom_list) 
    #self.tot_mols = len(self.mol_list)
    self.update_tot_charge()
    self.connectivity = {}
    for m in self.mol_info:
      self.connectivity[m['molname']] = m['connectivity']
    #print (self.connectivity)
    return
  
  #-----------------------------------------------------------------------------
  # Function closest_n_molecules_of_each_type
  #
  # atom_no_list- List of atom NUMBERS
  # mol_no_list - List of molecule NUMBERS
  #
  # Operation: Determines the closest molecules to the origin.
  #            mol_numbers_list specifies the number of each type of molecule
  #            to be included in the set of closest molecules.
  #-----------------------------------------------------------------------------
  def closest_n_molecules_of_each_type(self,mol_numbers_list,unpack=False,gen_atom_no_list=False,region_type=None):
    atom_no_list = []
    mol_no_list  = []
    region_charge = 0.0

    dist_list = []
    labels_list = []
    for m in self.mol_info:
      dist_list.append([])
      labels_list.append(m['molname'])
    for i,m in enumerate(self.mol_list):
      dist = m.dist_to_coords([0.0,0.0,0.0])
      for j,l in enumerate(labels_list):
        if m.label == l:
          dist_list[j].append((i,dist))
    
    for i,d in enumerate(dist_list):
      d = sorted(d,key=lambda x: x[1])
      if unpack == True:
        mol_no_list += [x[0]+1 for x in d[:mol_numbers_list[i]]]
        region_charge += self.mol_list[i].charge
        if gen_atom_no_list == True:
          atom_no_list += [a.atom_no for x in mol_no_list for a in self.mol_list[x-1].atom_list]
      # VC: Not currently used...
      else:
        mol_no_list.append([x[0]+1 for x in d[:mol_numbers_list[i]]])
        if gen_atom_no_list == True:
          atom_no_list.append([a.atom_no for x in mol_no_list[-1] for a in self.mol_list[x-1].atom_list])
          
    if region_type == 'QM':
      self.qm_atom_list = atom_list 
      self.qm_mol_list = mol_list 
      self.qm_charge = region_charge
      return
    elif region_type == 'Active':
      self.active_atom_list = atom_list
      self.active_mol_list = mol_list
      self.active_charge = region_charge
      return
    else:
      return (mol_no_list)

  #----------------------------------------------------------------------------- 
  # Function print_estimate_n_atoms_in_cluster
  #
  # Operation: Estimates and prints out the total number of molecules in 
  #            clusters of a specified range of radii cut from the system.
  #----------------------------------------------------------------------------- 
  def print_estimate_n_atoms_in_cluster(self,min_radius=5,max_radius=30,step=None):
    #density = len(self.atom_list)/(pow(self.unit_cell_length,3)) 
    if not step:
      step = int(self.unit_cell_length/15)
      if step < 1:
        step = 1
    if not max_radius:
      #max_radius = min_radius + 12*step
      max_radius = self.unit_cell_length*0.75
    print ('-'*80)
    print ('---' + 'Est. of # atoms in region for a range of radii'.center(74) + '---')
    print ('-'*80)
    print ('-' + ' '*78 + '-')
    print ('-' + ('Radius/' + u'\u212B').rjust(34) + ('Est. No. of Atoms').rjust(25) + ' '*19 + '-')

    r = min_radius
    while r < max_radius:
      r_str = Decimal(str(r)).quantize(Decimal('.1'))
      n_atoms = self.density*4.0/3.0*math.pi*r*r*r
      print ('-' + ' '*4 + str(r_str).rjust(30) + ('~' + str(int(n_atoms))).rjust(25) + ' '*19 + '-')
      r += step
    print ('-' + ' '*78 + '-')
    print ('-'*80)
      
    return

  #-----------------------------------------------------------------------------
  # Function update_tot_charge
  #
  # Operation: Initialises/updates the total system charge
  #-----------------------------------------------------------------------------
  def update_tot_charge(self):
    self.charge = 0.0
    for a in self.atom_list:
      self.charge += a.charge

  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def genOrderedListOfDistances(self):
  #  return
  
  #-----------------------------------------------------------------------------
  # Function: gen_chemshell_ff_file
  #
  # Operation: Generates a ChemShell force field file for the system.
  #-----------------------------------------------------------------------------
  def gen_chemshell_ff_file(self,ff_file):
    declared_atom_list = []
    with open(ff_file,'w') as f:
      for m in self.mol_info:
        for a in m['atominfo']:
          if a not in declared_atom_list:
            declared_atom_list.append(a)
            f.write('declare ' + a + '\n')
      for vdw_term in self.vdw_info:
        f.write('vdw ' + vdw_term['atom_type1'] + ' ' + vdw_term['atom_type2'] + ' ' + str(vdw_term['c6']) + ' ' + str(vdw_term['c12']) + '\n')
      for m in self.mol_info:
        for p_type in m['potential']:
          if p_type == 'bonds':
            for bond_term in m['potential'][p_type]:
              f.write('bond ' + m['atomlist'][bond_term['atom1']-1] + ' ' + m['atomlist'][bond_term['atom2']-1] + ' ' + str(bond_term['k']) + ' ' + str(bond_term['r0']) + '\n')
          elif p_type == 'angles':
            for angle_term in m['potential'][p_type]:
              f.write('angle ' + m['atomlist'][angle_term['atom1']-1] + ' ' + m['atomlist'][angle_term['atom2']-1] + ' ' + m['atomlist'][angle_term['atom3']-1] + ' ' + str(angle_term['k']) + ' ' + str(angle_term['theta0']) + '\n')
          elif p_type == 'dihedrals':
            for dihedral_term in m['potential'][p_type]:
              if dihedral_term['type'] == 'cos3':
                f.write('mm2ptor ' + m['atomlist'][dihedral_term['atom1']-1] + ' ' + m['atomlist'][dihedral_term['atom2']-1] + ' ' +  m['atomlist'][dihedral_term['atom3']-1] + ' ' + m['atomlist'][dihedral_term['atom4']-1] + ' ' + str(dihedral_term['A1']) + ' ' + str(dihedral_term['A2']) + ' ' + str(dihedral_term['A3']) + '\n')
              elif dihedral_term['type'] == 'cos':
                f.write('ptor ' + m['atomlist'][dihedral_term['atom1']-1] + ' ' + m['atomlist'][dihedral_term['atom2']-1] + ' ' +  m['atomlist'][dihedral_term['atom3']-1] + ' ' + m['atomlist'][dihedral_term['atom4']-1] + ' ' + str(dihedral_term['k']) + ' ' + str(dihedral_term['delta']) + ' ' + str(dihedral_term['m']) + 'i-j-k-l' +  '\n')
              else:
                print ('dihedral_term not recognised')
          elif p_type == 'constraints':
            pass
          else:
            print (p_type + ' potential type not recognised')
    return
  
  #-----------------------------------------------------------------------------
  # Function ret_ordered_atom_list 
  #
  # Operation: Returns the list of atoms in ascending order of distance to 
  # origin specified.
  #-----------------------------------------------------------------------------
  def ret_ordered_atom_list(self):
    for a in self.atom_list:
      a.dist_to_origin = a.dist_to_coords((0.0,0.0,0.0)) 
    ordered_atom_list = sorted(self.atom_list,key = lambda x: x.dist_to_origin)
    return (ordered_atom_list)
  
  #-----------------------------------------------------------------------------
  # Function choose_atom_origin
  #
  # Operation: This prompts the user to choose a specific atom to be an origin
  #-----------------------------------------------------------------------------
  def choose_atom_origin(self,qm_region_no):
    ordered_atom_list = self.ret_ordered_atom_list()
    atom_origin = None
    start_a = 1
    self.print_atom_list(start_a,40,ordered_atom_list) 
    atom_no = None
    while not atom_no:
      if start_a + 80 > len(ordered_atom_list):
        atom_no = raw_input('Please enter atom no. to be the origin of QM region no. ' + str(qm_region_no) + ' (1-' + str(len(ordered_atom_list)) + '). --> ').strip()
        if atom_no.isdigit():
          if int(atom_no) > 0 and int(atom_no) < len(ordered_atom_list):
            atom_no = int(atom_no)
          else:
            atom_no = None
            print (textwrap.fill('Please enter an integer between 1 and ' + str(len(ordered_omat_list)) + '.',80))
            print (dash_line)
        else:
          print (textwrap.fill(atom_no + ' is not valid. Please enter an integer between 1 and ' + str(len(ordered_atom_list)) + '.',80))
          print (dash_line)
          atom_no = None
      else:
        atom_no = raw_input('Please enter atom no. to be the origin or \'M\' to view more atoms --> ').strip()
        if atom_no in ('M','m'):
          start_a += 40*2
          self.print_atom_list(start_a,40,ordered_atom_list) 
          atom_no = None
        elif atom_no.isdigit():
          if int(atom_no) > 0 and int(atom_no) < len(ordered_atom_list):
            atom_no = int(atom_no)
          else:
            atom_no = None
            print (textwrap.fill('Please enter an integer between 1 and ' + str(len(ordered_omat_list)) + ', or \'M\' to view more atoms.',80))
            print (dash_line)
        else:
          print (textwrap.fill(atom_no + ' is not valid. Please enter an integer between 1 and ' + str(len(ordered_atom_list)) + ', or \'M\' to view more atoms.',80))
          print (dash_line)
          atom_no = None
    return (atom_no)

  #-----------------------------------------------------------------------------
  # Function print_atom_list
  #
  # Operation: Outputs list of atoms in the system to screen.
  #-----------------------------------------------------------------------------
  def print_atom_list(self,start=1,show_n_lines=40,ordered_atom_list=None):
    if not ordered_atom_list:
      ordered_atom_list = self.atom_list
    #if ascending_dist == False:
    #  work_atom_list = self.atom_list
    #elif ascending_dist == True:
    #  work_atom_list = self.order_atom_list()
    odd_flag = False
    if len(ordered_atom_list) < start + show_n_lines*2 - 1:
      show_n_lines = (len(ordered_atom_list) - start + 2)/2
      if len(ordered_atom_list)%2 == 1:
        odd_flag = True
      
    print (dash_line)
    print ('-' + (' LIST OF ATOMS IN INC. DIST. TO ORIGIN IN ' + self.name.upper() + ' ').center(78) + '-')
    print (dash_line)
    title_base = '  ' + ' No.' + ' ' + 'Label'.rjust(10) + ' ' + '  Molecule' + ' ' + ' d-to-O' + '  '
    print ('|' + title_base + '||' + title_base + '|')
    for i in range(show_n_lines):
      line = ''
      if i == show_n_lines - 1 and odd_flag == True:
        a_index = start + 2*i - 1
        a = ordered_atom_list[a_index]
        rounded_d_to_o = Decimal(str(a.dist_to_origin)).quantize(Decimal('.01'))
        line += '|  ' + str(a.atom_no).rjust(4) + ' ' + a.label.rjust(10) + ' ' + a.mol_label.rjust(10) + ' ' + str(rounded_d_to_o).rjust(7) + '  |'
        line += '|  ' + 'X'*34 + '  |'
      else:
        for j in range(2):
          a_index = start + 2*i + j - 1
          a = ordered_atom_list[a_index]
          rounded_d_to_o = Decimal(str(a.dist_to_origin)).quantize(Decimal('.01'))
          line += '|  ' + str(a.atom_no).rjust(4) + ' ' + a.label.rjust(10) + ' ' + a.mol_label.rjust(10) + ' ' + str(rounded_d_to_o).rjust(7) + '  |'
      print (line)
    print ('-'*80)
    print ('| ' + 'NOTE: d-to-O = distance to origin in Angstroms'.ljust(76) + ' |') 
    print ('-'*80)
    return

  #-----------------------------------------------------------------------------
  # Function print_mol_list
  #
  # Operation: Prints the list of moleulces in the system to screen.
  #            The molecules can be printed in ascending order of distance to
  #            origin or in ascending order of molecule label. A threshold 
  #            distance to origin can be specified to limit the number of 
  #            molecules printed.
  #-----------------------------------------------------------------------------
  def print_mol_list(self,ascending_dist=True,system_name='CLUSTER',threshold=None):
    if ascending_dist == False:
      work_mol_list = self.mol_list
    elif ascending_dist == True:
      work_mol_list = self.ordered_mol_list
    if threshold:
      work_mol_list = copy.deepcopy(work_mol_list)
      work_mol_list = [m for m in work_mol_list if m.dist_to_origin < threshold]  
    print (dash_line)
    print ('-' + (' CHOOSE MOLECULES FOR QM REGION OF ' + system_name.upper() + ' ').center(78) + '-')
    print (dash_line)
    title_base = '  ' + ' No.' + ' ' + 'Label'.rjust(10) + ' ' + 'd-to-O' + '  '
    print (title_base + '|' + title_base + '|' + title_base)

    for i in range(int(math.ceil(len(work_mol_list)/3.0))):
      line = ''
      if len(work_mol_list) - i*3 < 3:
        remainder = len(work_mol_list) - i*3
        for j in range(remainder):
          m_index = 3*i+j
          m = work_mol_list[m_index]
          rounded_d_to_o = Decimal(str(m.dist_to_origin)).quantize(Decimal('.01'))
          line += '  ' + str(m.mol_no).rjust(4) + ' ' + m.label.rjust(10) + ' ' + str(rounded_d_to_o).rjust(6) + '  '
          line += '|'
        for j in range(remainder,3):
          line += '  ' + 'X'*22 + '  ' 
          if j != 2:
            line += '|'
      else:
        for j in range(3):
          m_index = 3*i+j
          m = work_mol_list[m_index]
          rounded_d_to_o = Decimal(str(m.dist_to_origin)).quantize(Decimal('.01'))
          line += '  ' + str(m.mol_no).rjust(4) + ' ' + m.label.rjust(10) + ' ' + str(rounded_d_to_o).rjust(6) + '  '
          #print (str(m_index) +  ' ' + str(m.mol_no))
          
          if j != 2:
            line += '|'
      print (line)
    print ('-'*80)
    print ('NOTE: d-to-O = distance of centre of mass to origin in Angstroms') 
    print ('-'*80)
    return

  #-----------------------------------------------------------------------------
  # Function print_system_info
  #
  # Operation: prints information about molecules and atoms (number, type) of a 
  #            system
  #-----------------------------------------------------------------------------
  def print_system_info(self,title='SYSTEM INFORMATION',output_file=None,print_summary=True):
    lines = []
    lines.append('-'*80)
    lines.append((' ' + title + ' ').center(80,'-'))
    lines.append('-'*80)
    if self.unit_cell_length:
      lines.append(('-' + ('Periodic system ( unit cell length = ' + str(self.unit_cell_length) + ' ' + u'\u212B' + ' )').center(78) + '-').encode('utf-8'))
    else:
      lines.append('-' + 'Non-periodic system'.center(78) + '-')
    lines.append('-'*80)
    lines.append('-' + 'MOLECULES'.center(78) + '-')
    lines.append('-' + 'Label'.rjust(20) + 'Atoms'.rjust(12) + 'Number'.rjust(12) + 'Charge'.rjust(12) + 'Mass'.rjust(12) + ' '*10 + '-')
    for m in self.mol_info:
      lines.append('-' + ' '*78 + '-')
      lines.append('-' + ' '*3 + m['molname'].rjust(17) + str(m['natoms']).rjust(12) + str(m['nmols']).rjust(12) + str(round(sum(m['chargelist']),3)).rjust(12) + str(round(sum(m['masslist']),3)).rjust(12) + ' '*10 + '-')
      for a in sorted(m['atominfo'].keys()):
        lines.append('-' + '-->'.rjust(22) + a.rjust(10) + str(m['atominfo'][a]['natoms']).rjust(12) + str(round(m['atominfo'][a]['charge'],3)).rjust(12) + str(round(m['atominfo'][a]['mass'],3)).rjust(12) + ' '*10 + '-')
    lines.append('-' + ' '*78 + '-')
    lines.append('-' + 'Total Molecules'.rjust(35) + str(len(self.mol_list)).rjust(13) + ' '*30 +  '-')
    lines.append('-' + 'Total Atoms'.rjust(35) + str(len(self.atom_list)).rjust(13) + ' '*30 +  '-')
    lines.append('-' + 'Total Charge'.rjust(35) + str(round(self.charge,3)).rjust(13) + ' '*30 +  '-')
    lines.append('-'*80)
    if output_file:
      with open(output_file,'w') as f:
        for line in lines:
          f.write(line + '\n')
        if print_summary == True:
          print ('O'*80)
          output_lines = textwrap.wrap('System information written to ' + output_file,70)
          for l in output_lines:
            print ('OOO  ' + l.center(70) + '  OOO')
          #print (('  System information written to ' + output_file + '  ').center(80,'O'))
          print ('O'*80)
        
    else:
      for line in lines:
        print (line)
    return
  
  #-----------------------------------------------------------------------------
  # Function write_xyz
  #
  # Operation: Writes system configuration to an xyz file
  #-----------------------------------------------------------------------------
  def write_xyz(self,xyz_file='system.xyz',title=None,tiered=False,append=False):
    if append == True:
      f_option = 'a'
    else:
      f_option = 'w'
    with open (xyz_file,f_option) as f:
      f.write(str(len(self.atom_list)) + '\n')
      if not title:
        f.write('Generated by MolCluster.py\n')
      else:
        f.write(title + '\n')
      for a in self.atom_list:
        f.write(a.coords_to_string(8,tiered) + '\n')
    return

  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def write_oniom(self,gaussian_file='system.com'):
  #  return
  
  #-----------------------------------------------------------------------------
  # Function generate_connectivity_file
  #
  # Operation: Generates a file containing the connectivity information that
  #            can be used in the ChemShell fragment file.
  #-----------------------------------------------------------------------------
  def generate_connectivity_file(self,output_dir,output_file='conn.txt'):
    atom_index = 0
    with open(os.path.join(output_dir,output_file),'w') as f:
      for m in self.mol_list:
        conn_info = self.connectivity[m.label]
        for c in conn_info:
          atom_1 = atom_index + c[0] 
          atom_2 = atom_index + c[1]
          f.write(str(atom_1) + ' ' + str(atom_2) + '\n')
        atom_index += len(m.atom_list)
    return

  #-----------------------------------------------------------------------------
  # Function write_chemshell_fragment
  #
  # Operation: Outputs a ChemShell input file to generate a fragment
  #-----------------------------------------------------------------------------
  def write_chemshell_fragment(self,chemshell_file='system.chm',pun_file='bulk_system.pun'):
    #chemshell_file = os.path.join(output_dir,chemshell_file)
    output_dir = os.path.dirname(chemshell_file)
    with open(chemshell_file,'w') as f:
      f.write('c_create coords=' + pun_file + ' {\n')
      if self.unit_cell_length:
        f.write('space_group\n')
        f.write('1\n')
        f.write('cell_constants angstrom\n')
        ucl_dec = Decimal(str(self.unit_cell_length)).quantize(Decimal('1.000000'))
        f.write(str(ucl_dec) + ' ' + str(ucl_dec) + ' ' + str(ucl_dec) + ' 90.0 90.0 90.0\n')
        f.write('coordinates\n')
      else:
        f.write('coordinates angstrom\n')
      for a in self.atom_list:
        atom_line = a.label + ' ' 
        for c in a.coords:
          if self.unit_cell_length:
            atom_line += str(Decimal(str(c/self.unit_cell_length)).quantize(Decimal('1.000000'))) + ' '
          else:
            atom_line += str(Decimal(str(c)).quantize(Decimal('1.000000'))) + ' '
        atom_line += 'charge ' + str(a.charge)
        f.write(atom_line + '\n')
      f.write('}\n')
    self.generate_connectivity_file(output_dir)
    
    return
    
#===============================================================================
# Class BulkSystem
#
#       A BulkSystem represents a 3D periodic system.
#===============================================================================
# Length units assumed to be Angstrom
class BulkSystem(GenSystem):
  def __init__(self,coords_array,mol_info,unit_cell_length,vdw_info=None,time_step=None,time=None):
    self.system_type = 'Bulk'
    mol_list = []
    atom_list = []
    atom_info = {}
    # Loop over molecule types
    for t in mol_info:
      # Loop over number of molecules for each type
      for m in range(t['nmols']):
        mol_atom_list = []
        # Loop over number of atoms for each molecule
        for a in range(t['natoms']):
          atom_label = t['atomlist'][a]
          atom_obj = Atom(atom_label,t['masslist'][a],t['chargelist'][a],None,coords_array[len(mol_list)][a],len(atom_list)+1,t['molname'],unit_cell_length)
          mol_atom_list.append(atom_obj)  
          atom_list.append(atom_obj)
          if atom_label not in atom_info:
            atom_info[atom_label] = {'natoms':1,\
                                     'mass':t['masslist'][a],\
                                     'charge':t['chargelist'][a]}
          else:
            atom_info[atom_label]['natoms'] += 1
        #mol_count += 1
        mol_obj = Mol(t['molname'],mol_atom_list,len(mol_list)+1,unit_cell_length,t['potential'])
        mol_list.append(mol_obj)
    density = len(atom_list)/(pow(unit_cell_length,3))
    GenSystem.__init__(self,mol_info,mol_list,atom_info,atom_list,unit_cell_length,density,vdw_info,time_step,time)
    #self.print_system_info()
    return
  
  #-----------------------------------------------------------------------------
  # Function shift_origin_to_coords
  #
  # Operation: Shift co-ordinates such that the origin is (0.0,0.0,0.0)
  #-----------------------------------------------------------------------------
  def shift_origin_to_coords(self,origin_coords):
    origin_coords = copy.deepcopy(origin_coords)
    for m in self.mol_list:
      for i in range(3):
        for a in m.atom_list:
          a.coords[i] = a.coords[i] - origin_coords[i]
        m.centre_of_mass = m.retCentreOfMass()
        m.coords = m.centre_of_mass
        #print ('COM: ' + str(m.centre_of_mass))
        #print ('Coords: ' + str(m.coords))
        #print ('===')
        if self.unit_cell_length and m.centre_of_mass[i] > self.unit_cell_length/2.0:
          for a in m.atom_list:
            a.coords[i] -= self.unit_cell_length
          m.centre_of_mass = m.retCentreOfMass()
          m.coords = m.centre_of_mass
        elif self.unit_cell_length and m.centre_of_mass[i] < -self.unit_cell_length/2.0:
          for a in m.atom_list:
            a_coords[i] += self.unit_cell_length
          m.centre_of_mass = m.retCentreOfMass()
          m.coords = m.centre_of_mass
      m.update_dist_to_origin()
    return

  #----------------------------------------------------------------------------- 
  # Function print_estimate_active_and_frozen_atoms
  #
  # Operation: Estimates and prints to screen the number of active atoms 
  #            for various active radii within a cluster cut from the system
  #            with a specified size .
  #----------------------------------------------------------------------------- 
  def print_estimate_active_and_frozen_atoms(self,min_radius=5,cluster_radius=15):
    est_tot_atoms = self.density*4.0/3.0*math.pi*cluster_radius*cluster_radius*cluster_radius
    step = cluster_radius/12
    print ('-'*80)
    print ('---' + 'Est. of # active and frozen atoms for a range of active radii'.center(74) + '---')
    print ('-'*80)
    print ('-' + ' '*78 + '-')
    print ('-' + ('Active Radius/' + u'\u212B').rjust(25) + ('Active Atoms').rjust(19) + ('Frozen Atoms').rjust(19) + ' '*15 + '-')

    r = min_radius
    while r < cluster_radius:
      r_str = Decimal(str(r)).quantize(Decimal('.1'))
      n_active_atoms = self.density*4.0/3.0*math.pi*r*r*r
      n_frozen_atoms = est_tot_atoms - n_active_atoms
      print ('-' + ' '*4 + str(r_str).rjust(21) + ('~' + str(int(n_active_atoms))).rjust(19) + ('~' + str(int(n_frozen_atoms))).rjust(19) + ' '*15 + '-')
      r += step
    print ('-' + ' '*78 + '-')
    print ('-'*80)
      
    return

  
  #-----------------------------------------------------------------------------
  # Function det_origin_coords
  #
  # Operation: Initialises the origin coordinates based on the origin type.
  #-----------------------------------------------------------------------------
  def det_origin_coords(self,origin,origin_type):
    if origin_type == 'C':
      origin_coords = origin
    elif origin_type == 'MT':
      for m in self.mol_list:
        if origin == m.label:
          origin_coords = copy.deepcopy(m.centre_of_mass)
          break
    elif origin_type == 'AT':
      for a in self.atom_list:
        #print (a.label)
        if origin == a.label:
          origin_coords = copy.deepcopy(a.coords)
          break
    elif origin_type == 'MN':
      origin_coords = copy.deepcopy(self.mol_list[origin-1].centre_of_mass)
    elif origin_type == 'AN':
      origin_coords = copy.deepcopy(self.atom_list[origin-1].coords)
    return (origin_coords)

  #-----------------------------------------------------------------------------
  # Function cut_cluster_of_defined_radius
  # 
  # Operation: Cuts a cluster with a specified radius out of the system.
  #
  # origin_type can be:
  # --> C for coordinates (xyz coordinates list)
  # --> MT for molecule name (string) - the first molecule with that name will 
  #     be the origin
  # --> AT for atom name (string) - the first atom with that name will be the 
  #     origin
  # --> MN for molecule number (integer)
  # --> AN for atom number (integer)
  #
  # radius_type can be:
  # --> S for soft - atoms can be outside the boundary as long as the centre of 
  #     mass of the molecule is within the boundary
  # --> H for hard - no atoms whatsoever can be outside the boundary
  #-----------------------------------------------------------------------------
  def cut_cluster_of_defined_radius(self,radius,origin=[0.0,0.0,0.0],origin_type='C',radius_type='S',force_neutral=False,print_summary=False):

    origin_coords = self.det_origin_coords(origin,origin_type)
    
    # Shift the system so that the origin is at [0.0,0.0,0.0]
    self.shift_origin_to_coords(origin_coords)
    if radius*2 > self.unit_cell_length:
      supercell_dim = int(math.ceil(radius*2/self.unit_cell_length))
      print ('Generating ' + str(supercell_dim) + 'x' + str(supercell_dim) + 'x' + str(supercell_dim) + ' supercell...')
      self.supercell(supercell_dim)
    
    cluster_atom_list = []
    cluster_atom_info = {}
    cluster_mol_list = []
    cluster_mol_info = []
    cluster_dist_list = []
    cluster_mol_labels = []
    dist_list = []
    for i,m in enumerate(self.mol_list):
      if radius_type == 'S':
        dist = self.mol_list[i].dist_to_coords([0.0,0.0,0.0])
      elif radius_type == 'H':
        dist = self.mol_list[i].max_dist_to_coords([0.0,0.0,0.0])
      dist_list.append(dist)
      #print (str(i) + ': ' + str(dist))
      if dist < radius:
        m.dist_to_origin = dist
        if not cluster_mol_info or m.label not in cluster_mol_labels:
          cluster_mol_labels.append(m.label)

          cluster_mol_info.append({'molname':m.label,\
                                   'nmols':1,\
                                   'natoms':len(m.atom_list),\
                                   'atomlist':[a.label for a in m.atom_list],\
                                   'masslist':[a.mass for a in m.atom_list],\
                                   'chargelist':[a.charge for a in m.atom_list],\
                                   'connectivity':self.connectivity[m.label],\
                                   'potential':m.potential,\
                                   'atominfo':m.atominfo})
        else:
          cluster_mol_info[-1]['nmols'] += 1

        mol_copy = copy.deepcopy(self.mol_list[i])
        mol_copy.dist_to_origin = dist
        cluster_mol_list.append(mol_copy)
        for a in mol_copy.atom_list:
          cluster_atom_list.append(a)
          if a.label not in cluster_atom_info:
            cluster_atom_info[a.label] = {'natoms':1,\
                                  'mass':a.mass,\
                                  'charge':a.charge}
          else:
            cluster_atom_info[a.label]['natoms'] += 1
    #print (self.vdw_info)
    #print (self.density)
    cluster_obj = ClusterSystem(cluster_mol_info,cluster_mol_list,cluster_atom_info,cluster_atom_list,radius,radius_type,self.density,self.vdw_info)
    if force_neutral == True:
      cluster_obj.make_cluster_neutral()
    if print_summary == True:
      # Print summary
      print ('-'*80)
      print (' SUMMARY '.center(80,'-'))
      print ('-'*80)
      if radius_type == 'S':
        radius_type_str = 'soft'
      elif radius_type == 'H':
        radius_type_str = 'hard'

      print ('-' + ('A cluster of (' + radius_type_str + ') radius ' + str(radius) +  u'\u212B').center(78) + '-')
      print ('-' + 'has been cut from the bulk system.'.center(78) + '-')
    
      print ('-'*80)
    return (cluster_obj)

  #-----------------------------------------------------------------------------
  # Function cut_cluster_of_specified_n_of_each_molecule
  #
  # Operation: Cuts a cluster from the system with a specified number of each
  #            type of molecule
  #-----------------------------------------------------------------------------
  def cut_cluster_of_specified_n_of_each_molecule(self,cluster_mol_n_list,origin=[0.0,0.0,0.0],origin_type='C',print_summary=False):

    origin_coords = self.det_origin_coords(origin,origin_type)
    
    # Shift the system so that the origin is at [0.0,0.0,0.0]
    self.shift_origin_to_coords(origin_coords)
    max_supercell_dim = 1
    for i,m in enumerate(self.mol_info):
      if int(math.ceil(float(cluster_mol_n_list[i])/m['nmols'])) > max_supercell_dim*max_supercell_dim*max_supercell_dim:
        max_supercell_dim = int(math.ceil(math.pow(float(cluster_mol_n_list[i])/m['nmols'],1.0/3.0)))
    if max_supercell_dim > 1:
      self.supercell(max_supercell_dim)
      print ('Generating ' + str(max_supercell_dim) + 'x' + str(max_supercell_dim) + 'x' + str(max_supercell_dim) + ' supercell...')
    
    cluster_atom_list = []
    cluster_atom_info = {}
    cluster_mol_list = []
    cluster_mol_info = []
    cluster_dist_list = []
    cluster_mol_labels = []

    cluster_mol_number_list = self.closest_n_molecules_of_each_type(cluster_mol_n_list,False)
    for i,n in enumerate(cluster_mol_n_list):
      if n > 0:
        cluster_mol_info.append({'molname':self.mol_info[i]['molname'],\
                                 'nmols':n,\
                                 'natoms':self.mol_info[i]['natoms'],\
                                 'atomlist':self.mol_info[i]['atomlist'],\
                                 'masslist':self.mol_info[i]['masslist'],\
                                 'chargelist':self.mol_info[i]['chargelist'],\
                                 'connectivity':self.mol_info[i]['connectivity'],\
                                 'potential':self.mol_info[i]['potential'],\
                                 'atominfo':self.mol_info[i]['atominfo']})
        for j in cluster_mol_number_list[i]:
          mol_copy = copy.deepcopy(self.mol_list[j-1])
          mol_copy.update_dist_to_origin()

          cluster_mol_list.append(mol_copy)
          for a in mol_copy.atom_list:
            cluster_atom_list.append(a)
            if a.label not in cluster_atom_info:
              cluster_atom_info[a.label] = {'natoms':1,\
                                            'mass':a.mass,\
                                            'charge':a.charge}
            else:
              cluster_atom_info[a.label]['natoms'] += 1
    cluster_obj = ClusterSystem(cluster_mol_info,cluster_mol_list,cluster_atom_info,cluster_atom_list,None,None,self.density,self.vdw_info)
    if print_summary == True:
      # Print summary
      print ('-'*80)
      print (' SUMMARY '.center(80,'-'))
      print ('-'*80)
      #if radius_type == 'S':
      #  radius_type_str = 'soft'
      #elif radius_type == 'H':
      #  radius_type_str = 'hard'

      print ('-' + ('A cluster of XXX').center(78) + '-')
      print ('-' + 'has been cut from the bulk system.'.center(78) + '-')
      print ('-'*80)
    return (cluster_obj)

  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def cut_cluster_of_n_molecules(self,origin,origin_type,n_molecules,force_neutral=False):
  #  origin_coords = self.det_origin_coords(origin,origin_type)
  #  return
  
  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def cut_cluster_of_n_atoms(self,origin,origin_type,n_atoms,force_neutral=False):
  #  origin_coords = self.det_origin_coords(origin,origin_type)
  #  return
  
  #-----------------------------------------------------------------------------
  # Function rank_molecules_by_dist_to_origin
  #
  # Operation: Returns a list of molecules in the system ranked by ascending
  #            distance away from the origin.
  #-----------------------------------------------------------------------------
  def rank_molecules_by_dist_to_origin(self,origin=[0.0,0.0,0.0]):
    dist_list = []
    for m in self.mol_list:
      dist_list.append((m.dist_to_coords(origin),m))
    dist_list =  sorted(dist_list,key=lambda x: x[0], reverse=True)
    #print (dist_list)
    return (dist_list)
  
  #-----------------------------------------------------------------------------
  # Function supercell
  #
  # Operation: Creates a supercell of the system repeated in each dimension 
  #            by a specified number.
  #-----------------------------------------------------------------------------
  def supercell(self,repeats=2):
    tot_images = repeats*repeats*repeats
    new_mol_list = []
    orig_unit_cell_length = self.unit_cell_length
    self.unit_cell_length = self.unit_cell_length*repeats
    #print (orig_unit_cell_length) 
    for i in range(0,repeats):
      for j in range(0,repeats):
        for k in range(0,repeats):
          if i + j + k == 0:
            continue
          for m in self.mol_list:
            #if m.label == 'cop':
            #  print (m.atom_list[0].coords)
            m.unit_cell_length = self.unit_cell_length
            for a in m.atom_list:
              a.unit_cell_length = self.unit_cell_length
            #print (m.atom_list[0].unit_cell_length)
            m_image = copy.deepcopy(m)
            m_image.unit_cell_length = self.unit_cell_length
            # Change co-ordinates, centre of mass, atom_list, mol_list
            for a in m_image.atom_list:
              a.coords[0] += i*orig_unit_cell_length
              a.coords[1] += j*orig_unit_cell_length
              a.coords[2] += k*orig_unit_cell_length 
              a.atom_no = len(self.atom_list) + 1
              self.atom_list.append(a) 
              #print (self.atom_list[-1].coords)
            #print ('old' + str(a.coords))
            #if m.label == 'cop':
            #  print (m_image.atom_list[0].coords)
            m_image.centre_of_mass = m_image.retCentreOfMass()
            m_image.coords = m_image.centre_of_mass
            #print ('new' + str(a.coords))
            m_image.mol_no = len(self.mol_list) + len(new_mol_list) + 1
            new_mol_list.append(m_image)
    self.mol_list = self.mol_list + new_mol_list
    self.charge = self.charge*tot_images
    for a in self.atom_info:
      self.atom_info[a]['natoms'] = self.atom_info[a]['natoms']*tot_images
    for m in self.mol_info:
      m['nmols'] = m['nmols']*tot_images
    
    return


#===============================================================================
# Class ClusterSystem
#
#       A ClusterSystem object represents a cluster cut out from a BulkSystem
#===============================================================================
class ClusterSystem(GenSystem):

  def __init__(self,mol_info,mol_list,atom_info,atom_list,radius,radius_type,density=None,vdw_info=None):
    self.system_type = 'Cluster'
    GenSystem.__init__(self,mol_info,mol_list,atom_info,atom_list,None,density,vdw_info)
    self.radius = radius
    self.radius_type = radius_type
    #self.dist_list = dist_list
    
    self.reset_atom_and_mol_no()
    self.gen_ordered_mol_list()
    self.qm_atom_list = []
    self.qm_mol_list = []
    self.active_atom_list = []
    self.active_mol_list = []
    return
  
  #-----------------------------------------------------------------------------
  # Function gen_ordered_mol_list
  #
  # Operation: Generates a list of molecules ordered by distance to origin
  #-----------------------------------------------------------------------------
  def gen_ordered_mol_list(self):
    self.ordered_mol_list = sorted(self.mol_list,key = lambda x: x.dist_to_origin)
    return

  #-----------------------------------------------------------------------------
  # Function reset_atom_an_mol_no
  #
  # Operation: Updates the atom and molecules labels. This is needed as the 
  #            atoms and molcules of the ClusterSystem are a subset of the 
  #            BulkSystem
  #-----------------------------------------------------------------------------
  def reset_atom_and_mol_no(self):
    m_index = 1
    a_index = 1
    self.orig_to_new_atom_no_dict = {}
    self.new_to_orig_atom_no_dict = {}
    for m in self.mol_list:
      m.orig_mol_no = m.mol_no
      m.mol_no = m_index
      for a in m.atom_list:
        #a.orig_atom_no = a.atom_no
        a.atom_no = a_index
        a_index += 1
        self.orig_to_new_atom_no_dict[a.orig_atom_no] = a.atom_no
        self.new_to_orig_atom_no_dict[a.atom_no] = a.orig_atom_no
      m_index += 1
    return
  
  #-----------------------------------------------------------------------------
  # Function closes_n_molecules_of_single_type
  # 
  # Operation: Returns the INDEX of the n_mol molecules in mol_list closest to 
  #            the origin
  #-----------------------------------------------------------------------------
  def closest_n_molecules_of_single_type(self,mol_label,n_mol):
    dist_list = []
    for i,m in enumerate(self.mol_list):
      if m.label == mol_label:
        dist = m.dist_to_coords([0.0,0.0,0.0])
        dist_list.append((i,dist))
    dist_list = sorted(dist_list,key=lambda x: x[1])

    if n_mol > len(dist_list):
      print ('***WARNING: requested n_mol > number of ' + mol_label + ' molecules in system.***')
      n_mol = len(dist_list)

    return ([i[0] for i in dist[0:n_mol]])


  #-----------------------------------------------------------------------------
  # Function choose_molecules
  #
  # Operation: User interface for custom choice of molecules in a region
  #
  # atom_list  - List of atom NUMBERS to be included in the region
  # mol_list   - List of molecule NUMBERS to be included in the region
  #-----------------------------------------------------------------------------
  def choose_qm_molecules(self,system_name='CLUSTER',threshold=None,region_type='QM'):
    self.print_mol_list(True,system_name,threshold)
    #for m in self.mol_list:
    #  for a in m.atom_list:
    #    print (a.atom_no)
    user_input = None
    atom_list = []
    if region_type == 'QM':
      region_string = 'QM region'
    while not user_input:
      user_input = raw_input(textwrap.fill('Please choose molecules to include in the ' + region_string + '. Enter molecule numbers (integers between 1 and ' + str(len(self.mol_list)) + ') separated by space (e.g. "1 2 4 15"):',80) + ' ')
      mol_nos = user_input.split()
      mol_list = []
      mol_types = {}
      n_atoms = 0
      #self.qm_charge = 0.0
      region_charge = 0.0
      for n in mol_nos:
        if n.isdigit() == False:
          print (n + ' is not a valid molecule number.')
          print (dash_line)
          user_input = None
          break
        elif int(n) < 1 or int(n) > len(self.mol_list):
          print (n + ' is not within the range of molecule numbers (1-' + str(len(self.mol_list)) + ').')
          print (dash_line)
          user_input = None
          break
        else:
          if int(n) in mol_list:
            print ('Molecule no. ' + str(n) + ' already chosen, removing duplicate.')
          else:
            mol_list.append(int(n))
            if self.mol_list[int(n)-1].label not in mol_types:
              mol_types[self.mol_list[int(n)-1].label] = 1
            else:
              mol_types[self.mol_list[int(n)-1].label] += 1
            n_atoms += len(self.mol_list[int(n)-1].atom_list)
            region_charge += self.mol_list[int(n)-1].charge
      if user_input:
        print (dash_line)
        print ('You have included the following in the ' + region_string + ':')
        for m in mol_types:
          print ('--> ' + str(mol_types[m]) + ' x ' + m)
        print ('I'*80)
        print ('III' + ('  There are ' + str(n_atoms) + ' atoms in the QM region (tot. charge = ' + str(Decimal(region_charge).quantize(Decimal('.01'))) + ' e).  ').center(74) + 'III')
        print ('I'*80)
        confirm = None
        while not confirm:
          confirm = raw_input('Do you wish to continue [Y/N]? ')
          if confirm.upper() == 'Y':
            print (dash_line)
          elif confirm.upper() == 'N':
            user_input = None
            print (dash_line)
          else:
            confirm = None
            print ('Please enter \'Y\' or \'N\'')
            print (dash_line)

    for m in mol_list:
      for a in self.mol_list[m-1].atom_list:
        atom_list.append(a.atom_no)
    if region_type == 'QM':
      self.qm_atom_list = atom_list
      self.qm_mol_list = mol_list
      self.qm_charge = region_charge
    return 
  
  #-----------------------------------------------------------------------------
  # Function det_atom_list_for_set_of_qm_regions
  #
  # Operation: Generates the list of QM/Active atoms and molcules in the
  #            overall QM/Active region given a list of QM/Active sub regions
  #-----------------------------------------------------------------------------
  def det_atom_list_for_set_of_qm_regions(self,qm_regions,region_type=None):
    qm_atom_list = []
    qm_mol_list = []
    region_charge = 0.0
    for i,m in enumerate(self.mol_list):
      for r in qm_regions:
        if m.dist_to_coords(r[0]) < float(r[1]):
          qm_mol_list.append(m.mol_no)
          region_charge += m.charge
          for a in m.atom_list:
            qm_atom_list.append(a.atom_no)
          break
    if region_type == 'QM':
      self.qm_atom_list = qm_atom_list
      self.qm_mol_list = qm_mol_list
      self.qm_charge = region_charge
    elif region_type == 'Active':
      self.active_atom_list = qm_atom_list
      self.active_mol_list = qm_mol_list
      self.active_charge = region_charge
    return 

  #-----------------------------------------------------------------------------
  # Function det_atom_list_for_given_n_mol
  #
  # Operation: Generates the list of atoms and molcules in the QM/Active region 
  #            given a specified number of molecules to be included in the 
  #            region.
  #-----------------------------------------------------------------------------
  def det_atom_list_for_given_n_mol(self,n_mol,region_type=None):
    atom_list = []
    mol_list = []
    region_charge = 0.0
    for i in range(n_mol):
      mol_list.append(self.ordered_mol_list[i].mol_no)
      for a in self.ordered_mol_list[i].atom_list:
        atom_list.append(a.atom_no)
      region_charge += self.ordered_mol_list[i].charge
      
    if region_type == 'QM':
      self.qm_atom_list = atom_list
      self.qm_mol_list = mol_list
      self.qm_charge = region_charge
    elif region_type == 'Active':
      self.active_atom_list = atom_list
      self.active_mol_list = mol_list
      self.active_charge = region_charge
    return 

  #-----------------------------------------------------------------------------
  # Function det_atom_list_for_given_n_atoms
  #
  # Operation: Generates the list of atoms and molcules in the QM/Active region 
  #            given a specified number of molecules to be included in the 
  #            region
  #-----------------------------------------------------------------------------
  def det_atom_list_for_given_n_atoms(self,n_atoms,region_type=None):
    atom_list = []
    mol_list = []
    region_charge = 0.0
    if region_type == 'QM':
      region_string = 'QM region'
    elif region_type == 'Active':
      region_string = 'active region'
    if len(self.ordered_mol_list[0].atom_list) > n_atoms:
      print (textwrap.fill('The upper limit of atoms to be included in the ' + region_string + ' is greater than the number of atoms in the molecule closest to the origin. The ' + region_string + ' will consist of that molecule.',80))
      mol_list.append(self.mol_list[0].mol_no)
      region_charge += self.mol_list[0].charge
      for a in self.mol_list[0].atom_list:
        atom_list.append(a.atom_no)
    else:
      for m in self.ordered_mol_list:
        if len(atom_list) + len(m.atom_list) > n_atoms:
          break
        else:
          mol_list.append(m.mol_no)
          region_charge += m.charge
          for a in m.atom_list:
            atom_list.append(a.atom_no)
    if region_type == 'QM':
      self.qm_atom_list = atom_list
      self.qm_mol_list = mol_list
      self.qm_charge = region_charge
    elif region_type == 'Active':
      self.active_atom_list = atom_list
      self.active_mol_list = mol_list
      self.active_charge = region_charge
    return 

  #-----------------------------------------------------------------------------
  # Function det_atom_list_for_given_radius
  #
  # Operation: Generates the list of atoms and molcules in the QM/Active region 
  #            given the radius of the region 
  #-----------------------------------------------------------------------------
  def det_atom_list_for_given_radius(self,radius,region_type=None):
    atom_list = []
    mol_list = []
    atom_index = 0
    region_charge = 0.0
    for i,m in enumerate(self.mol_list):
      if m.dist_to_origin < radius:
        atom_list += range(atom_index+1,atom_index+len(m.atom_list)+1)
        mol_list.append(i+1)
      atom_index += len(m.atom_list)
      region_charge += m.charge
    if region_type == 'QM':
      self.qm_atom_list = atom_list 
      self.qm_mol_list = mol_list 
      self.qm_charge = region_charge
    elif region_type == 'Active':
      self.active_atom_list = atom_list
      self.active_mol_list = mol_list
      self.active_charge = region_charge
    return 

  #-----------------------------------------------------------------------------
  # Function make_cluster_neutral
  #
  # Operation: Force neutrality of the cluster by removing molecules
  #-----------------------------------------------------------------------------
  def make_cluster_neutral(self,random=False,log_file=None):
    #print ('make_cluster_neutral')
    deletion_order_dict = self.gen_deletion_order_dict(random)
    to_delete_list = self.gen_delete_list(deletion_order_dict)
    if len(to_delete_list) > 0:
      print ('To force neutrality the following molecules have been removed:')
       
      for m in to_delete_list:
        print ('-->Molecule ' + str(m.mol_no) + ' (type: ' + m.label + '; charge: ' + str(m.charge) + '; dist. to origin: ' + str(m.dist_to_origin) + u'\u212B' + ')')
        for i,mol_d in enumerate(self.mol_info):
          if mol_d['molname'] == m.label:
            self.mol_info[i]['nmols'] -= 1
        for a in m.atom_list:
          self.atom_info[a.label]['natoms'] -= 1
          self.atom_list.remove(a)
        self.mol_list.remove(m)
      print (dash_line)
    self.update_tot_charge()
    self.reset_atom_and_mol_no()
    #print (to_delete_list)
    return 
   
  #-----------------------------------------------------------------------------
  # Function gen_delete_list
  #
  # Operation: Algorithm to determine the list of molecules to be deleted to
  #            force neutrality of a cluster
  #-----------------------------------------------------------------------------
  def gen_delete_list(self,deletion_order_dict):
    avail_pos_chg = []
    avail_neg_chg = []
    delete_list = []
    tot_chg = self.charge
    #print ('tot_chg' + str(tot_chg))
    for c in deletion_order_dict:
      if c > 0.0:
        avail_pos_chg.append(c)
      else:
        avail_neg_chg.append(-c)
    #print ('pos' + str(avail_pos_chg))
    #print ('neg' + str(avail_neg_chg))
    avail_pos_chg = sorted(avail_pos_chg)
    avail_neg_chg = sorted(avail_neg_chg)
    while math.fabs(tot_chg) > 0.0001:
      chg_to_del = None
      if tot_chg > 0.0:
        for c in avail_pos_chg:
          # Test whether tot_chg can be divided exactly by c
          if math.fabs(tot_chg / float(c) - round(tot_chg/float(c))) < 0.00001:
            #
            if len(deletion_order_dict[c]) >= round(tot_chg/float(c)):
              n = int(round(tot_chg/float(c)))
              delete_list += deletion_order_dict[c][:n]
              tot_chg -= float(c)*n
              break
            else:
              n = len(deletion_order_dict[c])
              delete_list += deletion_order_dict[c][:n]
              tot_chg -= float(c)*n
          elif tot_chg > c:
            chg_to_del = c
          elif c == avail_pos_chg[-1] and not chg_to_del:
            chg_to_del = c
      else:
        abs_chg = math.fabs(tot_chg)
        for c in avail_neg_chg:
          if math.fabs(abs_chg / float(c) - round(abs_chg/float(c))) < 0.00001:
            if len(deletion_order_dict[c]) >= round(abs_chg/float(c)):
              n = int(round(abs_chg/float(c)))
              delete_list += deletion_order_dict[-c][:n]
              tot_chg += float(c)*n
              break
            else:
              n = len(deletion_order_dict[c])
              delete_list += deletion_order_dict[-c][:n]
              tot_chg += float(c)*n
              abs_chg -= float(c)*n
          elif tot_chg > c:
            chg_to_del = c
          elif c == avail_neg_chg[-1] and not chg_to_del:
            chg_to_del = c
      if chg_to_del:
        delete_list.append(deletion_order_dict[chg_to_del].pop(0))
        tot_chg -= float(chg_to_del)
        if not deletion_order_dict[chg_to_del]:
          del deletion_order_dict[chg_to_del]
      #print (tot_chg) 
    #print (delete_list)
    #print (len(delete_list))
    return (delete_list)

  #-----------------------------------------------------------------------------
  # Function gen_deletion_order_dict
  #
  # Operation: Generates the order of molecules to be considered for deletion
  #            to force neutrality of a cluster
  #-----------------------------------------------------------------------------
  def gen_deletion_order_dict(self,random=False):
    deletion_order_dict = {}
    if random == True:
      for m in self.mol_list:
        chg_dec = Decimal(m.charge.quantize(Decimal(1.00000)))
        if chg_dec not in deletion_order_dict:
          deletion_order_dict[chg_dec] = [m]
        else:
          deletion_order_dict[chg_dec].append(m)
      for c in deletion_order_dict:
        deletion_order_dict[c] = random.shuffle(deletion_order_dict[c])
    else:
      dist_list = [m.dist_to_origin for m in self.mol_list]
      #print (dist_list) 
      indexed_dist_list = zip(range(len(dist_list)),dist_list)
      sorted_dist_list = sorted(indexed_dist_list,key=lambda x: x[1], reverse=True)
      for i in sorted_dist_list:
        m = self.mol_list[i[0]]
        chg_dec = Decimal(m.charge).quantize(Decimal(1.00000))
        if chg_dec not in deletion_order_dict:
          deletion_order_dict[chg_dec] = [m]
        else:
          deletion_order_dict[chg_dec].append(m)
    #print (deletion_order_dict) 
    return (deletion_order_dict)
  
  #-----------------------------------------------------------------------------
  # Function gen_tiered_label
  #
  # Operation: Adds the tiered_label attribute to each atom representing the
  #            region the atom is located in
  #-----------------------------------------------------------------------------
  def gen_tiered_label(self,qm_suffix,active_suffix=None,other_suffix=None,rm_pre_numbers_in_label=True):
    for i,a in enumerate(self.atom_list):
      if rm_pre_numbers_in_label == False:
        a.tiered_label = a.label
      else:
        a.tiered_label = ''.join(c for c in a.label if not c.isdigit())
      if (i+1) in self.qm_atom_list:
        a.tiered_label += qm_suffix
      elif self.active_atom_list and (i+1) in self.active_atom_list:
        a.tiered_label += active_suffix
      elif other_suffix:
        a.tiered_label += other_suffix
    return

  #-----------------------------------------------------------------------------
  # Function gen_orig_to_new_atom_no_mapping
  #
  # Operation: Writes to atom_no_mapping.txt the mapping of atom numbers of
  #            atoms in the original BulkSystem to those in the ClusterSystem
  #-----------------------------------------------------------------------------
  def gen_orig_to_new_atom_no_mapping(self,output_dir):
    with open(os.path.join(output_dir,'atom_no_mapping.txt'),'w') as f:
      title_line = 'orig_atom_no new_atom_no\n'
      f.write(title_line)
      for a in sorted(self.orig_to_new_atom_no_dict.keys()):
        line = str(a).rjust(12) + ' ' + str(self.orig_to_new_atom_no_dict[a]).rjust(11) + '\n'
        f.write(line)
    return
  
  #-----------------------------------------------------------------------------
  # Function update_atom_and_mol_region_attribute
  #
  # Operations: Assigns the region attribute for each atom and molecule
  #-----------------------------------------------------------------------------
  def update_atom_and_mol_region_attribute(self):
    for i,m in enumerate(self.mol_list):
      if (i+1) in self.qm_mol_list:
        m.region = 'QM'
        for a in m.atom_list:
          a.region = 'QM'
      elif (i+1) in self.active_mol_list:
        m.region = 'Active'
        for a in m.atom_list:
          a.region = 'Active'
      else:
        m.region = 'Frozen'
        for a in m.atom_list:
          a.region = 'Frozen'
    return

  #-----------------------------------------------------------------------------
  # Function gen_atom_by_atom_info
  #
  # Operation: Generates a file containing key info of all atoms in the system
  #-----------------------------------------------------------------------------
  def gen_atom_by_atom_info(self,output):
    with open(output,'w') as f:
      f.write('AtomNo.'.ljust(9) +\
              'Label'.ljust(10) +\
              'MolNo.'.ljust(8) +\
              'MolLabel'.ljust(12) +\
              'Charge'.ljust(9) +\
              'Dist2Or'.ljust(10) +\
              'Region\n')
      for a in self.atom_list:
        line = str(a.atom_no).ljust(9) + \
               a.label.ljust(10) + \
               str(a.parent_mol_obj.mol_no).ljust(8) + \
               str(a.parent_mol_obj.label).ljust(12) + \
               str(round(a.charge,4)).ljust(9) + \
               str(round(a.dist_to_origin,4)).ljust(10) + \
               a.region + '\n'
        f.write(line)
    return

  #-----------------------------------------------------------------------------
  # Function gen_mol_by_mol_info
  #
  # Operation: Generates a file containing key info of all molecules in the 
  #            system
  #-----------------------------------------------------------------------------
  def gen_mol_by_mol_info(self,output):
    with open(output,'w') as f:
      f.write('MolNo.'.ljust(8) +\
              'Label'.ljust(12) +\
              'AtomNos'.ljust(10) +\
              'Charge'.ljust(9) +\
              'COM2Or'.ljust(10) +\
              'MinD2Or'.ljust(10) +\
              'MaxD2Or'.ljust(10) +\
              'Region\n')
      for m in self.mol_list:
        if len(m.atom_list) > 1:
          atom_no_str = str(m.atom_list[0].atom_no) + '-' + str(m.atom_list[-1].atom_no)
        elif len(m.atom_list) == 1:
          atom_no_str = str(m.atom_list[0].atom_no) 
        line = str(m.mol_no).ljust(8) + \
               m.label.ljust(12) + \
               atom_no_str.ljust(10) + \
               str(round(m.charge,4)).ljust(9) + \
               str(round(m.dist_to_origin,4)).ljust(10) + \
               str(round(m.min_dist_to_origin,4)).ljust(10) + \
               str(round(m.max_dist_to_origin,4)).ljust(10) + \
               m.region + '\n'
        f.write(line)

    return

  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------
  #def gen_region_by_region_info(self):
  #  return
  

#===============================================================================
# Class EmbeddedClusterSystem
#
#       An EmbeddedClusterSystem object represents a cluster embedded in a shell
#       of point charges.
#===============================================================================
class EmbeddedClusterSystem(ClusterSystem):
  def __init__(self,bare_cluster_obj,radius,pt_charge_info,pt_charge_label):
    self.system_type = 'EmbeddedCluster' 
    GenSystem.__init__(self,\
                       bare_cluster_obj.mol_info,\
                       bare_cluster_obj.mol_list,\
                       bare_cluster_obj.atom_info,\
                       bare_cluster_obj.atom_list,\
                       None,\
                       bare_cluster_obj.density,\
                       bare_cluster_obj.vdw_info)
    self.radius = radius
    #self.radius_type = radius_type
    #self.dist_list = bare_cluster_obj.dist_list
    
    self.reset_atom_and_mol_no()
    self.gen_ordered_mol_list()
    self.qm_atom_list = bare_cluster_obj.qm_atom_list
    self.qm_mol_list = bare_cluster_obj.qm_mol_list 
    self.active_atom_list = bare_cluster_obj.active_atom_list
    self.active_mol_list = bare_cluster_obj.active_mol_list

    self.pt_charge_list = []
    for i,c in enumerate(pt_charge_info):
      self.pt_charge_list.append(PointCharge(pt_charge_label,c[-1],c[:3],i+1))

    return

  #-----------------------------------------------------------------------------
  # Function gen_tiered_label
  #
  # Operation: Adds a suffix to the label of each atom according to which region
  #            it is in.
  #-----------------------------------------------------------------------------
  def gen_tiered_label(self,qm_suffix,active_suffix=None,other_suffix=None,pt_charge_suffix=None,rm_pre_numbers_in_label=True):
    for i,a in enumerate(self.atom_list):
      if rm_pre_numbers_in_label == False:
        a.tiered_label = a.label
      else:
        a.tiered_label = ''.join(c for c in a.label if not c.isdigit())
      if (i+1) in self.qm_atom_list:
        a.tiered_label += qm_suffix
      elif self.active_atom_list and (i+1) in self.active_atom_list:
        a.tiered_label += active_suffix
      elif other_suffix:
        a.tiered_label += other_suffix
    for i,c in enumerate(self.pt_charge_list):
      c.tiered_label = c.label + pt_charge_suffix
    return

  #-----------------------------------------------------------------------------
  # Function write_chemshell_fragment
  #
  # Operation: Outputs a ChemShell input file to generate a fragment
  #-----------------------------------------------------------------------------
  def write_chemshell_fragment(self,chemshell_file='system.chm',pun_file='embedded_cluster.pun'):
    #chemshell_file = os.path.join(output_dir,chemshell_file)
    output_dir = os.path.dirname(chemshell_file)
    with open(chemshell_file,'w') as f:
      f.write('c_create coords=' + pun_file + ' {\n')
      if self.unit_cell_length:
        f.write('space_group\n')
        f.write('1\n')
        f.write('cell_constants angstrom\n')
        ucl_dec = Decimal(str(self.unit_cell_length)).quantize(Decimal('1.000000'))
        f.write(str(ucl_dec) + ' ' + str(ucl_dec) + ' ' + str(ucl_dec) + ' 90.0 90.0 90.0\n')
        f.write('coordinates\n')
      else:
        f.write('coordinates angstrom\n')
      for a in self.atom_list:
        atom_line = a.label + ' ' 
        for c in a.coords:
          if self.unit_cell_length:
            atom_line += str(Decimal(str(c/self.unit_cell_length)).quantize(Decimal('1.000000'))) + ' '
          else:
            atom_line += str(Decimal(str(c)).quantize(Decimal('1.000000'))) + ' '
        atom_line += 'charge ' + str(a.charge)
        f.write(atom_line + '\n')

      for p in self.pt_charge_list:
        charge_line = p.label + ' '
        for c in p.coords:
          if self.unit_cell_length:
            charge_line += str(Decimal(str(c/self.unit_cell_length)).quantize(Decimal('1.000000'))) + ' '
          else:
            charge_line += str(Decimal(str(c)).quantize(Decimal('1.000000'))) + ' '
        charge_line += 'charge ' + str(p.charge)
        f.write(charge_line + '\n')
      f.write('}\n')
    self.generate_connectivity_file(output_dir)
    
    return

  #-----------------------------------------------------------------------------
  # Function write_xyz
  #
  # Operation: Writes system configuration to an xyz file
  #-----------------------------------------------------------------------------
  def write_xyz(self,xyz_file='system.xyz',title=None,tiered=False,append=False):
    if append == True:
      f_option = 'a'
    else:
      f_option = 'w'
    with open (xyz_file,f_option) as f:
      f.write(str(len(self.atom_list)+len(self.pt_charge_list)) + '\n')
      if not title:
        f.write('Generated by MolCluster.py\n')
      else:
        f.write(title + '\n')
      for a in self.atom_list:
        f.write(a.coords_to_string(8,tiered) + '\n')
      for c in self.pt_charge_list:
        f.write(c.coords_to_string(8,tiered) + '\n')
    return
