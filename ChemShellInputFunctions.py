import decimal

def gen_atom_list_str(atom_list):
  atom_list_str = '{ '
  for a in atom_list:
    atom_list_str += str(a) + ' '
  atom_list_str += '}'
  return (atom_list_str)

def gen_active_str(active_list):
  active_str = 'active_atoms = ' + gen_atom_list_str(active_list)

  return (active_str)

def gen_qm_atoms_str(qm_list):
  qm_atoms_str = 'qm_region = ' + gen_atom_list_str(qm_list)
  return (qm_atoms_str)

def gen_groups_str(n_atoms,n_pt_charges=0):
  groups_str = 'groups = { { '
  for i in range(n_atoms+n_pt_charges):
    groups_str += str(i+1) + ' '
  groups_str += '} }'
  return (groups_str)

def gen_atom_types_str(atom_list,pt_charge_list = None):
  atom_types_str = 'atom_types = { '
  for a in atom_list:
    atom_types_str += a.label + ' '
  if pt_charge_list:
    for p in pt_charge_list:
      atom_types_str += p.label + ' '
  atom_types_str += '}'
  return (atom_types_str)

def gen_constraints_str(system_obj,exclude_list=None,log_file=None):
  atom_index = 0
  constraints_str = '{ '
  constraints_flag = False
  #print ('VC: exclude_list: ' + str(exclude_list))
  if log_file:
    with open(log_file,'a') as f:
      f.write('-->Determining constraints in the system - constraints removed between atoms in the QM region.\n')
  for m in system_obj.mol_list:
    if 'constraints' in m.potential:
      if constraints_flag == False:
        constraints_flag = True
      for c in m.potential['constraints']:
        if not exclude_list:
          constraints_str += '{bond ' + str(int(c['atom_1']) + atom_index) + ' ' + str(int(c['atom_2']) + atom_index) + '} '
        else:
          a1_index = int(c['atom_1']) + atom_index
          a2_index = int(c['atom_2']) + atom_index
          if a1_index in exclude_list and a2_index in exclude_list:
            if log_file:
              with open(log_file,'a') as f:
                f.write('---->Atoms ' + str(a1_index) + '(' + system_obj.atom_list[a1_index-1].label + ') and ' + str(a2_index) + '(' + system_obj.atom_list[a2_index-1].label + ') in QM region, constraint removed.\n')
            #pass
          elif a1_index not in exclude_list and a2_index not in exclude_list:
            constraints_str += '{bond ' + str(a1_index) + ' ' + str(a2_index) + '} '
            if constraints_flag == False:
              constraints_flag = True
          elif a1_index not in exclude_list and a2_index in exclude_list:
            print ('Atom ' + str(a1_index) + ' not in the QM region but atom ' + str(a2_index) + ' is in the QM region even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!')
            if log_file:
              with open(log_file,'a') as f:
                f.write('---->Atom ' + str(a1_index) + ' not in the QM region but atom ' + str(a2_index) + ' is in the QM region even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!\n')
          elif a1_index in exclude_list and a2_index not in exclude_list:
            print ('Atom ' + str(a1_index) + ' not in the QM region but atom ' + str(a2_index) + ' is in the QM region even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!')
            if log_file:
              with open(log_file,'a') as f:
                f.write('---->Atom ' + str(a1_index) + ' not in the QM region but atom ' + str(a2_index) + ' is in the QM region even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!\n')
    atom_index += len(m.atom_list)
  constraints_str += '}'
  if constraints_flag == False:
    constraints_str = None
  return (constraints_str)

def gen_residues_str(system_obj,active_list):
  residues_str = 'residues= {{'
  if active_list:
    frozen_res_str = '{'
  else:
    frozen_res_str = ''
   
  single_res_str_list = []
  for m in system_obj.mol_list:
    if active_list:
      if m.atom_list[0].atom_no in active_list:
        single_res_str = '{'
        for a in m.atom_list:
          single_res_str += str(a.atom_no) + ' '
          if a.atom_no not in active_list:
            print ('-->Atom ' + str(a.atom_no) + ' not in the active list but atom ' + str(m.atom_list[0].atom_no) + ' is even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!\n')
        single_res_str = single_res_str[:-1] + '}'
        single_res_str_list.append(single_res_str)
      else:
        for a in m.atom_list:
          if a.atom_no in active_list:
            print ('-->Atom ' + str(a.atom_no) + ' is in the active list but atom ' + str(m.atom_list[0].atom_no) + ' is not even though they both belong to molecule ' + str(m.mol_no) + ': must be a BUG!\n')
          frozen_res_str += str(a.atom_no) + ' '
    else:
      single_res_str = '{'
      for a in m.atom_list:
        single_res_str += str(a.atom_no) + ' '
      single_res_str = single_res_str[:-1] + '}'
      single_res_str_list.append(single_res_str)
  
  if frozen_res_str != '{':
    frozen_res_str = frozen_res_str[:-1] + '}'
    for i in range(len(single_res_str_list) + 1):
      residues_str += 'res' + str(i+1) + ' '
    residues_str = residues_str[:-1] + '} {' 
    for r in single_res_str_list:
      residues_str += r + ' ' 
    residues_str += frozen_res_str + '}}'
  else:
    for i in range(len(single_res_str_list)): 
      residues_str += 'res' + str(i+1) + ' '
    residues_str = residues_str[:-1] + '} {'
    for r in single_res_str_list:
      residues_str += r + ' ' 
    residues_str = residues_str[:-1] + '}}'
  #print (residues_str)
  return (residues_str)


def qmmm_opt_input(system_obj,input_file,pun_file='cluster.pun',final_opt_pun_file='cluster_opt.pun',charge=0,mult=1,qm_constraints=False,log_file=None,coord_system='C'):
  n_atoms = len(system_obj.atom_list)
  qm_atoms_str = gen_qm_atoms_str(system_obj.qm_atom_list)
  if system_obj.system_type == 'EmbeddedCluster':
    atom_types_str = gen_atom_types_str(system_obj.atom_list,system_obj.pt_charge_list)
  else:
    atom_types_str = gen_atom_types_str(system_obj.atom_list)
  if 'pt_charge_list' in system_obj.__dict__:
    groups_str = gen_groups_str(n_atoms,len(system_obj.pt_charge_list))
  else:
    groups_str = gen_groups_str(n_atoms,0)
  if coord_system in ('D','H'):
    if qm_constraints == True:
      #print ('VC: ' + str(qm_list))
      constraints_str = gen_constraints_str(system_obj,system_obj.qm_atom_list,log_file)
    elif qm_constraints == False:
      constraints_str = gen_constraints_str(system_obj)
  else:
    constraints_str = None

  if system_obj.active_atom_list:
    active_str = gen_active_str(system_obj.active_atom_list)
  residues_str = gen_residues_str(system_obj,system_obj.active_atom_list)
  with open(input_file,'w') as f:
    f.write('dl-find coords=' + pun_file + '\\\n')
    if constraints_str:
      f.write('      constraints = ' + constraints_str + '\\\n')
    if coord_system == 'H':
      f.write('      coordinates=hdlc \\\n')
      f.write('      ' + residues_str + '\\\n')
    elif coord_system == 'D':
      f.write('      coordinates=dlc \\\n')
    f.write('      theory=hybrid: { coupling=shift \\\n')
    f.write('      cutoff = 20 \\\n')
    f.write('      ' + qm_atoms_str + '\\\n')
    f.write('      qm_theory=gaussian : { nproc=7 maxcyc=200 scfconv=5 basis=631gdp g98_mem=80000000 charge=' + str(int(decimal.Decimal(str(charge)).quantize(decimal.Decimal('0.01')+decimal.Decimal('0.01')))) + ' mult=' + str(mult) + ' hamiltonian=b3lyp } \\\n')
    f.write('      mm_theory=dl_poly : { mm_defs=ff.dat \\\n')
    f.write('      ' + atom_types_str + ' } \\\n')
    f.write('      ' + groups_str + ' } \\\n')
    if system_obj.active_atom_list:
      f.write('      ' + active_str + '\\\n')
    f.write('      list_option = full \\\n')
    f.write('      maxcycle = 600 \\\n')
    f.write('      dump = 1 \\\n')
    f.write('      result = ' + final_opt_pun_file + '\n')
    f.write('\n')
    #f.write('write_xyz coords=' + coord_file_base + '_opt.pun file=' + coord_file_base + '_opt.xyz')
  return
