import re

basename_p = re.compile('^[A-Za-z0-9_]*?$')
dash_line = '-' * 80 
double_dash_line = '=' * 80
ANGSTROM_STR = u'\u212B'.encode('utf-8')
alphanumeric_p = re.compile('^[A-Za-z0-9]*?$')
exit_line = '\n' + 'E'*80 + '\nEEE  ' + 'Exiting molcluster'.center(70) + '  EEE\n' + 'E'*80 + '\n'
#bugs_line = 'Please report any bugs to Vincent Chen (vhc08@ic.ac.uk) or Gabriel Lau (gvl07@ic.ac.uk). The log file can be found in \'' + str(log_file) + '\'.'
