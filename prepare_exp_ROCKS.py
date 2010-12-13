#! /usr/bin/env python3.1

import os
import random
import tarfile

# ****************************** SETTINGS ******************************
nb_seeds = 5
mutation_rates = ['1e-6', '2e-6', '5e-6', '1e-5', '5e-5', '1e-4']
# The following are needed to create the scripts for the BSMC cluster.
# If you don't need these scripts, comment the code section number 1
executable = '~/src/aevol-2.0.0/src/aevol'
base_dir = '~/sim/aevol-2.0.0/IPG_20-11-2009/'
script_name = 'IPG_20-11-2009'
nb_gener = 20000
# **********************************************************************




# **********************************************************************
# 1) Build the scripts that will be launched on the cluster
# **********************************************************************
scripts_archive = tarfile.open( script_name + '__scripts.tar.gz', 'w:gz')
for mut_rate in mutation_rates:
  for i in range(1, nb_seeds+1):
    script_file_name = script_name + '__mut_' + mut_rate + '_seed' + str(i) + '.sh'
    script_file = open(script_file_name, 'w')
    script_file.write('#!/bin/bash\n')
    script_file.write('#\n')
    script_file.write('#$ -m beas\n') # Send mail on beginning, end, abortion or suspension of job
    script_file.write('#$ -M david.parsons@insa-lyon.fr\n') # Mailto 
    script_file.write('#$ -cwd\n') # 
    script_file.write('#$ -S /bin/bash\n')
    script_file.write('#\n')
    script_file.write('\n')
    script_file.write('cd ' + base_dir + 'mut_' + mut_rate + '/seed' + str(i) + '\n')
    script_file.write(executable + ' -n ' + str(nb_gener) + ' > /dev/null\n')
    script_file.close()
    scripts_archive.add(script_file_name)
    os.remove(script_file_name)

scripts_archive.close()



# **********************************************************************
# 2) Create tree structure
# **********************************************************************
for mut_rate in mutation_rates:
  if not os.path.exists ('mut_' + mut_rate):
    os.makedirs ('mut_' + mut_rate)
  for i in range(1, nb_seeds+1):
    path = 'mut_' + mut_rate + '/seed' + str(i)
    if not os.path.exists (path):
      os.makedirs (path)
      

# **********************************************************************
# 3) Create parameter files using template
# **********************************************************************
random.seed()


for mut_rate in mutation_rates:
  for i in range(1, nb_seeds+1):
    path = 'mut_' + mut_rate + '/seed' + str(i) + '/'
    template_param_file = open('param.in', 'r')
    new_param_file = open( path + 'param.in', 'w' )

    for cur_line in template_param_file:
      if cur_line.find('SEED') == 0:
        new_param_file.write('SEED ' + str(random.randint(1, 10000000)) + '\n')
      elif cur_line.find('POINT_MUTATION_RATE') == 0:
        new_param_file.write('POINT_MUTATION_RATE   ' + mut_rate + '\n')
      elif cur_line.find('SMALL_INSERTION_RATE') == 0:
        new_param_file.write('SMALL_INSERTION_RATE  ' + mut_rate + '\n')
      elif cur_line.find('SMALL_DELETION_RATE') == 0:
        new_param_file.write('SMALL_DELETION_RATE   ' + mut_rate + '\n')
      elif cur_line.find('DUPLICATION_RATE') == 0:
        new_param_file.write('DUPLICATION_RATE      ' + mut_rate + '\n')
      elif cur_line.find('DELETION_RATE') == 0:
        new_param_file.write('DELETION_RATE         ' + mut_rate + '\n')
      elif cur_line.find('TRANSLOCATION_RATE') == 0:
        new_param_file.write('TRANSLOCATION_RATE    ' + mut_rate + '\n')
      elif cur_line.find('INVERSION_RATE') == 0:
        new_param_file.write('INVERSION_RATE        ' + mut_rate + '\n')
      else:
        new_param_file.write(cur_line)
        
    new_param_file.close()
    template_param_file.close()
    
    
    
    

# **********************************************************************
# 4) Put the freshly created tree structure in a tarball (and delete them)
# **********************************************************************
dirs_archive = tarfile.open( script_name + '__dirs.tar.gz', 'w:gz')
for mut_rate in mutation_rates:
  dirs_archive.add('mut_' + mut_rate)
dirs_archive.close()


for mut_rate in mutation_rates:
  for i in range(1, nb_seeds+1):
    os.remove('mut_' + mut_rate + '/seed' + str(i) + '/param.in')
    os.rmdir('mut_' + mut_rate + '/seed' + str(i))
  os.rmdir('mut_' + mut_rate)


    