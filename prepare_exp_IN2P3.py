#! /usr/bin/env python3.1

import os
import stat
import random
import tarfile

# ****************************** SETTINGS ******************************
nb_seeds = 5
mutation_rates = ['2e-4', '1e-4', '5e-5', '2e-5', '1e-5', '5e-6', '2e-6', '1e-6']
# The following are needed to create the scripts for the IN2P3 computing center.
# If you don't need these scripts, comment the code section number 1
aevol_dir     = '/sps/liris/aevol/aevol-2.2.0'
aevol_bin     = aevol_dir + '/src/aevol_IN2P3'
script_backup = aevol_dir + '/make_big_backup_IN2P3.py'
base_dir      = '/sps/liris/dparsons/sim/neutral_zones_2'
script_name   = 'NZ_2'
nb_gener      = 50000
# **********************************************************************




# **********************************************************************
# 1) Build the scripts that will be launched on the cluster
# **********************************************************************
scripts_archive = tarfile.open( script_name + '__scripts.tar.gz', 'w:gz')
submit_all_sh   = open( 'submit_all.sh', 'w' )
for mut_rate in mutation_rates:
  for i in range(1, nb_seeds+1):
    sim_dir           = 'mut_' + mut_rate + '/seed' + str(i)
    sps_sim_dir       = base_dir + '/' + sim_dir
    sim_name          = 'mut_' + mut_rate + '_seed' + str(i)
    job_name          = script_name + '__' + sim_name
    script_file_name  = 'mut_' + mut_rate + '_seed' + str(i) + '.sh'
    
    script_file = open(script_file_name, 'w')
    script_file.write('#!/bin/sh\n')
    script_file.write('#PBS -q T                             # Execution Class\n')
    script_file.write('#PBS -l platform=LINUX                # Execution Platform\n')
    script_file.write('#PBS -l u_sps_liris                   # We use the Semi Permanent Storage\n')
    script_file.write('#PBS -N ' + job_name + ' # Job name\n')
    script_file.write('#PBS -l T=4286000                     # Duration (in normalized units)\n')
    script_file.write('#PBS -l M=4096MB                      # Virtual Memory in MB\n')
    script_file.write('#PBS -l scratch=30GB                  # Scratch size in MB (disk space on calculator)\n')
    script_file.write('#PBS -l spool=500KB                   # Spool size in KB (space for stdout and stderr)\n')
    script_file.write('#PBS -l model=Xeon                    # Processor\n')
    script_file.write('#PBS -eo                              # Redirect stderr to stdout\n')
    script_file.write('#PBS -mb                              # Send Mail on Begin\n')
    script_file.write('#PBS -me                              # Send Mail on End\n')
    script_file.write('\n\n\n')
    script_file.write('# Simulation directories -> Where you have you input files and where you want your output files to be copied\n')
    script_file.write('SPS_BASE_DIR="' + base_dir + '"\n')
    script_file.write('SIM_DIR="' + sim_dir + '"\n')
    script_file.write('SIM_NAME="' + sim_name + '"\n')
    script_file.write('SPS_SIM_DIR="' + sps_sim_dir + '"\n')
    script_file.write('\n')
    script_file.write('# Binary program path.\n')
    script_file.write('EXEC="' + aevol_bin + '"\n')
    script_file.write('\n')
    script_file.write('# Create the simulation dir and cd into it\n')
    script_file.write('mkdir -p $SIM_DIR             # Create the simulation subdirectory on the calculator\n')
    script_file.write('cd $SIM_DIR                   # cd on the calculator\n')
    script_file.write('\n')
    script_file.write('# Copy the param.in file\n')
    script_file.write('cp $SPS_SIM_DIR"/param.in" .  # Copy param.in from SPS to the calculator\n')
    script_file.write('\n')
    script_file.write('# Write the sim_name.txt, SPS_dir.txt and cpu_info.txt files\n')
    script_file.write('echo ' + sim_name + ' > sim_name.txt\n')
    script_file.write('echo ' + sps_sim_dir + ' > SPS_dir.txt\n')
    script_file.write('cp /proc/cpuinfo cpu_info.txt\n')
    script_file.write('\n')
    script_file.write('# Create a symbolic link to the big_backup script\n')
    script_file.write('ln -s ' + script_backup + ' make_big_backup.py\n')
    script_file.write('\n')
    script_file.write('# Run the program "$EXEC" with some options\n')
    script_file.write('$EXEC -n ' + str(nb_gener) + ' > /dev/null\n')
    script_file.write('\n')
    script_file.write('# Make a .tar.gz from the result (with the directory structure $SIM_DIR)\n')
    script_file.write('cd ../..\n')
    script_file.write('tar zcf $SIM_NAME.tgz $SIM_DIR --remove-files\n')
    script_file.write('mv  $SIM_NAME.tgz $SPS_BASE_DIR/\n')
    script_file.close()
    
    submit_all_sh.write('qsub ' + script_file_name + '\n')
    
    os.chmod(script_file_name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
    scripts_archive.add(script_file_name)
    os.remove(script_file_name)

submit_all_sh.close()
os.chmod('submit_all.sh', stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
scripts_archive.add('submit_all.sh')
os.remove('submit_all.sh')

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
    
    
    
    

# ************************************************************************
# 4) Put the freshly created tree structure in a tarball (and delete them)
# ************************************************************************
dirs_archive = tarfile.open( script_name + '__dirs.tar.gz', 'w:gz')
for mut_rate in mutation_rates:
  dirs_archive.add('mut_' + mut_rate)
dirs_archive.close()


for mut_rate in mutation_rates:
  for i in range(1, nb_seeds+1):
    os.remove('mut_' + mut_rate + '/seed' + str(i) + '/param.in')
    os.rmdir('mut_' + mut_rate + '/seed' + str(i))
  os.rmdir('mut_' + mut_rate)


    