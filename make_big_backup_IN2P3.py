#!/usr/local/bin/python

import os
import sys
import stat
import glob
import random
import tarfile
import shutil

#~ if len(sys.argv) != 2:
  #~ sys.stdout.write( 'ERROR : missing \'generation number\' operand\n' )
#~ else:
  #~ num_gener = sys.argv[1]
  
  
if 1 == 1:
  ################################################################################################################
  ## Retreive the simulation name, the number of the last saved generation and the SPS directory of the experiment
  ################################################################################################################
  # Get simulation name
  sim_name_file = open( 'sim_name.txt', 'r' )
  sim_name = (sim_name_file.readline()).rstrip( '\n' )
  sim_name_file.close()
  
  # Get number of the last saved generation
  last_gener_file = open( 'last_gener.txt', 'r' )
  last_gener = (last_gener_file.readline()).rstrip( '\n' )
  last_gener_file.close()
  
  # Get SPS dir
  sps_dir_file = open( 'SPS_dir.txt', 'r' )
  sps_dir = (sps_dir_file.readline()).rstrip( '\n' )
  sps_dir_file.close()
  
  
  
  ##################################################################################
  ## Save a copy of all the stat files in a targzipped "stat_<last_gener>" directory
  ##################################################################################
  # mkdir stat_<last_gener>
  os.makedirs( 'stat_' + last_gener )
  
  # cp stat_*.out stat_<last_gener>
  for filename in glob.glob( os.path.join( '.', 'stat_*.out' ) ):
    shutil.copy( filename, 'stat_' + last_gener + '/' + filename )
  
  # tar zcf stat_<last_gener>.tgz stat_<last_gener>
  stat_archive = tarfile.open( 'stat_' + last_gener + '.tgz', 'w:gz' )
  stat_archive.add( 'stat_' + last_gener )
  stat_archive.close()
  
  # rm -rf stat_<last_gener>
  shutil.rmtree( 'stat_' + last_gener )
  
  
  
  #######################################################################################
  ## Save the backup directory in "backup_<last_gener>.tgz" and empty it (the backup dir)
  #######################################################################################
  # tar zcf backup_<last_gener>.tgz backup
  backup_archive = tarfile.open( 'backup_' + last_gener + '.tgz', 'w:gz' )
  backup_archive.add( 'backup' )
  backup_archive.close()
  
  # rm -f backup/*
  shutil.rmtree( 'backup' )
  os.makedirs( 'backup' )
  
  
  
  #################################################################################
  ## Save the tree directory in "tree_<last_gener>.tgz" and empty it (the tree dir)
  #################################################################################
  # tar zcf tree_<last_gener>.tgz tree
  tree_archive = tarfile.open( 'tree_' + last_gener + '.tgz', 'w:gz' )
  tree_archive.add( 'tree' )
  tree_archive.close()
  
  # rm -f tree/*
  shutil.rmtree( 'tree' )
  os.makedirs( 'tree' )
  
  
  
  ####################################################################################
  ## Put these files (stat, backup and tree) into a tar archive and move it to the SPS
  ####################################################################################
  # tar zf <sim_name>_<last_gener>.tgz (stat, backup, tree)_<last_gener>.tgz
  general_archive = tarfile.open( sim_name + '_' + last_gener + '.tar', 'w' )
  general_archive.add( 'stat_' + last_gener + '.tgz' )
  general_archive.add( 'backup_' + last_gener + '.tgz' )
  general_archive.add( 'tree_' + last_gener + '.tgz' )
  # add *.txt to general archive
  for filename in glob.glob( os.path.join( '.', '*.txt' ) ):
    general_archive.add( filename )
  general_archive.close()
  
  # rm (stat, backup, tree)_<last_gener>.tgz
  os.remove( 'stat_' + last_gener + '.tgz' )
  os.remove( 'backup_' + last_gener + '.tgz' )
  os.remove( 'tree_' + last_gener + '.tgz' )
  
  # move general archive to SPS
  #os.rename( sim_name + '_' + last_gener + '.tar', sps_dir + '/' + sim_name + '_' + last_gener + '.tar')
  shutil.copy( sim_name + '_' + last_gener + '.tar', sps_dir )
  os.remove( sim_name + '_' + last_gener + '.tar' )
  