#!/usr/bin/env python
# encoding=utf8

import os, time, datetime
import shutil
from shutil import copyfile

from execo import Local
from execo.log import style
from execo import logger as ex_log
from execo.time_utils import timedelta_to_seconds
from execo_engine import Engine, logger, ParamSweeper, sweep, slugify
import pickle

class raevol_matrix(Engine):

    def __init__(self):
        """Overloading class initialization with parent and adding options"""
        super(raevol_matrix, self).__init__()
        self.options_parser.set_usage("usage: %prog ")
        self.options_parser.set_description("Execo Engine that can be used to" + \
                "perform automatic R-Aevol campaign")
        self.options_parser.add_option("-w", dest="walltime",
                    help="walltime for the reservation",
                    type="string",
                    default="120:00:00")
        self.options_parser.add_option("-u", dest="selected_cluster",
                    help="run on a specific cluster.",
                    type="string",
                    default="thoth")

    def run(self):
        """ """
        self.define_parameters()
	
	self.working_dir = '/services/beagle/rouzaudc/large_pop_ltisee/experiments/'
	self.aevol_binary_directory = '/services/beagle/rouzaudc/large_pop_ltisee/binary/'
	self.template_param_file = self.working_dir+'/param_tmpl.in'
	self.binding_matrix_file = self.working_dir+'/binding_matrix.rae'
	self.lookup_table_file = self.working_dir+'/lookup_table.ae'
	self.nb_last_generation = 1000000
	
	self.cluster = self.options.selected_cluster
	
        self.oarjob_dict = {}
        self.oarjob_dict_file = self.working_dir + '/dict_backup_beagle.p'
        
        # If dict backup exist, restore it
        if os.path.isfile(self.oarjob_dict_file):
            self.oarjob_dict = pickle.load(open(self.oarjob_dict_file,'rb'))
            
        while len(self.sweeper.get_remaining()) > 0 or len(self.sweeper.get_inprogress()) > 0:
            comb = self.sweeper.get_next()
            if not comb:
                combToRelaunch = 0
                
                # Check if some experiment reservation has finished
                for kcomb in self.oarjob_dict.keys():

                    print "checking "+slugify(kcomb)+ " oar "+self.oarjob_dict[kcomb]
                    
                    frontend = Local('oarstat -f -j '+self.oarjob_dict[kcomb]+'  | grep state | cut -d\'=\' -f2  | cut -d\' \' -f2',process_args = { 'shell': True })
                    frontend.run()
                    
                    a_message = ''
                    for p in frontend.processes:
                        a_message = p.stdout.strip('\n')
                        
                    # To do it, check if some comb has a oar job id that is not Waiting or Running and relaunch workflow for this combination
                    if a_message != 'Waiting' and a_message != 'Running':
                        combToRelaunch += 1
                        self.sweeper.cancel(kcomb)
                        self.oarjob_dict.pop(kcomb, None)
                        
                # Checkpoint dict
                pickle.dump(self.oarjob_dict, open(self.oarjob_dict_file,'wb'))
                
                # If no comb to relaunch , sleep for 10 min
                if combToRelaunch == 0:
                    logger.info("Parameter sweeper : remaining "+str(len(self.sweeper.get_remaining()))+" -- in progress "+str(len(self.sweeper.get_inprogress()))
                                +" -- done "+str(len(self.sweeper.get_done())))
                    logger.info("[-----------------------------] Combination in progress [-----------------------------]")
                    
                    for kcomb in self.oarjob_dict:
                        bucketname = self.working_dir+'/'+slugify(kcomb)+'/'
                        gen_file = open(bucketname+'/last_gener.txt', 'r')
                        last_gen = gen_file.read().strip('\n')
                        gen_file.close()
                        logger.info("Combination "+slugify(kcomb)+" is at generation "+last_gen+"/"+str(self.nb_last_generation))
                    
                    logger.info("[-------------------------------------------------------------------------------------]")
                    
                    time.sleep(600)
            else:
                # Launch workflow for this combination
                self.workflow(comb)


    def define_parameters(self):
        """ """
        parameters = {
          'seed' : [51456165, 33263658, 7158785, 456847894, 1223144, 878944, 121145, 3587842, 2784564, 68984554,45564564,8789789,213145016,455445121,564691587],
          'origins' : [ 'wild_type_1','wild_type_2','wild_type_3'],
          'grid_size_height' : [ 16, 64 ],
          'grid_size_width' : [ 16, 64 ],
        }

        sweeps = sweep(parameters)
        self.sweeper = ParamSweeper(os.path.join(self.result_dir, "sweeps"), sweeps)
        logger.info('Number of parameters combinations %s', len(self.sweeper.get_remaining()))
        
    def write_run_script(self, bucketname, resume=False, resumeGeneration = 0):
        script_file = bucketname+'/run.sh'
        
        if os.path.isfile(script_file):
            os.remove(script_file)
            
        launch_cmd = 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/services/beagle/intel_lib_2019/ '+self.aevol_binary_directory+'/pt.sh '+bucketname
        
        launch_cmd += '&> output.log'
        

        script_file_fd = open(script_file, 'w')
        script_file_fd.write('#!/bin/bash'+os.linesep)
        script_file_fd.write(launch_cmd+os.linesep)
        script_file_fd.close()
        
        os.chmod(script_file,0777)
    
    def workflow(self, comb):
        """ """
        if comb['grid_size_width'] != comb['grid_size_height']:
            self.sweeper.done(comb)
            return
       
        bucketname = self.working_dir+'/'+slugify(comb)+'/'
        
        # If directory is empty
        if not os.path.exists(bucketname+'/pt_has_run'):
            # Generate run.sh file
            self.write_run_script(bucketname,True,0)
                
            # Launch aevol_run
            frontend = Local('cd '+bucketname+'; oarsub -l /nodes=1/core=8,walltime=120:00:00 ./run.sh -p \"(cluster=\'SIC\' OR cluster=\'beagle\' OR cluster=\'kinovis\' OR cluster=\'thoth\' OR cluster=\'mistis\')\" | grep OAR_JOB_ID | cut -d\'=\' -f2',process_args = { 'shell': True })
            frontend.run()
                
            a_message = ''
            for p in frontend.processes:
                a_message = p.stdout.strip('\n')
                        
            # Fill dict with key: comb, value: oar job id
            self.oarjob_dict[comb] = a_message
                
            logger.info("Launching post treatment for  "+slugify(comb)+" from (OAR_JOB_ID: "+a_message+')')
        else:
            logger.info("R-Aevol combination "+slugify(comb)+" is DONE")
            # Mark comb has done
            self.sweeper.done(comb)

if __name__ == "__main__":
    engine = raevol_matrix()
    engine.start()
