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
	
	self.working_dir = '/services/beagle/rouzaudc/phenotypic_noise/experiments/'
	self.origins = '/services/beagle/rouzaudc/phenotypic_noise/origins/'

        self.aevol_binary_directory = '/services/beagle/rouzaudc/phenotypic_noise/binary/'
	self.template_param_file = self.origins+'/param_tmpl.in'
	self.binding_matrix_file = self.working_dir+'/binding_matrix.rae'
	self.lookup_table_file = self.working_dir+'/lookup_table.ae'
	self.nb_last_generation = 100000
	
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
          'seed' : [51456165, 33263658, 7158785, 456847894, 1223144, 878944, 121145, 3587842, 2784564, 68984554],
          'origins' : [ 'WT1','WT2','WT3','WT4','WT5','WT6'],
          'env' : [ 'constant','single' ]
        }
        sweeps = sweep(parameters)
        self.sweeper = ParamSweeper(os.path.join(self.result_dir, "sweeps"), sweeps)
        logger.info('Number of parameters combinations %s', len(self.sweeper.get_remaining()))
        
    def write_run_script(self, bucketname, resume=False, resumeGeneration = 0):
        script_file = bucketname+'/run.sh'
        
        if os.path.isfile(script_file):
            os.remove(script_file)
            
        launch_cmd = 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/services/beagle/intel_lib_2019/ '+self.aevol_binary_directory+'/aevol_run -e '+str(self.nb_last_generation)+' -p 16 '
        
        launch_cmd += '&> output.log'
        

        script_file_fd = open(script_file, 'w')
        script_file_fd.write('#!/bin/bash'+os.linesep)
        script_file_fd.write(launch_cmd+os.linesep)
        script_file_fd.close()
        
        os.chmod(script_file,0777)
    
    def write_create_script(self, bucketname,chromosome_file):
        script_file = bucketname+'/run_create.sh'
        
        if os.path.isfile(script_file):
            os.remove(script_file)
            
        launch_cmd = 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/services/beagle/intel_lib_2019/ '+self.aevol_binary_directory+'/aevol_create -C '+chromosome_file+' '
            
        launch_cmd += '&> create_output.log'
        
        script_file_fd = open(script_file, 'w')
        script_file_fd.write('#!/bin/bash'+os.linesep)
        script_file_fd.write(launch_cmd+os.linesep)
        script_file_fd.close()     
        
        os.chmod(script_file,0777)
        
    def write_param_file(self, comb, bucketname):
        param_file = bucketname+'/param.in'
        
        f_template = open(self.template_param_file)
        f = open(param_file, 'w')
        
        for line in f_template:
            line = line.replace('SEED_NUMBER',str(comb['seed']))
            if comb['env'] == 'constant':
                line = line.replace('ENV_VARIATION_VALUE','switch_environment 0.0')
                line = line.replace('ENV_GAUSSIAN_VALUE','NB_ENVIRONMENTS 1\nENV_ADD_GAUSSIAN 1 1.2   0.52  0.12\nENV_ADD_GAUSSIAN 1 -1.4   0.5   0.07\nENV_ADD_GAUSSIAN 1 0.3   0.8   0.03\n')
            else:
                line = line.replace('ENV_VARIATION_VALUE','switch_environment 50000.0')
                line = line.replace('ENV_GAUSSIAN_VALUE','NB_ENVIRONMENTS 2\nENV_ADD_GAUSSIAN 1 1.2   0.52  0.12\nENV_ADD_GAUSSIAN 1 -1.4   0.5   0.07\nENV_ADD_GAUSSIAN 1 0.3   0.803   0.03\nENV_ADD_GAUSSIAN 2 1.2   0.52  0.12\nENV_ADD_GAUSSIAN 2 -1.4   0.5   0.07\nENV_ADD_GAUSSIAN 2 0.3   0.803   0.03\n')
            f.write(line)
                                                
        f_template.close()
        f.close()
        

    def workflow(self, comb):
        """ """
        bucketname = self.working_dir+'/'+slugify(comb)+'/'
        
        # If directory is empty
        if os.path.isdir(bucketname) and os.path.exists(bucketname+'/last_gener.txt'):
            gen_file = open(bucketname+'/last_gener.txt', 'r')
            last_gen = gen_file.read().strip('\n')
            gen_file.close()
            logger.info("Managing combination "+slugify(comb))

            if int(last_gen) < self.nb_last_generation:
                # Generate run.sh file
                self.write_run_script(bucketname,True,last_gen)
                
                # Launch aevol_run
                frontend = Local('cd '+bucketname+'; oarsub -l /nodes=1/core=16,walltime=120:00:00 ./run.sh -p \"cluster=\'beagle\'\" | grep OAR_JOB_ID | cut -d\'=\' -f2',process_args = { 'shell': True })
                frontend.run()
                
                a_message = ''
                for p in frontend.processes:
                    a_message = p.stdout.strip('\n')
                        
                # Fill dict with key: comb, value: oar job id
                self.oarjob_dict[comb] = a_message
                
                logger.info("Resuming R-Aevol combination "+slugify(comb)+" from "+last_gen.strip('\n')+' (OAR_JOB_ID: '+a_message+')')
            else:
                logger.info("R-Aevol combination "+slugify(comb)+" is DONE")
                # Mark comb has done
                self.sweeper.done(comb)
        else:
            frontend = Local('mkdir -p '+bucketname,process_args = { 'shell': True })
            frontend.run()
            
            # Generate param file
            self.write_param_file(comb,bucketname)
            
            # Launch aevol_create
            chromosome_file = self.origins + '/genome_' + comb['origins'] + '.txt'
            self.write_create_script(bucketname,chromosome_file)
                        
            time.sleep(10)
            
            frontend = Local('cd '+bucketname+'; oarsub -l /nodes=1/core=1,walltime=120:00:00 ./run_create.sh -p \"(cluster=\'SIC\' OR cluster=\'beagle\' OR cluster=\'kinovis\' OR cluster=\'thoth\' OR cluster=\'mistis\')\" | grep OAR_JOB_ID | cut -d\'=\' -f2',process_args = { 'shell': True })
            frontend.run()
            
            a_message = ''
            for p in frontend.processes:
                a_message = p.stdout.strip('\n')
            
            oarjobcreate_id = a_message
            
            frontend = Local('oarstat -f -j '+oarjobcreate_id+'  | grep state | cut -d\'=\' -f2  | cut -d\' \' -f2',process_args = { 'shell': True })
            frontend.run()
            
            a_message = ''
            for p in frontend.processes:
                a_message = p.stdout.strip('\n')
            
            while a_message == 'Waiting' or a_message == 'Running':
                time.sleep(30)
                frontend = Local('oarstat -f -j '+oarjobcreate_id+'  | grep state | cut -d\'=\' -f2  | cut -d\' \' -f2',process_args = { 'shell': True })
                frontend.run()
                
                a_message = ''
                for p in frontend.processes:
                    a_message = p.stdout.strip('\n')
                
                
            # Generate run.sh file
            self.write_run_script(bucketname)
                
            time.sleep(10)
            
            # Launch aevol_run
            frontend = Local('cd '+bucketname+'; oarsub -l /nodes=1/core=16,walltime=120:00:00 ./run.sh -p \"(cluster=\'SIC\' OR cluster=\'beagle\' OR cluster=\'kinovis\' OR cluster=\'thoth\' OR cluster=\'mistis\')\" | grep OAR_JOB_ID | cut -d\'=\' -f2',process_args = { 'shell': True })
            frontend.run()
            
            a_message = ''
            for p in frontend.processes:
                a_message = p.stdout.strip('\n')
                
                
            # Fill dict with key: comb, value: oar job id
            self.oarjob_dict[comb] = a_message
                
            logger.info("Start R-Aevol combination "+slugify(comb)+' (OAR_JOB_ID: '+a_message+')')


if __name__ == "__main__":
    engine = raevol_matrix()
    engine.start()
