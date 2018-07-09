# -*- coding: utf-8 -*-
"""
@author: jenovencio
"""

import subprocess
import os

class AnsysWorkbench():
    ''' This class can open Ansys Workbench,
        do some introspection to see input and output variables
        change variables, run a work bench case, and also send some 
        APDL commands
    ''' 
    
    def __init__(self, workbench_wbpj = '', working_directory = '', workbench_path= u'C:\Program Files\ANSYS Inc\v181'):
        
        self.workbench_wbpj = workbench_wbpj
        self.workbench_path = workbench_path
        self.working_directory = working_directory
        self.workbench_exe = os.path.join(self.workbench_path,u'\Framework\bin\Win64\RunWB2.exe')
        
    def set_workbench_wbpj(self, workbench_wbpj):
        ''' This function set the workbench project
            
        arguments
            workbench_path : str
                string with workbench execuble path
        
        returns
            a new string with workbench executable path
        
        '''
        self.workbench_wbpj = workbench_wbpj
        print(self.workbench_wbpj)
        return self.workbench_wbpj 
        
    def set_workbench_path(self, workbench_path):
        ''' This function set the workbench executable path
        
        arguments
            workbench_path : str
                string with workbench execuble path
        
        returns
            a new string with workbench executable path
        '''
        self.workbench_path = workbench_path
        print(self.workbench_path)
        return self.workbench_path
        
    def set_working_directory(self, working_directory):
        ''' This function sets the the working directory where Ansys will be 
        run
        
        arguments
            working_directory : str
                string with working_directory
        
        returns
            a new string with wworking_directory
        '''
        self.working_directory = working_directory
        print(self.working_directory)
        return self.working_directory
        
    def openGUI(self):
        ''' This function open a workbench project
        '''
        
        try:
            if self.workbench_wbpj: 
                command = self.workbench_exe + '-F '  + self.workbench_wbpj 
            else:
                command = self.workbench_exe
            print(command)
            subprocess.call(command, shell=True)
        except:    
            print("Please make sure that Ansys Workbench path was set correctly")
        
    def intro(self):
        pass 
        
if '__name__' == '__main__':
    
    workbench_wbpj = ''
    workbench_path = u'C:\Program Files\ANSYS Inc\v181'
    working_directory = ''
    
    wb_case = AnsysWorkbench(workbench_wbpj, working_directory, workbench_path)
    wb_case.openGUI()