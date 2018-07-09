import os

Open(FilePath="C:/Projects/data/workbench_files/test1.wbpj")

os.chdir(r'C:\Projects\wrappers')
# Getting paramenter and write it to xml
base_design_point = Parameters.GetBaseDesignPoint()

#all_paramenters = Parameters.GetAllParameters()


param_dict = Parameters.GetInputParameterValues(base_design_point)

with open('intro.json','w') as f:
    for param, value in param_dict:
        f.write(str(param) + '=' + value + '\n')


