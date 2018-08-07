import shutil
import glob
import os


def del_files_list(files_list):
     for file in files_list:
         print('Deleting file : ' + file)
         os.remove(file)
    
def get_files(folder,exception_ext = 'txt'):
    files_list = []
    for root, dirs, files in os.walk(folder):
        for f in files:
            if f is not files_list:
                if f[-3:] != exception_ext[-3:]:
                    local_folder = os.path.join(folder, root)
                    files_list.append(os.path.join( root,f))
        #for d in dirs:
        #    files_ = get_files(os.path.join(folder,d),exception_ext)
        #    files_list.extend(files_)
    
    return files_list    
        

if __name__ == '__main__':
    import sys
    folder_path = sys.argv[1]
    print(folder_path)
    files_list =  get_files(folder_path)
    del_files_list(files_list)

