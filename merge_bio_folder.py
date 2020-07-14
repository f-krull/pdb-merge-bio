# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 18:13:35 2012

@author: tmeyer
"""

#------------------------------------------------------------------------------

import os
import sys
from multiprocessing import Pool, Value, Lock, Array

#------------------------------------------------------------------------------

def merge_biounits(pdb_file_list):
    import pdb_analysis

    myjobid = 0
    with lock:
        myjobid = jobid.value
        jobid.value += 1
        
    errors = ''
    MergeError  = pdb_analysis.pdb_from_biopython.MergeError
    IO_Error    = pdb_analysis.pdb_from_biopython.IO_Error
    HeaderError = pdb_analysis.pdb_from_biopython.HeaderError

    for i, pdb_file_names in enumerate(pdb_file_list):
        (pdb_file, source_filename, target_filename) = pdb_file_names
        if os.path.exists(target_filename):
            with lock:
                print('skipping: ' + pdb_file)
                sys.stdout.flush()
            continue
        else:
            with lock:
                print(str(myjobid) + ' processing: ' + pdb_file)
                sys.stdout.flush()

        filesize = os.path.getsize(source_filename)
        # filesize*100 = mem-usage -> 3M PDB-gz needs about 300M memory
        if (filesize * 100) > g_threadMaxmem_kb * 1e3:
            errors += 'error ' + source_filename + ' '
            errors += 'This structure seems to be too large ('
            errors += str(filesize / 1e3) + ' kB'  + ')\n'
            continue

        pdb = pdb_analysis.pdb_from_biopython()
        try:
            # Merge models in biological assambly.
            pdb.parse_pdb_file(str(source_filename), quiet=True)
            pdb.merge_models(overwrite=False)
            pdb.write_pdb_file(str(target_filename) + "_tmp")
            os.rename(str(target_filename) + "_tmp", str(target_filename))
        except (ValueError, IndexError, IO_Error, HeaderError, MergeError) as inst:
            #f_err.write('### ERROR: pdb-file could not be parsed.\n')
            errors += 'error ' + source_filename + ' '
            errors += 'message: "' + inst.__str__() + '" \n'
            continue
        except:
            errors += 'error ' + source_filename + ' '
            e = sys.exc_info()[1]
            e_str = "Error: %s" % e
            errors += 'message: ' + e_str + '\n'
            continue
    # write any error messages to log-file
    if errors != "":
        with lock:
            f_err = open(logFilename[:], 'a')
            f_err.write(errors)
            f_err.flush()
            f_err.close()


logFilename = Array('u','merge_errors.log')
jobid = Value('i', 0)
lock = Lock()
g_threadMaxmem_kb = 0

#------------------------------------------------------------------------------

def initializer(*args):
    global jobid, lock, logFilename
    jobid, lock, logFilename = args
    
#------------------------------------------------------------------------------



import argparse



if __name__ == '__main__':

    param_numthreads = 1
    param_dotest = False


    parser = argparse.ArgumentParser()
    parser.add_argument("--numthreads", help="specify the number of threads to use")
    parser.add_argument("--maxmem", help="specify the maximum amount of memory to use (kB)")
    parser.add_argument("--test", action="store_true", help="only process the first 10 PDBs for testing")
    parser.add_argument("--src", help="path to pdb_bio folder")
    parser.add_argument("--dst", help="path to output folder")
    
    args = parser.parse_args()
    
    if args.numthreads:
        param_numthreads = int(args.numthreads)
    else:
        param_numthreads = 1

    if args.maxmem:
        g_threadMaxmem_kb = int(args.maxmem) / param_numthreads
    else:
        # filesize*100 = mem-usage -> 20M PDB-gz needs about 2000M memory
        g_threadMaxmem_kb = 20 * 1e3;

    if args.test:
        param_dotest = True
        
        
        
    
    print("Using " + str(param_numthreads) + " threads")
    print("Using " + str(g_threadMaxmem_kb * param_numthreads) + " kB of memory")
    
    

    source_folder = args.src
    target_folder = args.dst

    # This file contains all error messages that occured during parsing and
    # merging the biological assamblies.

    subfolder = os.listdir(source_folder)

    pdb_file_list = []
    for sd in subfolder:
        current_source_subfolder = source_folder + '/' + sd
        current_target_subfolder = target_folder + '/' + sd
        if os.path.isdir(current_source_subfolder) and len(sd) == 2:

            # Get all pdbs in the folder.
            pdb_files = os.listdir(current_source_subfolder)

            # Create target folder exist.
            if not os.path.isdir(current_target_subfolder):
                os.mkdir(current_target_subfolder)

            # Create a list of pdb files.
            for pdb_file in pdb_files:
                source_filename = current_source_subfolder + '/' + pdb_file
                ### for bio folder
                target_filename = current_target_subfolder + '/' + pdb_file[:-3]
                # skip existing files
                if not os.path.exists(target_filename):
                    pdb_file_list.append( (pdb_file, source_filename, target_filename) )

    num_of_pdbs = len(pdb_file_list)
    print('Found ' + str(num_of_pdbs) + ' pdb files.')

    count = 0
    max_files = 1
    files_splitted = []
    current_files = []
    for pdb_file_names in pdb_file_list:
        current_files.append(pdb_file_names)
        count += 1
        if count >= max_files:
            count = 0
            files_splitted.append(current_files)
            current_files = []
    # if test: use only a small sample 
    if param_dotest:
        del files_splitted[10:]
    
    pool = Pool(param_numthreads, initializer, (jobid, lock, logFilename)) 
    # clear log file
    f_err = open(logFilename[:], 'w')
    f_err.write("")
    f_err.flush()
    f_err.close()
    # start threads
    pool.map(merge_biounits, files_splitted)
    

