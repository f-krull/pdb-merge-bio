# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 18:13:35 2012

@author: tmeyer
"""

import os

import sys
#sys.path.append('/user/tmeyer/workspace/script/python/bivalent_ligands')


def merge_biounits(pdb_file_list):

    import pdb_analysis

    errors = ''

    MergeError  = pdb_analysis.pdb_from_biopython.MergeError
    IO_Error    = pdb_analysis.pdb_from_biopython.IO_Error
    HeaderError = pdb_analysis.pdb_from_biopython.HeaderError

    for i, pdb_file_names in enumerate(pdb_file_list):

        (pdb_file, source_filename, target_filename) = pdb_file_names

        if os.path.exists(target_filename):
            #print 'skipping: ' + pdb_file
            continue
        else:
            #print '# processing ' + str(i) + ' of ' + str(num_of_pdbs) + ' : ' + pdb_file
            print 'processing: ' + pdb_file

        filesize = os.path.getsize(source_filename)
        if filesize / 1e6 > 20:
            # File is too large.
            errors += 'filename: ' + source_filename + '\n'
            errors += 'message: "' + 'This structure seems to be too large: '
            errors += str(filesize / 1e6) + ' MB'  + '"\n\n'
            continue

        pdb = pdb_analysis.pdb_from_biopython()

        try:

            # Merge models in biological assambly.
            pdb.parse_pdb_file(str(source_filename), quiet=True)

            print 'subunits: ' + str(len(pdb.all_structs))

#            if len(pdb.all_structs) > 100:
#                errors += 'filename: ' + source_filename + '\n'
#                errors += 'message: "' + 'This structure seems to be too large: ' + str(len(pdb.all_structs)) + ' subunits!'  + '"\n\n'
#                continue

            #if not pdb.is_nmr:
            pdb.merge_models(overwrite=False)

            pdb.write_pdb_file(str(target_filename))

        except (ValueError, IndexError, IO_Error, HeaderError, MergeError) as inst:
            #f_err.write('### ERROR: pdb-file could not be parsed.\n')
            errors += 'filename: ' + source_filename + '\n'
            errors += 'message: "' + inst.__str__() + '"\n\n'
            continue
        except:
            errors += 'filename: ' + source_filename + '\n'
            e = sys.exc_info()[1]
            e_str = "Error: %s" % e
            errors += 'message: "' + e_str + '"\n\n'
            continue


    return str(errors)





if __name__ == '__main__':

    source_folder = './pdb_bio/'
    target_folder = './pdb_bio_merged/'

    # This file contains all error messages that occured during parsing and
    # merging the biological assamblies.
    f_err = open(target_folder + '/errors.log', 'w')



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
                ### for pdb folder
                # pdb_code = pdb_file[3:-7]
                # assert len(pdb_code) == 4
                # target_filename = current_target_subfolder + '/' + pdb_code + '.pdb0'

    #            if pdb_file.find('.ent') > 0:
    #                print source_filename
    #                continue

                if not os.path.exists(target_filename):
                    pdb_file_list.append( (pdb_file, source_filename, target_filename) )


    num_of_pdbs = len(pdb_file_list)
    print 'I found ' + str(num_of_pdbs) + ' pdb files.'


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


    ### submitt jobs ###
    #jobs_submitted = []
    files_left = num_of_pdbs
    for file_chunk in files_splitted:
        errors = merge_biounits(file_chunk)
        if errors != '':
            f_err.write(errors + '\n\n')
        f_err.flush()
        files_left -= max_files

    f_err.close()

