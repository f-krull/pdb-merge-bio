# -*- coding: utf-8 -*-

import os.path
import sys
sys.path.append('/user/tmeyer/git/')
sys.path.append('/user/tmeyer/workspace/script/python/bivalent_ligands')
sys.path.append('/scratch/scratch/tmeyer/lib')

#from kniopython.PDB.XStructure import XStructure
from nem  import CleverFile


#import pyximport; pyximport.install()
from all_vs_all_atoms import find_collisions


import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select

from Bio.PDB.StructureBuilder import StructureBuilder

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

import re

def is_std_aa(name):
    std_aa = {'ARG' : 'R', \
              'HIS' : 'H', \
              'LYS' : 'K', \
              'ASP' : 'D', \
              'GLU' : 'E', \
              'SER' : 'S', \
              'THR' : 'T', \
              'ASN' : 'N', \
              'GLN' : 'Q', \
              'CYS' : 'C', \
              'SEC' : 'U', \
              'GLY' : 'G', \
              'PRO' : 'P', \
              'ALA' : 'A', \
              'ILE' : 'I', \
              'LEU' : 'L', \
              'MET' : 'M', \
              'PHE' : 'F', \
              'TRP' : 'W', \
              'TYR' : 'Y', \
              'VAL' : 'V'}
    
    if std_aa.has_key(name) > 0:
        return 1
    else:
        return 0
      
        
class Warnings():
    def __init__(self):
        self.messages = []
        self.old_messages = []
        
    def warn(self, value):
        self.messages.append( str(value) )
        
    def has_warnings(self):
        if len(self.messages) > 0:
            return 1
        else:
            return 0
        
    def get_warnings(self):
        for wrn in self.messages:
            self.old_messages.append( wrn )
        
        w = list(self.messages)
        self.messages = []
        return w
    


class pdb_from_biopython(object):
    
    class HeaderError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)  

    class MergeError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)  

    class IO_Error(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)  
            
    
    def __init__(self):        
        self.c_file_parsed = 0
        
        self.all_structs        = None
        self.all_structs_merged = None
        
        self.warnings = Warnings()
        
        self.too_many_chains = False
        
        
    # This function overwrites the default function in the PDBParser class
    def _handle_PDB_exception(self, message, line_counter):
        """
        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE), or raises it again, this time adding the 
        PDB line number to the error message.
        """
        message="%s at line %i." % (message, line_counter)
#        if self.PERMISSIVE:
#            # just print a warning - some residues/atoms may be missing
#            warnings.warn("PDBConstructionException: %s\n"
#                          "Exception ignored.\n"
#                          "Some atoms or residues may be missing in the data structure."
#                          % message, PDBConstructionWarning)
#        else:
#            # exceptions are fatal - raise again with new message (including line nr)
#            raise PDBConstructionException(message)



    def __analyse_header(self):
        # Analyse LINK entries if available
        self.header_link = []
        
        if self.header.has_key("LINK"):
            for line in self.header["LINK"]:
                
                # use columwise definition
                entries = []
                # 1 - 6 Record name "LINK "
                entries.append( line[0:6].strip() )
                # 13 - 16 Atom name1 Atom name.
                entries.append( line[12:16].strip() )
                # 17 Character altLoc1 Alternate location indicator.
                entries.append( line[16] )
            # 3 # 18 - 20 Residue name resName1 Residue name.
                entries.append( line[17:20].strip() )
            # 4 # 22 Character chainID1 Chain identifier.
                entries.append( line[21] )
            # 5 # 23 - 26 Integer resSeq1 Residue sequence number.
                entries.append( line[22:26].strip() )                
                # 27 AChar iCode1 Insertion code.
                entries.append( line[26] )
                # 43 - 46 Atom name2 Atom name.
                entries.append( line[42:46].strip() )
                # 47 Character altLoc2 Alternate location indicator.
                entries.append( line[46] )
            # 9 # 48 - 50 Residue name resName2 Residue name.
                entries.append( line[47:50].strip() )
            #10 # 52 Character chainID2 Chain identifier.
                entries.append( line[51] )
            #11 # 53 - 56 Integer resSeq2 Residue sequence number.
                entries.append( line[52:56].strip() )
                # 57 AChar iCode2 Insertion code.
                entries.append( line[56] )
                # 60 - 65 SymOP sym1 Symmetry operator atom 1.
                entries.append( line[59:65].strip() )
                # 67 - 72 SymOP sym2 Symmetry operato
                entries.append( line[66:72].strip() )

                # entries for resname and reidue number are required
                if entries[3] != '' and entries[5] != '' and \
                   entries[9] != '' and entries[11] != '':
                    entries[5] = int(entries[5])
                    entries[11] = int(entries[11])
                else:
                    s =  'ERROR in \"pdb_from_biopython.__analyse_header\": ' 
                    s += 'invalid LINK entry'
                    s += ' "' + line + '"'                    
                    #print '# ERROR in \"pdb_from_biopython.__analyse_header\": ' 
                    #print '  invalid LINK entry'
                    #print '  "' + line + '"'
                    self.header_link = []
                    raise self.HeaderError(s)
                
                
#                # use regular expressions
#                #LINK        MN    MN   200                 OD1 ASP A  64                
#                #LINK         SG  CYS B 388                CU   CUB B3921
#                re_link = re.compile(r'(\w+)\s+([^\s]+)\s*([^\s]+)\s*([^\s\d]?)\s*([-\d]+)\s+([^\s]+)\s*([^\s]+)\s*([^\s\d]?)\s*([-\d]+)[\s\w]*')
#                mo = re_link.match(line)
#                if mo != None:
#                    entries = mo.groups()
#                else:
#                    s =  'ERROR in \"pdb_from_biopython.__analyse_header\": ' 
#                    s += 'invalid LINK entry'
#                    s += ' "' + line + '"'                    
#                    #print '# ERROR in \"pdb_from_biopython.__analyse_header\": ' 
#                    #print '  invalid LINK entry'
#                    #print '  "' + line + '"'
#                    self.header_link = []
#                    raise self.HeaderError(s)
#                    #warn(s)
#                    #continue
                   
                # use split
                #entries = line.split()                
#                if len(entries) != 9:
#                    # trying to fix LINK entry:
#                    # Is chain ID missing?
#
#                    #if isinstance(entries[3], str) == False:
#                    if entries[3].isdigit():
#                        entries.insert(3, "-")
#                    if entries[7].isdigit():
#                        entries.insert(7, "-")
              
                # Skip entry if chain id does not exist in pdb file. This can
                # occur if biounit files are used, where parts of the protein
                # have been deleted.
                chain_list = {' ' : 1}
                for c in self.struct.get_list():
                    chain_list[c.get_id()] = 1
                    
                if not chain_list.has_key(entries[4]) or \
                        not chain_list.has_key(entries[10]):
                    continue
                
                # Skip entry if "altLoc Alternate location indicator" is not
                # " " or A
                if entries[2] != " " and entries[2] != "A":
                    continue
                if entries[8] != " " and entries[8] != "A":
                    continue
              
                pep_found = 0
                for i, item1 in enumerate(self.header_link):
                    for j, item2 in enumerate(item1):
                        
                        if item2["name"] == entries[3] and \
                           item2["chain"] == entries[4] and \
                           item2["resid"] == entries[5]:
                               
                            new_monomer = {}
                            new_monomer["name"] = entries[9]
                            new_monomer["chain"] = entries[10]
                            new_monomer["resid"] = entries[11]
                            self.header_link[i].append(new_monomer)
                            pep_found = 1
                            break
                    
                    if pep_found == 1:
                        break
                
                if pep_found == 0:
                    new_monomer1 = {}
                    new_monomer1["name"] = entries[3]
                    new_monomer1["chain"] = entries[4]
                    new_monomer1["resid"] = entries[5]
                    new_monomer2 = {}
                    new_monomer2["name"] = entries[9]
                    new_monomer2["chain"] = entries[10]
                    new_monomer2["resid"] = entries[11]
                    self.header_link.append([new_monomer1, new_monomer2])
                    
#        for entry in self.header_link:
#            print "#"
#            for res in entry:
#                print res["name"] + "  " + res["chain"] + " " + res["resid"]
        
        
    # model:
    #    0 : if the structure is an NMR structure the first frame is used
    # name : An arbitrary name used to describe the structure.
    def parse_pdb_file(self, filename, name='PROTEIN', quiet=False,\
            model=-1):
                
        reg = re.compile(r'.+([\d\w]{4})\.pdb\d')
        reg_m = reg.match(filename)
        if reg_m != None:
            self.pdb_code_from_filename = reg_m.groups()[0].upper()
        else:
            self.pdb_code_from_filename = 'unknown'

        
        self.parser = PDBParser(QUIET=quiet)
        # to prevent the parser to print unnecessary warnings
        self.parser._handle_PDB_exception = self._handle_PDB_exception
        
        # to use the functionalliy of XStructure
        #self.all_structs = XStructure.from_file(f)
        # If files are not zipped.
        #f = open(filename)
        f = CleverFile(filename)
        self.all_structs = self.parser.get_structure(name, f)
        
        
        if model == -1:
            self.struct = self.all_structs[0]
        else:
            self.struct = self.all_structs[model]
        
        
        f.seek(0)
        
        self.header = {}
        self.header_fields = []
        for line in f:
            fields = line.split()
            if fields[0] != "ATOM" and fields[0] != "TER" and fields[0] != "HETATM":
                # bio files have no REMARK entries anymore.
                if fields[0] == 'REMARK':
                    entry_name = fields[0] + ' ' + fields[1]
                else:
                    entry_name = fields[0]

                if self.header.has_key(entry_name):
                    self.header[entry_name].append( line )
                else:
                    self.header[entry_name] = []
                    self.header[entry_name].append( line )
                self.header_fields.append(entry_name)
        
        self.__analyse_header()
        
        self.c_file_parsed = 1
        f.close()


#        self.is_nmr = False
#        if self.header.has_key('EXPDTA'):
#            if self.header['EXPDTA'][0].find('NMR') != -1:
#                #print "NMR structure found: " + filename
#                self.is_nmr = True
            
            
#        if model == -1:
#            if len(self.all_structs) > 1 and self.header['EXPDTA'][0].find('NMR') == -1:
#                self.all_structs_merged = self.merge_models()
#                
#                self.struct = self.all_structs_merged[0]
#                
##                # save structure
##                io = PDBIO()
##                io.set_structure(self.all_structs_clean)
##                io.save('/user/tmeyer/bio_struct.pdb')

        
        return 1
    
    def write_pdb_file(self, filename, merged=True):
        #PDBIO.save = pdb_from_biopython.save
        io = PDBIO()
        # replace original function with modified one.
        io.save = self.save
        
        if self.c_file_parsed == 1:
            
            if merged and self.all_structs_merged != None:
                io.set_structure(self.all_structs_merged)
                
                # Add some information from the original header.
                org_header = []
                org_header_fields = ['HEADER', 'TITLE', 'COMPND', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'JRNL', 'HETNAM',
                                     'FORMUL', 'CAVEAT']
                #org_header.append('')
                #org_header.append('The following remarks are an excerpt from the header of the original PDB file.')
                for key in org_header_fields:
                    if self.header.has_key(key):
                        for line in self.header[key]:
                            org_header.append( line.strip('\n') )


                header = []
                header.append('')
                #header.append('')
                # Add status.
                header.append("This is a modified biological assembly. All models have been merged ")
                header.append("into a single one. Chain and residue IDs have been changed if necessary.")
                header.append("The header above is just an excerpt from the header of the original")
                header.append("biological assembly.")
                

                header.append('')
                header.append('The file was created by \'pdb-bio-merger\'.')
                header.append('version: 1.10')
                
                header.append('PDB-code of the original structure: ' + self.pdb_code_from_filename)
                
                if not self.too_many_chains:
                    header.append('status: ' + 'ok')
                else:
                    header.append('status: ' + 'Too many chains, chain ID is set to \' \'')
                header.append('')
                
                if self.too_many_chains:
                    message = 'WARNING: There are too many chains in the '
                    message += 'structure to represent them with a single character.'
                    header.append(message)
                    message = '         Chain ID is set to " "! Please use segname instead.'
                    header.append(message)
                    header.append('')
                
                header += self.change_log

                header.append('')
                header.append('HET-atoms deleted (distanced < 0.5A to an other atom): '\
                        + str(len(self.deleted_atoms)) )
                        
                if len(self.deleted_atoms) > 0:
                    for del_atm in self.deleted_atoms:
                        
                        id_str = str( del_atm['model'] )       + ' ' \
                                 + del_atm['resname'].rjust(4) + ' ' \
                                 + del_atm['chain']            + ' ' \
                                 + '{0:.0f}'.format( del_atm['resid'] ).rjust(5) + ' '\
                                 + del_atm['name']
                        header.append(id_str)
                    header.append('')
                else:
                    #header.append('None')
                    header.append('')
            
            elif merged and self.all_structs_merged == None:
                s =  'ERROR in \"pdb_from_biopython.write_pdb_file\": '
                s += 'Merge flag is set, but no merged structure is available.'
                raise self.IO_Error(s)
                
            else:
                io.set_structure(self.all_structs)
                
                org_header = []
                for key in self.header_fields:
                    for line in self.header[key]:
                        org_header.append( line.strip('\n') )

                header = []                
                
            io.save(io, filename, header_list=header, org_header_list=org_header)
            
        else:
            s =  'ERROR in \"pdb_from_biopython.write_pdb_file\": '
            s += 'No structure loaded that could be written.'
            raise self.IO_Error(s)
    
    # This is a copy of the 'PDBIO.save' function. An option to store a header is added.
    def save(pdb_from_biopython, self, file, select=Select(), write_end=0, header_list=[], org_header_list=[]):
        """
        @param file: output file
        @type file: string or filehandle 

        @param select: selects which entities will be written.
        @type select: 
            select hould have the following methods:
                - accept_model(model)
                - accept_chain(chain)
                - accept_residue(residue)
                - accept_atom(atom)
            These methods should return 1 if the entity
            is to be written out, 0 otherwise.

            Typically select is a subclass of L{Select}.
        """
        get_atom_line=self._get_atom_line
        if isinstance(file, basestring):
            fp=open(file, "w")
            close_file=1
        else:
            # filehandle, I hope :-)
            fp=file
            close_file=0

        # write original header
        for line in org_header_list:
            fp.write(line + '\n')
            
        # write custom header
        for line in header_list:
            fp.write('REMARK 99  ' + line + '\n')
        
        # multiple models?
        if len(self.structure)>1 or self.use_model_flag:
            model_flag=1
        else:
            model_flag=0
        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            # necessary for ENDMDL 
            # do not write ENDMDL if no residues were written
            # for this model
            model_residues_written=0
            atom_number=1
            if model_flag:
                fp.write("MODEL      %s\n" % model.serial_num)
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id=chain.get_id()
                # necessary for TER 
                # do not write TER if no residues were written
                # for this chain
                chain_residues_written=0
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode=residue.get_id()
                    resname=residue.get_resname()  
                    segid=residue.get_segid()
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            chain_residues_written=1
                            model_residues_written=1
                            # The following lines have been modified by Tim: 
                            #s=get_atom_line(atom, hetfield, segid, atom_number, resname,
                            #    resseq, icode, chain_id[0])
                            if atom_number > 99999:
                                index = 99999
                            else:
                                index = atom_number
                            if len(chain_id) != 1:
                                c_id = ' '
                            else:
                                c_id = chain_id
                            s=get_atom_line(atom, hetfield, segid, index, resname,
                                resseq, icode, c_id)
                                
                            fp.write(s)
                            atom_number=atom_number+1
                if chain_residues_written:
                    fp.write("TER\n")
            if model_flag and model_residues_written:
                fp.write("ENDMDL\n")
            if write_end:
                fp.write('END\n')
        if close_file:
            fp.close()    
    
    
    def incr_chainid(self, cid):
        
        if len(cid) != 1:
            s =  'ERROR in \"pdb_from_biopython.incr_chainid\": '
            s += 'invalid chain ID (wrong length):'
            s += 'Chainid to be increased: ' + cid
            raise self.MergeError(s)
            
        uc = ord( cid.decode() )
        
        
        if uc == 90:
            # Z -> a
            uc = 97
        elif uc == 122:
            # z -> 0
            uc = 48
        else:
            uc += 1
        
        if not ( (65 <= uc <= 90) or (97 <= uc <= 122) or (48 <= uc <= 57)):
            s =  'ERROR in \"pdb_from_biopython.incr_chainid\": '
            s += 'invalid chain ID (unicode not in range):'
            s += 'Chainid to be increased: ' + cid
            raise self.MergeError(s)
            
        return unichr(uc).encode()

    
    def merge_models(self, overwrite=False):
        """
        Merges all models found in self.all_structs. Usefull to convert a\
        multi-model biological assembly pdb from the PDB-database into a
        single-model structure.
        
        Parameter:
            overwrite: True -> Overwrite self.all_structs_merged if it exits.
                       False-> (Default) do nothing, if self.all_structs_merged exits.
        
        -> Returns a Bio.PDB.Structure object, containing a single model with 
           the merged structure.pdb.merge
        
        
        ---------------------------------------------------------------------
        Chains are renamed obeying the following rules:
            - The chain IDs of protein chains in the first model are kept.
            - Starting with the second model the chain IDs are increased as 
              long, until an unused ID is found. (Chain IDs used for HETATOM 
              are overwritten.)
              
            - Chain IDs of pure HETATOM chains in the original pdb-file are
              ignored. Chain IDs higher than the last protein-chain ID are 
              used for them.
              
        The chain assignment is changed obeying the following rules:
            - A HEATATOM in a protein chain, that is covalent bound to the 
              protein chain remain in this chain. It can be bound indirect
              through other HEATOMs.
            - A HEATATOM in a protein chain, that is NOT bound to the protein
              chain are moved to a new chain.
            - A HETATOM that is NOT in a protein chain is placed in a new 
              chain.
            - If HEATATOMs in chains without protein are covalent bound, their
              chains are merged.
            - All monomeric HETATOM residues are merged into one chain.
            - All water molecules are placed in one chain. This chain gets the 
              highest unsused chain ID.
             
        Resids are changed obeying the following rules:
            - ATOM resids are not changed.
            - HETATOM resids are increased if the RESID exists already in the 
              same chain. This is needed when chains are merged.
            - Waters are completely renumbered, starting with 1.
              
        Residues are completely removed, if all of their atoms do overlap 
        (distance < 0.5 Angtröm). Residues with higer model IDs are deleted
        preferably. Only the first altloc atoms are compared (' ' or 'A'), if 
        there is more than one conformation available.

        If only some atoms of two residue overlap, an error is raised.
        
        """
        
        
        if self.all_structs == None:
            s =  'ERROR in \"pdb_from_biopython.merge_models\": '
            s += 'Cannot merge structures if no structures have been loaded.'
            raise self.MergeError(s)
        
        # Only x-ray structures are merged.
        is_xray = False
        if self.header.has_key('EXPDTA'):
            if self.header['EXPDTA'][0].find('X-RAY') == -1:
                s =  'ERROR in \"pdb_from_biopython.merge_models\": '
                s += 'Structure is not resolved by x-ray diffraction: '
                s += self.header['EXPDTA'][0].strip('EXPDTA').strip(' ')
                raise self.MergeError(s)
        
            
        if self.all_structs_merged != None and not overwrite:
            # Merged structure exits already and overwrite parameter is not set.
            return -1


            
        for m, model in enumerate(self.all_structs):
            for chainid in model.child_dict.keys():
                for res in model.child_dict[chainid]:
                    resid   = res.get_id()[1]
                    resname = res.get_resname()

                    id_str = str(m) + ' ' \
                             + resname.rjust(4) + ' ' \
                             + str(chainid)                              + ' ' \
                             + '{0:.0f}'.format( resid ).rjust(5)
                    res.old_identity = id_str
            

        # Count the number of atoms to make sure, that no atom is lost.
        noa_start = 0
        for i,m in enumerate(self.all_structs):
            x = [x for x in m.get_atoms()]
            noa_start += len(x)
       
       
        #print "### starting ###"
       
        ###############################################################
        # preparation step: Generate a List of bonds between HETATOM residues 
        #                   and other HETATOM or ATOM residues.
        ##############################################################
       
        all_het_connections = []
        for struct in self.all_structs:
            model_connections = {}
            
            # create residue list
            res_list = []
            for chainid in struct.child_dict.keys():
                for i, res in enumerate( struct.child_dict[chainid].child_list ):
                    res_list.append( (chainid, res) )
              
            for i, (chainid, res) in enumerate( res_list[:-1] ):
                # Compare "res" with next residue in the chain.
                (chainid_comp, res_comp) = res_list[i+1]
                
                # Check if the two residues are connected.
                is_bond = False
                # Do not compare ATOM with ATOM residues.
                if res.get_id()[0] != ' ' or res_comp.get_id()[0] != ' ':
                    for atom in res.child_list:
                        if atom.element == 'H':
                            continue                        
                        
                        coord = atom.get_coord()
                        
                        for atom_comp in res_comp.child_list:
                            if atom.element == 'H':
                                continue
                            
                            coord_comp = atom_comp.get_coord()
                            
                            dist = coord_comp - coord
                            
                            if 0.5 < np.sqrt(np.dot( dist, dist )) < 1.7:
                                # Bond between the atoms found.
                                is_bond = True
                                break
                            
                        if is_bond:
                            break
                        
                
                # If they are bond add an connection entry.
                # A residue is define by: (chainid and residue ID)
                res_description      = ( chainid, res.get_id() )
                if not model_connections.has_key( res_description ):
                    model_connections[ res_description ] = {}
                    
                if is_bond:
                    res_comp_description = ( chainid_comp, res_comp.get_id() )
                    model_connections[ res_description ][ res_comp_description ] = True
                    
                    if not model_connections.has_key( res_comp_description ):
                        model_connections[ res_comp_description ] = {}
                    model_connections[ res_comp_description ][ res_description ] = True
            
            
            # To be complete add the last residue to the connection dictionary.
            (chainid, res) = res_list[-1]
            res_description      = ( chainid, res.get_id() )
            model_connections[ res_description ] = {}
            
            all_het_connections.append( model_connections )
       
        #print '### preparation step done ###'
       
       
       
        #######################################################################
        # first step: Split ATOM and HETATOM chains.
        #             Store alle chains seperatly, they will be merged later.
        #             Only chain IDs of ATOMs are kept.
        #######################################################################
        
        # One entry for every model
        all_p_chains = []
        all_h_chains = []
        
        # each model
        for struct_id, struct in enumerate( self.all_structs ):
            p_chains = []
            h_chains = []

            chain_list = [c[0] for c in struct.child_dict.keys()]
            chain_list.sort()
            
            # each chain
            for cid in chain_list:
                # Each chain is split into ATOM and HETATOM residues
                
                # This will be used for ATOMS (chain id is kept for now)
                atom_residues = Chain(cid)
                # This will be used for HEATOMS (chain id is lost)
                het_residues  = Chain(cid)
                
                # each residue
                for res in struct.child_dict[cid].child_list:
                    
                    if res.get_id()[0] == ' ':
                    #if is_std_aa(res.get_resname()):
                        # ATOM residue found
                        atom_residues.add(res)
                    else:
                        # HETATOM residue found
                        
                        res_description = ( cid, res.get_id() )
                        for (x, bound_res_id) in all_het_connections[ struct_id ][ res_description ]:
                            if bound_res_id[0] == ' ':
                                # HETATOM residue is bound to ATOM residue
                                atom_residues.add(res)
                                break
                        else:
                            # Will be deleted from this list later, if neccesary.
                            # This is required to keep the residue order.
                            atom_residues.add(res)
                            # Residues in here will be checked later.
                            het_residues.add(res)
                         


                # Check if the HETATOMs are covalent bound to another ATOM
                # residue in the same chain. -> If yes: Keep it in the 
                # protein chain.
                something_changed = True
                while something_changed:
                    something_changed = False
                    
                    for i, res in enumerate(het_residues.child_list):
                        res_description = ( cid, res.get_id() )
                        
                        for res_comp in atom_residues.child_list:
                            res_comp_description = ( cid, res_comp.get_id() )

                            # 'res_comp' is only a valid connection partner, if 
                            # it is not contained in 'het_residues'. That means
                            # it is ether a ATOM Residue or a resdidue that is
                            # directly or indirectly connected to an ATOM residue.
                            if not het_residues.has_id( res_comp.get_id() ):
                            
                                if all_het_connections[ struct_id ]\
                                        [ res_description ].has_key( res_comp_description ):
                                            
                                    het_residues.detach_child( res.get_id() )
                                    something_changed = True
                                    break
                                
                        if something_changed:
                            break

                # remove all entries from atom_residues that are still in het_residues
                something_changed = True
                while something_changed:
                    something_changed = False
                    
                    for res_atom in atom_residues.child_list:
                        if het_residues.has_id( res_atom.get_id() ):
                            atom_residues.detach_child( res_atom.get_id() )
                            something_changed = True
                            break
                       
                if len(atom_residues.child_list) > 0:
                    p_chains.append(atom_residues)
                if len(het_residues.child_list) > 0:
                    h_chains.append(het_residues)

                
            all_p_chains.append(p_chains)
            all_h_chains.append(h_chains)

        #print "### first step done ###"


        

        ######################################################################
        # second step: Look for covalent bound HETATOM residues for each model.
        #              Make sure, they are in the same chain.
        ######################################################################
        
        # merge chains that are covalent connected through at least one residue pair
        all_h_chains_new = []
        
        # All monomeric HETATOMs are merged into one chain.
        monomer_chain_list = []
        
        # All waters are merged into one chain
        water_chain = Chain('')
        used_water_residues = {}        
        
        # each model
        for m_id, model in enumerate( all_h_chains ):
            
            # generate list with all residues in the model
            res_descr_list = []
            for chain in model:
                chainid = chain.get_id()
                
                for res in chain.child_list:
                    res_description = ( chainid, res.get_id() )
                    res_descr_list.append( res_description )
            
            res_done = {}
            
            # first get all waters
            for res_description in res_descr_list:
                if res_description[1][0] == 'W':
                    res_done[res_description] = True
                    
                    #resid = res_description[1][1]
                    resid = 1
                    while used_water_residues.has_key(resid):
                        resid += 1
                    used_water_residues[resid] = True
                    
                    res = self.all_structs[m_id]\
                            .child_dict[ res_description[0] ]\
                            .child_dict[ res_description[1] ]
                    
                    res.id = (res.id[0],\
                              resid,\
                              res.id[2])
                              
                    # The attribute "chainid" is not standart for a
                    # chiain-object. But it will be neded later.
                    res.chainid_old = res_description[0]
                    res.resid_old   = res_description[1][1]
                    
                    water_chain.add( res )
                            
                        
                        
            while True:
                new_chain = Chain('')
                used_resid = {}
            
                # each residue
                # Find a new residue to start with
                for res_description in res_descr_list:
                    if not res_done.has_key( res_description ):
                        break
                else:
                    # All residues have been added to the new chain list. -> stopping
                    break
                        
                to_look_at_list = [res_description]
                res_done[res_description] = True
                
                
                # Make sure, that no resid is used twice.
                resid = res_description[1][1]
                while used_resid.has_key(resid):
                    resid += 1
                used_resid[resid] = True
                
                res = self.all_structs[m_id]\
                        .child_dict[ res_description[0] ]\
                        .child_dict[ res_description[1] ]
                
                res.id = (res.id[0],\
                          resid,\
                          res.id[2])
                
                # The attribute "chainid" is not standart for a
                # chiain-object. But it will be neded later.
                res.chainid_old = res_description[0]
                res.resid_old   = res_description[1][1]
                
                new_chain.add( res )
                
                while len(to_look_at_list) > 0:
                    
                    res_description = to_look_at_list.pop()
                    
                    for res_desc_partner in all_het_connections[m_id][ res_description ].keys():
                        if not res_done.has_key( res_desc_partner ):
                            
                            to_look_at_list.append( res_desc_partner )
                            res_done[res_desc_partner] = True
                            
                            res = self.all_structs[m_id]\
                                    .child_dict[ res_desc_partner[0] ]\
                                    .child_dict[ res_desc_partner[1] ]
                            
                            while used_resid.has_key(resid):
                                resid += 1
                            used_resid[resid] = True
                            
                            res.id = (res.id[0],\
                                      resid,\
                                      res.id[2])
                            
                            # The attribute "chainid" is not standart for a
                            # chiain-object. But it will be needed later.
                            res.chainid_old = res_desc_partner[0]
                            res.resid_old   = res_desc_partner[1][1]
                            
                            new_chain.add( res )

                if len(new_chain.child_list) == 1:
                    monomer_chain_list.append( new_chain.child_list[0] )
                else:
                    all_h_chains_new.append(new_chain)
        
        # monomeric HETATOMs
        monomer_chain = Chain('')
        used_monomer_residues = {}

        for res in monomer_chain_list:
            resid = res.get_id()[1]
            while used_monomer_residues.has_key(resid):
                resid += 1
                    
            res.id = (res.id[0],\
                      resid,\
                      res.id[2])
                     
            monomer_chain.add(res)
            used_monomer_residues[resid] = True
                      
        # The chain with monomeric HETATMS should appear befor the water 
        # chainin the pdb file.
        all_h_chains_new.append(monomer_chain)
        
        # Waters should appear ath the end of the pdb file.
        all_h_chains_new.append(water_chain)



        #######################################
        # third step: Delete overlapping atoms. 
        #######################################

        ### ATOM collisions ###
        # Todo
        
        ### HETATOM collisions ###
        
        # Stores the deleted atoms.
        deleted_atoms = []

        noa = 0
        atom_list = []
        for c_nr, chain in enumerate(all_h_chains_new):
            for residue in chain.child_list:
                res_id = residue.get_id()
                for atom in residue.child_list:
                    # Add the non-standart attribute "chainid" to the atom.
                    # It is taken from the residue attribute "chainid", that
                    # that was set previously.
                    atom.chainid_old = residue.chainid_old
                    atom.resid_old   = residue.resid_old
                    atom.resname     = residue.get_resname()
                    atom_list.append( (c_nr, res_id ,atom) )
                    noa += 1
                    
        atom_delete_dict = {}
        
      

        
#        #for i, (c_nr, res_id, atom) in enumerate(atom_list[:-1]):
#        for i  in range( len(atom_list) - 1 ):
#            (c_nr, res_id, atom) = atom_list[i]
#            
#            coord = atom.get_coord()
#            
#            print str(i) + "/" + str(len(atom_list))
#            print coord
#            
#            #for (c_comp_nr, res_comp_id, atom_comp) in atom_list[i+1:]:
#            for j  in range( i+1, len(atom_list) ):
#                #(c_comp_nr, res_comp_id, atom_comp) = atom_list[j]
#                                
#                #coord_comp = atom_comp.get_coord()
#                
#                #dist = coord_comp - coord
#                
#                dist = atom_list[j][2].get_coord() - coord
#                
#                if abs(dist[0]) > 0.6:
#                    continue
#                if abs(dist[1]) > 0.6:
#                    continue
#                if abs(dist[2]) > 0.6:
#                    continue                
#                
#                if np.sqrt(np.dot( dist, dist )) < 0.5:
#                    (c_comp_nr, res_comp_id, atom_comp) = atom_list[j]
#                    
#                    if atom_delete_dict.has_key( (c_comp_nr, res_comp_id, atom_comp.get_id()) ):
#                        continue
#
#                    print "atom to delete: " + str( (c_comp_nr, res_comp_id, atom_comp.get_id()) )
#                    
#                    # Bond between the atoms found.
#                    del_atom = {}
#                    del_atom['name']  = atom_comp.get_id()
#                    del_atom['chain'] = atom_comp.chainid_old
#                    del_atom['resid'] = atom_comp.resid_old
#                    del_atom['resname'] = atom_comp.resname
#                    del_atom['model'] = c_comp_nr
#                    
#                    deleted_atoms.append(del_atom)
#                    atom_delete_dict[ (c_comp_nr, res_comp_id, atom_comp.get_id()) ] = True
#                        

        # Find overlapping atoms.
        atom_coord_list = [x[2].get_coord() for x in atom_list]
        atom_coord_list = np.array(atom_coord_list, dtype=np.double)
        
        #del_atoms = all_vs_all_atoms.find_collisions(atom_coord_list, 0.5)
        if len(atom_coord_list) > 0:
            del_atoms = find_collisions(atom_coord_list, 0.5)
        else:
            del_atoms = []
        
        # Store the found atoms, that will not occur in the merged structure.
        for i in del_atoms:
            (c_comp_nr, res_comp_id, atom_comp) = atom_list[i]

            #print "atom to delete: " + str( (c_comp_nr, res_comp_id, atom_comp.get_id()) )
                    
            del_atom = {}
            del_atom['name']  = atom_comp.get_id()
            del_atom['chain'] = atom_comp.chainid_old
            del_atom['resid'] = atom_comp.resid_old
            del_atom['resname'] = atom_comp.resname
            del_atom['model'] = c_comp_nr
            
            if not atom_delete_dict.has_key( (c_comp_nr, res_comp_id, atom_comp.get_id()) ):
                deleted_atoms.append(del_atom)
                atom_delete_dict[ (c_comp_nr, res_comp_id, atom_comp.get_id()) ] = True


       
        #print "### third step done ###"
        
        
        
        
        ####################################################################
        # fourth step: Merge chains again.
        #              First all ATOM chains, then the HETATOM chains
        #              Chain IDs of the protein are kept as far as possible
        ####################################################################
        
        change_log = []
        change_log.append("")
        change_log.append("Single chain and residue IDs have been change as described below.")
        change_log.append("1) Model number in biological assembly.")
        change_log.append("2) Residue description in biological assembly.")
        change_log.append("3) Residue description in this file.")
        change_log.append("1  -----2-----      -----3-----")
        change_log.append("")
        
        structure_builder = StructureBuilder()
        self.structure = structure_builder
        structure_builder.init_structure('protein')
        structure_builder.init_model(0)
        
        c_id = ''
        index = 1
        
        # Keeps track of the original chain id of a chain.
        prot_chain_changes = {}
        
        
        num_chains = 0
#        for model in all_p_chains:
#            num_chains += len(all_p_chains)
#        num_chains += len(all_h_chains)
        for model in all_p_chains:
            num_chains += len(model)
        num_chains += len(all_h_chains_new)
        
        # A-Z + a-z + 0-9 minus one for the water chain.
        if num_chains > 26 * 2 + 10 - 1:
            self.too_many_chains = True
        else:
            self.too_many_chains = False
        
        
        # add protein atoms
        for m_id, model in enumerate(all_p_chains):
            for i, chain in enumerate(model):
                c_id = chain.get_id()
                
                if not self.too_many_chains:
                    while structure_builder.structure[0].has_id( c_id ):
                        c_id = self.incr_chainid( c_id )
                else:
                    # If there are too many chains, only segid is used
                    # The chain ID set here will not be saved.
                    c_id = chain.get_id() + '_' + str(m_id)

                segid = c_id
                
                structure_builder.init_chain( c_id )
                #segid = chain.get_id() + '_' + str(m_id)
                structure_builder.init_seg(segid)


                
                # transfer residues
                for res in chain.child_list:
                    res_id = res.get_id()

                    #hetflag = res_id[0]
                    field   = res_id[0]
                    resid   = res_id[1]
                    icode   = res_id[2]
                    resname = res.get_resname()
                    
                        
                    structure_builder.init_residue(resname, field, resid, icode)

                    # Log the change in chainid and resid.
                    id_str = resname.rjust(4) + ' ' \
                             + str(c_id)      + ' ' \
                             + '{0:.0f}'.format( resid ).rjust(5)
                    change_log.append( res.old_identity + '     ' + id_str )
                    


                    reg = re.compile(r'\d+\s+\w+\s+(\w)\s+[-\d]+')
                    reg_m = reg.match(res.old_identity)
                    if reg_m == None:
                        s =  'ERROR in \"pdb_from_biopython.merge_models\": '
                        s += 'Residue identity could not be parsed:\n'
                        s += res.old_identity
                        raise self.MergeError(s)

                    if not prot_chain_changes.has_key(c_id):
                        if not self.too_many_chains:
                            prot_chain_changes[c_id] = reg_m.groups()[0]
                        else:
                            prot_chain_changes[segid] = reg_m.groups()[0]
                    else:
                        inconsitence = False
                        if not self.too_many_chains:
                            if prot_chain_changes[c_id] != reg_m.groups()[0]:
                                inconsitence = True
                        else:
                            if prot_chain_changes[segid] != reg_m.groups()[0]:
                                inconsitence = True
                                
                        if inconsitence:
                            s =  'ERROR in \"pdb_from_biopython.merge_models\": '
                            s += 'Renaining scheme for protein chains is inconsistent.'
                            raise self.MergeError(s)



                    # transfer atoms
                    for atm in res.child_list:
                        atm_id = atm.get_id()
                        
                        structure_builder.init_atom(\
                                atm.get_name(),\
                                atm.get_coord(),\
                                atm.get_bfactor(),\
                                atm.get_occupancy(),\
                                atm.get_altloc(),\
                                atm.get_fullname(),\
                                serial_number=index,\
                                element=atm.element)
                                
                        index += 1
                
        counter = 0
        
        
        # add HETAOMS
        het_chains_descr = {}
        
        num_of_chains = len(all_h_chains_new)        
        
        #for m_id, model in enumerate(all_h_chains):
        for c_nr, chain in enumerate(all_h_chains_new):
            
            if len(chain.child_list) == 0:
                # Chain is emtpy -> no water / monomers
                continue
            
            if not self.too_many_chains:
                if c_id == '':
                    # No Protein chains -> set first HETATM chain to A.
                    c_id = 'A'
                else:
                    while structure_builder.structure[0].has_id( c_id ):
                        c_id = self.incr_chainid( c_id )
                    
                segid = c_id + 'HET'
            else:
                # If there are too many chains, only segid is used.
                # The chain ID set here will not be saved.
                c_id = 'H' + str(counter)
                segid = 'HET' + str(counter)
            
            structure_builder.init_chain(c_id)
            structure_builder.init_seg(segid)
            
            if chain.child_list[0].get_id()[0] == 'W':
                if not self.too_many_chains:
                    het_chains_descr[c_id] = 'water'
                else:
                    het_chains_descr[segid] = 'water'
            elif counter == num_of_chains - 2:
                if not self.too_many_chains:
                    het_chains_descr[c_id] = 'HETATM monomers'
                else:
                    het_chains_descr[segid] = 'HETATM monomers'
            else:
                if not self.too_many_chains:
                    het_chains_descr[c_id] = 'HETATM chain'
                else:
                    het_chains_descr[segid] = 'HETATM chain'

            counter += 1
                    
            used_resid = {}
            
            # transfer residues
            too_many_resid = False
            for res in chain.child_list:
                res_id = res.get_id()
                
                
                #hetflag = res_id[0]
                field   = res_id[0]
                
                if not too_many_resid:
                    resid = res_id[1]
                else:
                    resid = 0
                    
                icode   = res_id[2]
                resname = res.get_resname()
                
                
                number_of_atoms_is_residue = len(res.child_list)
                atoms_left_in_residue      = len(res.child_list)
                
                # Make sure, that no resid is used twice in one chain.
                while used_resid.has_key(resid):
                    resid += 1
                used_resid[resid] = True
                
                
                # Not more than 4 digits are allowed in the standard pdb format.
                if resid > 9999:
                    # start a new chain
                    if not self.too_many_chains:
                        while structure_builder.structure[0].has_id( c_id ):
                            c_id = self.incr_chainid( c_id )
                    else:
                        # If there are too many chains, only segid is used
                        c_id = ' ' + str(counter)
                    
                    structure_builder.init_chain( c_id )
                    segid = 'HET' + str(counter)
                    structure_builder.init_seg(segid)
        
                    counter += 1
                    
                    used_resid = {}
                    resid = 0
                    used_resid[resid] = True
                    
                    too_many_resid = True
                    

                
                
                
                structure_builder.init_residue(resname, field, resid, icode)

                # Log the change in chainid and resid.
                id_str = resname.rjust(4) + ' ' \
                         + str(c_id)      + ' ' \
                         + '{0:.0f}'.format( resid ).rjust(5)
                change_log.append( res.old_identity + '     ' + id_str )

                
                # transfer atoms (ordered and disordered)
                for atom in res.child_list:
                    
                    if atom.is_disordered() == 0:
                        altloc_atoms = [atom]
                    else:
                        altloc_atoms = atom.disordered_get_list()
                    
                    for atm in altloc_atoms:
                        atm_id = atm.get_id()
                        
                        if atom_delete_dict.has_key( (c_nr, res_id, atm_id) ):
                            atoms_left_in_residue -= 1
                            continue
                        
                        structure_builder.init_atom(\
                                atm.get_name(),\
                                atm.get_coord(),\
                                atm.get_bfactor(),\
                                atm.get_occupancy(),\
                                atm.get_altloc(),\
                                atm.get_fullname(),\
                                serial_number=index,\
                                element=atm.element)
                                
                        index += 1
                
                # Check that only complete residues are deleted.
                if atoms_left_in_residue != number_of_atoms_is_residue and \
                       atoms_left_in_residue != 0:
                    # Add exceptions for common uncritical residues.
                    if resname not in ['SO4', 'HOH', 'PO4', 'GOL', 'MPD', 'TRS', 'EDO']:
                        s =  'ERROR in \"pdb_from_biopython.merge_models\": '
                        s += 'Residue ' + resname + '_' + res.chainid_old + '-' + str(res.resid_old)
                        s += ' is about to be deleted partially.'
                        raise self.MergeError(s)
              
        
        change_log.append('')
        
        
        #print "### fourth step done ###"
        
        
        # Add chain translation into header.
        chain_table = []
        chain_table.append('')
        chain_table.append('Chain IDs have been changed as follows. The origin of HETATM chains')
        chain_table.append('can be mixed and is therefore set to \'-\'.')
        chain_table.append('new chain ID <-- original chain ID')
        if self.too_many_chains:
            chain_table.append('NO CHAIN IDS HAVE BEEN ASSIGNED!')
        else:
            for chain in sorted(prot_chain_changes.keys()):
                #if not self.too_many_chains:
                chain_table.append(chain + ' <-- ' + prot_chain_changes[chain])
                #else:
                #    chain_table.append('-' + ' <-- ' + prot_chain_changes[chain])
            for chain in sorted(het_chains_descr.keys()):
                #if not self.too_many_chains:
                chain_table.append(chain + ' <-- - (' + het_chains_descr[chain] + ')')
                #else:
                #    chain_table.append('-' + ' <-- - (' + het_chains_descr[chain] + ')')
        chain_table.append('')
        

        pos = 0
        for entry in chain_table:
            change_log.insert(pos, entry)
            pos += 1
        
        
        
        #####################################
        # Final check that no atoms are lost.
        #####################################
        noa_end = 0
        x = [x for x in structure_builder.structure.get_atoms()]
        noa_end = len(x)
        
        if noa_end + len(deleted_atoms) != noa_start:
            s =  'ERROR in \"pdb_from_biopython.merge_models\": ' 
            s += 'Some atoms are lost!'
            s += '\n Atoms in the start structure: ' + str(noa_start)
            s += '\n Deleted atoms: '                + str(len(deleted_atoms))
            s += '\n Atoms in the final structure: ' + str(noa_end)
            raise self.MergeError(s)
        
        
#        if len(deleted_atoms) > 0:
#            print "Some atoms have been deleted: (" + str(len(deleted_atoms)) + ")"
#        for atom in deleted_atoms:
#            print atom
        
        self.all_structs_merged = structure_builder.structure
        
        self.struct = structure_builder.structure[0]
        
        self.change_log    = change_log
        self.deleted_atoms = deleted_atoms
        
        return structure_builder.structure

        
        
   
    def count_resname_occurence(self, name, check_for_errors=0):
        
        if len(name.split()) == 1:
            
            # ligand is a monomer
            
            occurence = []
            
            for resid in self.struct.get_residues():
                resname = resid.get_resname()                

                if resname.strip() == name.strip():                    
                    
                    # if ligand is a standart amino acid, it has to be a HETATOM
                    # or the only residue in its chain
                    if is_std_aa(name):
                        if resid.id[0] == ' ' and len(resid.parent) > 1:
                            # residue seems to be part of the protein chain
                            continue
                            
                    
                    new_occurence = {}
                    new_occurence["name"] = name
                    new_occurence["chain"] = str( resid.get_full_id()[2] )
                    new_occurence["resid"] = \
                            int( resid.get_full_id()[3][1] )
                    new_occurence["atom"] = resid

                    if check_for_errors == 1:
                        # Check number of atoms in neighbourhood, if the numbers
                        # differs to much between different occurancies of the
                        # ligand, there is probably an error in the pdb.
                        #
                        # altloc = 'B' is ignored
                        
                        # iterate over all atoms (children) of residue 'resid'
                        ref_coord = np.array([0,0,0])
                        num = 0
                        for atom in resid.get_list():
                            if atom.get_altloc() != 'B':
                                ref_coord += atom.get_coord()
                                num += 1
                        ref_coord /= num

                        new_occurence["atoms_in_sphere"] = 0                     
                        
                        for atom in self.struct.get_atoms():
                            if atom.get_altloc() != 'B':
                                coord = atom.get_coord()
                                dist = coord - ref_coord
                                if np.sqrt(np.dot( dist, dist )) < 10:
                                    new_occurence["atoms_in_sphere"] += 1
                            
                    # entry is a list to be compatible with peptide notation
                    occurence.append([new_occurence])


            # Check if two Monomers are covalent bound.
            # -> If they are, merge them.
            restart = True
            while restart:
                restart = False
            
                for i, occ in enumerate(occurence):
                    for mon in occ:
                        
                        # look in the occurances that follow.
                        for j, occ_comp in enumerate( occurence[i+1:] ):
                            j += 1 # the for loop is always start with 0 for j
                            for mon_comp in occ_comp:
                                
                                # compare the two monomers
                                for atom in mon['atom'].get_list():
                                    coord = atom.get_coord()
    
                                    for atom_comp in mon_comp['atom'].get_list():
                                        coord_comp = atom_comp.get_coord()
                                        
                                        dist = coord_comp - coord
                                        #print dist
                                        #print np.sqrt(np.dot( dist, dist ))
                                        
                                        if np.sqrt(np.dot( dist, dist )) < 1.5:
                                            # Bond between the monomers found!
                                            # Merge the whole occurence
                                            occurence[i] += occurence[j]
                                            occurence.pop(j)

                                            restart = True
                                            break
                                            
                                    if restart:
                                        break
                                if restart:
                                    break
                            if restart:
                                break
                        if restart:
                            break
                    if restart:
                        break



                                    
                                    
                        
#                            for atom in self.struct.get_atoms():
#                                if atom.get_altloc() != 'B':
#                                    coord = atom.get_coord()
#                                    dist = coord - ref_coord
#                                    if np.sqrt(np.dot( dist, dist )) < 10:
#                                        new_occurence["atoms_in_sphere"] += 1
                        
                    
                
                

            
            if check_for_errors == 1:
                # find highest number of "atoms_in_sphere" attribute
                max_ais = -1;
                for num in occurence:
                    if num[0]["atoms_in_sphere"] > max_ais:
                        max_ais = num[0]["atoms_in_sphere"]
    
                # If "atoms_in_sphere" is smaller then 0.25*max_ais, the entry is
                # rejected.
                for i,item in enumerate(occurence):
                    if item[0]["atoms_in_sphere"] < 0.25*max_ais:
                        # only works for monomer!
                        lig = occurence.pop(i)[0]
                        s = "### WARNING: Ligand entry rejected due to small " \
                             +"\"atoms_in_sphere\" value (ref.value is " \
                             + str(max_ais) + ").\n"                    
                        for key, value in lig.items():
                            s += '    ' + key + ' : ' + str(value) + '\n'
                        s += '\n'                    
                        self.warnings.warn(s)
            
            for o in occurence:
                for res in o:
                    res.pop('atom')
            
            return occurence
        else:
            ###########################
            ### ligand is a peptide ###
            ###########################
            
            pep_sequence = name.split()
                        
            # maximum number of residues that are allowed to be missing in 
            # pdb-file
            max_missing = 2
            # minimum number of a peptide monomers that must exist in the pdb
            min_length = 2
            
            # If the peptide consist only of standart amino acids, at least 4
            # residues are required.
#            is_std = 1
#            for res in pep_sequence:
#                if is_std_aa(res) == 0:
#                    is_std = 0
#                    break

            
            occurence = []
            
#            ########################################
#            ### FIRST STEP: look at link entries ###
#            ########################################
#            for link_entry in self.header_link:
#                
#                skip_link_entry = False
#                
#                for res in link_entry:
#                    
#                    if res['chain'] != ' ':
#                        ### check that entry is consistent with ATOM or HETATM entries                
#                        name = res['name']
#                        chain = res['chain']
#                        resid = res['resid']
#                        # get chain (will exist, since checked in "analyse_header")
#                        chain_dict = self.struct.child_dict[ chain ].child_dict
#                        is_inconsistent = False
#                        
#                        fe = chain_dict.keys()[0][0]
#                        le = chain_dict.keys()[0][2]
#                        k  = (fe, resid, le)
#                        if chain_dict.has_key( k ):
#                            if not name != chain_dict[ k ].get_resname():
#                                is_inconsistent = True
#                        else:
#                            is_inconsistent = True
#                        
#                        if is_inconsistent:
#                            s = "### WARNING: Link entry of ligand monomer "\
#                                + res['name'] + ":" + res['chain'] \
#                                + ":" + str(res['resid'])\
#                                + " is inconsistent with ATOM/HETATM entries! "\
#                                + "Attemting to find the correct ID.\n"
#                            self.warnings.warn(s)
#                            res['chain'] = ' '
#                            #skip_link_entry = True
#                    
#                    
#                    ### COMPLEMENT CHAIN IDENTIFIER
#                    # check that all chain ids are defined. They can be missing, if the
#                    # ligand entry was obtained from LINK entries in the pdb-header.
#                    # If the correct chain id can not be determined unambiguously, a 
#                    # warning will be send.
#                    # Todo: move this to "analyse_header"?
#                    if res['chain'] == ' ':
#                        hits = 0
#                        # try to find the corresponding residue
#                        for r in self.struct.get_residues():
#                            resname = r.get_resname()
#                            chain = str( r.get_full_id()[2] )
#                            resid = int( r.get_full_id()[3][1] )
#                                
#                            if resname == res['name'] and resid == res['resid']:
#                                hits += 1
#                                res['chain'] = chain
#                            
#                        if hits > 1:
#                            s = "### WARNING: Chain ID of ligand monomer "\
#                            + res['name'] + ":" + str(res['resid'])\
#                            + " could not be determined unambiguously! "\
#                            + "The LINK entry will be ignored.\n"
#                            self.warnings.warn(s)
#                            skip_link_entry = True
#                        if hits == 0:
#                            s = "### WARNING: Chain ID of ligand monomer "\
#                            + res['name'] + ":" + str(res['resid'])\
#                            + " could not be determined! "\
#                            + "The LINK entry will be ignored.\n"
#                            self.warnings.warn(s)
#                            skip_link_entry = True
#
#                   
#                if skip_link_entry:
#                    # skip this sequence from LINK entries
#                    continue
#                
#                
#                hetatm_sequence = []
#                for resid in link_entry:
#                    hetatm_sequence.append(resid['name'])
#
#                
#                # try with original and reversed HETATOM entry
#                reverse = 0
#                for reverse in range(0,2):
#                    if reverse == 1:
#                        hetatm_sequence.reverse()
#                
#                    # To prevent that the peptide sequence matches a part of a
#                    # hetatom sequene.
#                    if len(hetatm_sequence) > len(pep_sequence):
#                        continue
#                    
#                    # compare sequence from link entry with peptide sequence
#                    missing_left = max_missing
#                    diverging  = 0
#                    length_left = len(pep_sequence)
#                    j = 0
#                    for i, item in enumerate(pep_sequence):
#                        if j >= len(hetatm_sequence) or \
#                                hetatm_sequence[j] != pep_sequence[i]:
#                                    
#                            missing_left -= 1
#                            length_left -= 1
#                            j -= 1
#                            if missing_left < 0 or length_left < min_length:
#                                diverging = 1
#                                break
#                            
#                        j += 1
#                        
#                    if diverging == 0:
#                        occurence.append(link_entry)
                    
            
            
            ####################################        
            ### SECOND STEP: look at protein ###            
            ####################################

            # iterate over chains
            for chain in self.struct.child_list:
                chain_sequence = []
                
                # iterate over residues
                for resid in chain.get_list():
                    chain_sequence.append( resid.get_resname().strip() )
                    
                # compare sequence with ligand sequence
                # To prevent that the peptide sequence matches a part of a
                # chain sequene.
                if len(chain_sequence) > len(pep_sequence):
                    continue
                # no chance for hit
                if len(chain_sequence) + max_missing < len(pep_sequence):
                    continue
                
                missing_left = max_missing
                diverging  = 0
                length_left = len(pep_sequence)
                j = 0
                for i, item in enumerate(pep_sequence):
                    if j >= len(chain_sequence) or \
                            chain_sequence[j] != pep_sequence[i]:
#                        print chain_sequence[j]
                        #print pep_sequence[i]
                                
                        missing_left -= 1
                        length_left -= 1
                        j -= 1
                        if missing_left < 0 or length_left < min_length:
                            diverging = 1
                            break
                        
                    j += 1
                    
                if diverging == 0:
                    new_occurence = []
                    for resid in chain.get_list():
                        new_entry = {}
                        new_entry["name"] = resid.get_resname()
                        new_entry["chain"] = str( resid.get_full_id()[2] )
                        new_entry["resid"] = \
                                int( resid.get_full_id()[3][1] )
                                
                        new_occurence.append(new_entry)
                        
                    occurence.append(new_occurence)
            
            
            
            ############################################################
            ### THIRD STEP: Look at protein and allow ligand to have ###
            ###             several chains.                          ###
            ############################################################
            
            # This Method becomes a mess if the peptide consists of only one
            # type of monomer (e.g. 3GXR and "NAG NAG" / "NAG NAG NAG")
            is_monotonous = 0
            for m in pep_sequence:
                if m != pep_sequence[0]:
                    break
            else:
                # skip third step
                is_monotonous = 1
            
            if not is_monotonous:
            
                max_len = len(pep_sequence)
                min_len = len(pep_sequence) - max_missing
                if min_len < min_length:
                    min_len = min_length
                
                # get chains
                chains = self.struct.get_list()
                
                # iterate over possible start chains
                for c_start, chain in enumerate(chains):                
                    # peptide sequence is stored as residue sequence and as list of chain objects
                    chain_sequence = []                
                    chain_list = []
                    # iterate over residues in start chain
                    for resid in chain.get_list():
                        chain_sequence.append( resid.get_resname().strip() )
                        
                    chain_list.append(chain)
                    
                    # make sure that the first chain is not too long already
                    if len(chain_sequence) > max_len:
                        continue
                    
                    # concatenate chains unless the length is >max_len
                    for c_add in range(c_start+1, len(chains)):
                        length_added = len( chain_sequence ) + len( chains[c_add] )
                        
                        if length_added <= max_len:
                            # add chain to sequence and chain list
                            for resid in chains[c_add].get_list():
                                chain_sequence.append( resid.get_resname().strip() )
                            chain_list.append(chains[c_add])
                        else:
                            break
                                
                    # continue to add chains, if found sequence is too short
                    if len(chain_sequence) < min_len:
                        continue
                    
    #                print chain_sequence
    #                print pep_sequence
    
    
                    # The length of the sequence could match, so check if the sequence
                    # is the same.                    
                    missing_left = max_missing
                    diverging  = 0
                    length_left = len(pep_sequence)
                    j = 0
                    last_matching_residue = -1
                    
                    for i, item in enumerate(pep_sequence):
                        if (j >= len(chain_sequence)) or \
                                (chain_sequence[j] != pep_sequence[i]):
                            missing_left -= 1
                            length_left -= 1
                            j -= 1
                            if missing_left < 0 or length_left < min_length:
                                diverging = 1
                                break
                        else:
                            last_matching_residue = j
                            
                        j += 1
                        
                    if last_matching_residue + 1 < len(chain_sequence):
                        # Residue list has to be rebuild. Only chains at the end 
                        # are removed. If the chain at the beginning is too much
                        # it will be found in the next step of the main loop.
                        tmp_list = []
                        num = 0
                        
                        for chain in chain_list:
                            #if num + len(chain.get_list()) <= last_matching_residue + 1:
                            tmp_list.append(chain)
                            num += len(chain.get_list())
                            if num >= last_matching_residue + 1:
                                break
                        
                        chain_list = tmp_list
                        
                    if diverging == 0:
                        new_occurence = []
                        for chain in chain_list:
                            for resid in chain.get_list():
                                new_entry = {}
                                new_entry["name"] = resid.get_resname()
                                new_entry["chain"] = str( resid.get_full_id()[2] )
                                new_entry["resid"] = \
                                        int( resid.get_full_id()[3][1] )
                                        
                                new_occurence.append(new_entry)
                        
                        occurence.append(new_occurence)      
            
            #END: if not is_monotonous:



            print "---"
            for o in occurence:
                print " next:"
                for res in o:
                    print res['name'] + ' : ' \
                            + res['chain'] + ' : ' + str(res['resid'])
            print "---"
            
            
            ############################
            ### CHECK FOR UNIQUENESS ###
            ############################
            # Check that no occurence is used more than once.
            # If no chain is defined, but name and resid matches the entries are
            # considered to be the same and the entry with chain id is kept.            
            restart = 1
            while restart:
                restart = 0
                # loop over every entry and residue in "occurence"
                for i_ref,o_ref in enumerate(occurence):
                    for j_ref, res_ref in enumerate(o_ref):
                        
                        for i_comp in range(i_ref+1, len(occurence)):
                            for res_comp in occurence[i_comp]:
                                
                                if res_comp == res_ref:
                                    # delete the shorter entry
                                    if len(occurence[i_comp]) > len(occurence[i_ref]):
                                        occurence.pop(i_ref)
                                        restart = 1
                                        break
                                    if len(occurence[i_comp]) <= len(occurence[i_ref]):
                                        occurence.pop(i_comp)
                                        restart = 1
                                        break
                                else:
                                    # if chain id is missing name and resid are enough
                                    if res_comp['chain'] == ' ' or res_ref['chain'] == ' ':
                                        if res_comp['name'] == res_ref['name'] and\
                                                res_comp['resid'] == res_ref['resid']:
                                           
                                            # keep the entry with chain id if possible
                                            if res_comp['chain'] == '':
                                                occurence.pop(i_comp)
                                                restart = 1
                                                break
                                            else:
                                                occurence.pop(i_ref)
                                                restart = 1
                                                break

                                            
                                        
                            
                            if restart:
                                break
                        if restart:
                            break
                    if restart:
                        break    
            # while
            
            


            print "filtered:"
            for o in occurence:
                print "next:"
                for res in o:
                    print res['name'] + ' : ' \
                            + res['chain'] + ' : ' + str(res['resid'])
        
            
            return occurence
            
            



# ToDo:
# kann der Lingand an das Protein gebunden sein?
#    -> Warnung ausgeben!
#
# Kann man irgendwie sicherstellen, dass der Ligand auch in der Bindungstasche 
# sitzt und nicht irgendwo anders cokristallisiert wurde?



if __name__ == '__main__':

    filename = '/scratch/scratch/pdb/pdb_bio/hi/1hiq.pdb1.gz'

    pdb = pdb_from_biopython()
    pdb.parse_pdb_file(filename)
    pdb.merge_models()
    pdb.write_pdb_file('/user/tmeyer/temp/merged1.pdb')
