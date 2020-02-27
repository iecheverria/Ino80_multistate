############################################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter
import IMP.pmi.restraints.occams

import random
import numpy as np
import glob
from sys import exit
from sys import argv
import sys
#from occams_restraint import *
#from npc_util import *

###################### SYSTEM SETUP #####################

include_sys1 = True
include_Occams = True
include_XLs_spec1 = True

data_dir = 'data/'
top_dir = './'

mdl = IMP.Model()

###############################
# Species 1
###############################
top_spec1 = top_dir+'topology_02122019.dat'
reader_spec1 = IMP.pmi.topology.TopologyReader(top_spec1,
                                               pdb_dir = data_dir+'yeast/',
                                               fasta_dir = data_dir+'yeast/')

bs_spec1 = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1,10])
bs_spec1.add_state(reader_spec1)


hier_S1, dof_S1 = bs_spec1.execute_macro(max_rb_trans=4.0,
                                          max_rb_rot=1.0)
mols_S1 = bs_spec1.get_molecules()[0]

output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_S1])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

###############################
# Species 2
###############################
if include_sys1:
    top_spec2 = top_dir+'top_6fml_nuc.dat'
    reader_spec2 = IMP.pmi.topology.TopologyReader(top_spec2,
                                                   pdb_dir = data_dir+'/human',
                                                   fasta_dir = data_dir+'/human')
                                                   

    bs_spec2 = IMP.pmi.macros.BuildSystem(mdl,
                                          resolutions=[1,10])
    bs_spec2.add_state(reader_spec2,
		       keep_chain_id=True)

    hier_S2,  dof_S2 = bs_spec2.execute_macro(max_rb_trans=2.8,
                                              max_rb_rot=0.08)
    mols_S2 = bs_spec2.get_molecules()[0]


output = IMP.pmi.output.Output()
output.init_rmf("ini_all_human.rmf3", [hier_S2])
output.write_rmf("ini_all_human.rmf3")
output.close_rmf("ini_all_human.rmf3")


##############################
# Combined hierarchy
##############################
p = IMP.Particle(mdl)
hier_all = IMP.atom.Hierarchy.setup_particle(p)
hier_all.add_child(hier_S1)
hier_all.set_name('System')

states = IMP.atom.get_by_type(hier_all,IMP.atom.STATE_TYPE)


print(type(hier_all) == IMP.atom.Hierarchy)
print(hier_all.get_name(), hier_S1.get_name())
states = IMP.atom.get_by_type(hier_all,IMP.atom.STATE_TYPE)

#Write initial configuration
output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_all])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []
display_restraints = []


crs = []
for molname in mols_S1:
    for mol in mols_S1[molname]:
        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)
        
##############################
# Excluded Volume
##############################
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_S1.values(),
                                                              resolution=10)
evr1.add_to_model()
evr1.set_label('sys0')
evr1.set_weight(1.0)
output_objects.append(evr1)


##############################
# XLs
##############################
if include_XLs_spec1:
    cldbkc=CrossLinkDataBaseKeywordsConverter()
    cldbkc.set_protein1_key("protein1")
    cldbkc.set_protein2_key("protein2")
    cldbkc.set_residue1_key("residue1")
    cldbkc.set_residue2_key("residue2")
    #cldbkc.set_unique_id_key("UniqueID")
    cldbkc.set_psi_key("Psi")
    
    cldb_1=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb_1.create_set_from_file(data_dir+"crosslinks/LYS_NUCL.csv")

    xl1_1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier_S1,
                                                                                  CrossLinkDataBase=cldb_1,
                                                                                  resolution=1,
                                                                                  length=26.0,
                                                                                  slope=0.04,
                                                                                  label="XLs_NUC")


    xl1_1.add_to_model()
    xl1_1.set_weight(1.0)
    rmf_restraints.append(xl1_1)
    output_objects.append(xl1_1)
    dof_S1.get_nuisances_from_restraint(xl1_1)
    rmf_restraints.append(xl1_1)
    xl1_1.set_psi_is_sampled(True)


##############################
# Occams Spring restaint
##############################

if include_Occams:
    occ = IMP.pmi.restraints.occams.OccamsRestraint(hier_S1,
                                                    hier_S2,
                                                    top_dir+'equiv_assig_Ino80_em.dat',
                                                    data_dir+'/alns/',
                                                    sample_sys_1 = False,
                                                    sigma_init=10.0,
                                                    slope = 0.0005,
                                                    psi_nuisances_are_optimized=True,
                                                    sigma_nuisances_are_optimized=True)
    occ.add_to_model()
    occ.set_weight(15.0)
    output_objects.append(occ)
    dof_S1.get_nuisances_from_restraint(occ)
    print(occ.get_output())
    occ.write_distances('dist_nucleosome')

##############################
# Shuffle
##############################    
IMP.pmi.tools.shuffle_configuration(hier_S1)
                                    #max_translation=200,
                                    #bounding_box=((-80, -80, -80), (80, 80, 80)))
dof_S1.optimize_flexible_beads(200)

############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling
#taskid = int(sys.argv[1])
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier_S1,                           
                                    crosslink_restraints=rmf_restraints,          
                                    monte_carlo_sample_objects=dof_S1.get_movers(),   
                                    replica_exchange_maximum_temperature=5.0,
                                    global_output_directory='output',
                                    output_objects=output_objects,
                                    monte_carlo_steps=20,
                                    number_of_frames=50000,
                                    number_of_best_scoring_models=0)

rex.execute_macro()
exit()


















