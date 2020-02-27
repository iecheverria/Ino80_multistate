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

include_sys2 = True
include_sys3 = True
include_Occams = True
include_XLs_spec1 = True
include_XLs_spec2 = True
include_XLs_spec3 = True
include_Occams_12 = True
include_Occams_23 = True


data_dir = 'data/'
top_dir = './'

mdl = IMP.Model()

###############################
# State 1
###############################
top_spec1 = top_dir+'topology_small.dat'
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
output.init_rmf("ini_state1.rmf3", [hier_S1])
output.write_rmf("ini_state1.rmf3")
output.close_rmf("ini_state1.rmf3")

###############################
# State 2
###############################
if include_sys2:
    top_spec2 = top_dir+'topology_small.dat'
    reader_spec2 = IMP.pmi.topology.TopologyReader(top_spec2,
                                                   pdb_dir = data_dir+'/yeast',
                                                   fasta_dir = data_dir+'/yeast')
                                                   

    bs_spec2 = IMP.pmi.macros.BuildSystem(mdl,
                                          resolutions=[1,10])
    bs_spec2.add_state(reader_spec2,
		       keep_chain_id=True)

    hier_S2,  dof_S2 = bs_spec2.execute_macro(max_rb_trans=2.8,
                                              max_rb_rot=0.08)
    mols_S2 = bs_spec2.get_molecules()[0]


output = IMP.pmi.output.Output()
output.init_rmf("ini_state2.rmf3", [hier_S2])
output.write_rmf("ini_state2.rmf3")
output.close_rmf("ini_state2.rmf3")

###############################
# State 3
###############################
if include_sys3:
    top_spec3 = top_dir+'topology_small.dat'
    reader_spec3 = IMP.pmi.topology.TopologyReader(top_spec3,
                                                   pdb_dir = data_dir+'/yeast',
                                                   fasta_dir = data_dir+'/yeast')
                                                   

    bs_spec3 = IMP.pmi.macros.BuildSystem(mdl,
                                          resolutions=[1,10])
    bs_spec3.add_state(reader_spec3,
		       keep_chain_id=True)

    hier_S3,  dof_S3 = bs_spec3.execute_macro(max_rb_trans=2.8,
                                              max_rb_rot=0.08)
    mols_S3 = bs_spec3.get_molecules()[0]


output = IMP.pmi.output.Output()
output.init_rmf("ini_state3.rmf3", [hier_S2])
output.write_rmf("ini_state3.rmf3")
output.close_rmf("ini_state3.rmf3")


##############################
# Combined hierarchy
##############################
p = IMP.Particle(mdl)
hier_all = IMP.atom.Hierarchy.setup_particle(p)
print(hier_S1.get_children()[0])

S1  = hier_S1.get_children()[0]
S1.set_name('State_1')
S2  = hier_S2.get_children()[0]
S2.set_name('State_2')
S3  = hier_S3.get_children()[0]
S3.set_name('State_3')
hier_all.add_child(S1)
hier_all.add_child(S2)
hier_all.add_child(S3)
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

all_mols = [mols_S1, mols_S2, mols_S3]

##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []
display_restraints = []


crs = []
for i, mols in enumerate(all_mols):
    for molname in mols:
        for mol in mols[molname]:
            copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
            cr.set_label(str(i)+mol.get_name()+'.'+str(copy_n))
            cr.add_to_model()
            output_objects.append(cr)
            crs.append(cr)
        
##############################
# Excluded Volume
##############################
for i, mols in enumerate(all_mols):
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols.values(),
                                                              resolution=10)
    evr.add_to_model()
    evr.set_label(str(i)+'_EV')
    evr.set_weight(1.0)
    output_objects.append(evr)


##############################
# XLs
##############################

cldbkc=CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("protein1")
cldbkc.set_protein2_key("protein2")
cldbkc.set_residue1_key("residue1")
cldbkc.set_residue2_key("residue2")
#cldbkc.set_unique_id_key("UniqueID")
cldbkc.set_psi_key("Psi")


if include_XLs_spec1:
    cldb_1=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb_1.create_set_from_file(data_dir+"crosslinks/LYS_sol.csv")

    xl1_1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier_S1,
                                                                                  CrossLinkDataBase=cldb_1,
                                                                                  resolution=1,
                                                                                  length=26.0,
                                                                                  slope=0.04,
                                                                                  label="XLs_state1")


    xl1_1.add_to_model()
    xl1_1.set_weight(1.0)
    rmf_restraints.append(xl1_1)
    output_objects.append(xl1_1)
    dof_S1.get_nuisances_from_restraint(xl1_1)
    rmf_restraints.append(xl1_1)
    xl1_1.set_psi_is_sampled(True)


if include_XLs_spec2:
    cldb_2=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb_2.create_set_from_file(data_dir+"crosslinks/LYS_DNA.csv")

    xl1_2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier_S2,
                                                                                  CrossLinkDataBase=cldb_2,
                                                                                  resolution=1,
                                                                                  length=26.0,
                                                                                  slope=0.04,
                                                                                  label="XLs_state2")


    xl1_2.add_to_model()
    xl1_2.set_weight(1.0)
    rmf_restraints.append(xl1_2)
    output_objects.append(xl1_2)
    dof_S1.get_nuisances_from_restraint(xl1_2)
    rmf_restraints.append(xl1_2)
    xl1_2.set_psi_is_sampled(True)

if include_XLs_spec3:
    cldb_3=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb_3.create_set_from_file(data_dir+"crosslinks/LYS_NUCL.csv")

    xl1_3 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier_S3,
                                                                                  CrossLinkDataBase=cldb_3,
                                                                                  resolution=1,
                                                                                  length=26.0,
                                                                                  slope=0.04,
                                                                                  label="XLs_state3")


    xl1_3.add_to_model()
    xl1_3.set_weight(1.0)
    rmf_restraints.append(xl1_3)
    output_objects.append(xl1_3)
    dof_S1.get_nuisances_from_restraint(xl1_3)
    rmf_restraints.append(xl1_3)
    xl1_2.set_psi_is_sampled(True)

##############################
# Occams Spring restaint
##############################

if include_Occams_12:
    occ1 = IMP.pmi.restraints.occams.OccamsRestraint(hier_S1,
                                                    hier_S2,
                                                    top_dir+'equiv_assig_S12.dat',
                                                    data_dir+'alns/',
                                                    sample_sys_1 = True,
                                                    sigma_init=10.0,
                                                    slope = 0.0005,
                                                    psi_nuisances_are_optimized=True,
                                                    sigma_nuisances_are_optimized=True)
    occ1.add_to_model()
    occ1.set_weight(10.0)
    output_objects.append(occ1)
    dof_S1.get_nuisances_from_restraint(occ1)
    dof_S2.get_nuisances_from_restraint(occ1)
    occ1.write_distances('dist_spec12')

if include_Occams_23:
    occ2 = IMP.pmi.restraints.occams.OccamsRestraint(hier_S2,
                                                    hier_S3,
                                                    top_dir+'equiv_assig_S23.dat',
                                                    data_dir+'alns/',
                                                    sample_sys_1 = True,
                                                    sigma_init=10.0,
                                                    slope = 0.0005,
                                                    psi_nuisances_are_optimized=True,
                                                    sigma_nuisances_are_optimized=True)
    occ2.add_to_model()
    occ2.set_weight(10.0)
    output_objects.append(occ2)
    dof_S2.get_nuisances_from_restraint(occ2)
    dof_S3.get_nuisances_from_restraint(occ2)
    occ1.write_distances('dist_spec23')

##############################
# Shuffle
##############################    
IMP.pmi.tools.shuffle_configuration(hier_S1)
IMP.pmi.tools.shuffle_configuration(hier_S2)
IMP.pmi.tools.shuffle_configuration(hier_S3)

                                    
dof_S1.optimize_flexible_beads(200)
dof_S2.optimize_flexible_beads(200)
dof_S3.optimize_flexible_beads(200)

dof_all = dof_S1.get_movers() + dof_S2.get_movers() + dof_S3.get_movers()

############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling
#taskid = int(sys.argv[1])
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier_all,                            
                                    crosslink_restraints=rmf_restraints,           
                                    monte_carlo_sample_objects=dof_all,   
                                    replica_exchange_maximum_temperature=5.0,
                                    global_output_directory='output',
                                    output_objects=output_objects,
                                    monte_carlo_steps=20,
                                    number_of_frames=50000,
                                    number_of_best_scoring_models=0)

rex.execute_macro()
exit()


















