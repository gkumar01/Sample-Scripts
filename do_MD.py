#!/bin/env python
#
# $Id: do_MD.py,v 1.3 2006/01/26 17:50:41 unab Exp $
#
# Gaurav's MD driver
#


def main():
    #
    # Get the name of the pdb file
    #
    import sys
    pdbfile=sys.argv[1]
    root=pdbfile[:-3]

    if len(sys.argv)>2:
        topfile=sys.argv[2]
    #
    # Run pdb2gmx
    #
    pdb2gmx='/usr/local/bin/pdb2gmx'
    grofile=root+'gro'
    topologyfile=root+'top'
    fd=open('script','w')
    fd.write(pdb2gmx+' -f '+pdbfile+' -o '+grofile+' -p '+topologyfile+' << EOF \n')
    fd.write('7\n')
    fd.close()
    import os
    os.system('source ./script')
    #
    #Run editconf
    #
    afterbox_grofile='box.gro'
    os.system('/usr/local/bin/editconf  -f ' +grofile+ ' -o ' +afterbox_grofile+ ' -d 0.85  -c ')

    #
    # Run genbox
    #
    water_grofile='b4em.gro'
    os.system('/usr/local/bin/genbox -cp ' +afterbox_grofile+ ' -cs -p ' +topologyfile+  ' -o ' +water_grofile)
    
    #
    # Define the standard em.mdp file
    #
    params={'cpp':'/lib/cpp',
            'define':'-DFLEX_SPC',
            'constraints':'none',
            'integrator':'steep',
            'nsteps':'5000',
            'emtol':'100',
            'emstep':'0.01',
            'nstcomm':'1',
            'ns_type':'grid',
            'rlist':'1',
            'rcoulomb':'1.0',
            'epsilon_r':'1000.0',
            'rvdw':'1.0',
            'Tcoupl':'no',
            'Pcoupl':'no',
            'gen_vel':'no',
            'nstenergy':'1'}
    fd=open('em.mdp','w')
    for key in params.keys():
        fd.write('%s\t\t=  %s\n' %(key,params[key]))
    fd.close()


    #
    #Run grompp
    #
    emfile_mdpfile='em.mdp'
    emout='em_out.mdp'
    structure_file='em.tpr'
    os.system('/usr/local/bin/grompp -f '+emfile_mdpfile+ ' -po ' +emout+ ' -c ' +water_grofile+ ' -o ' +structure_file+ ' -p ' +topologyfile)

    #
    #Run mdrun
    #
    trajectory_file='em.trr'
    after_em_grofile='em.gro'
    after_em_logfile='em.log'
    after_em_edrfile='em.edr'
    os.system('/usr/local/bin/mdrun -s ' +structure_file+ ' -o ' +trajectory_file+ ' -c ' +after_em_grofile+ ' -g ' +after_em_logfile+ ' -e ' +after_em_edrfile)
    #
    # Define the standard pr.mdp file
    #
    params={'cpp':'/lib/cpp',
            'define':'-DPOSRES',
            'include':'-I../top',
            'constraints':'all-bonds',
            'integrator':'md',
            'dt':'0.001',
            'nsteps':'100000',
            'lincs_iter':'2',
            'coulombtype':'PME',
            'tau_p':'1',
            'tc-grps':'Protein OTHER',
            'tau_t':'0.1 0.1',
            'gen_temp':'300',
            'ref_t':'300 300',
            'ref_p':'1',
            'energygrps':'Protein OTHER',
            'Tcoupl':'berendsen',
            'gen_vel':'yes',
            'compressibility':'4.5e-5'}
    fd=open('pr.mdp','w')
    for key in params.keys():
        fd.write('%s\t\t=  %s\n' %(key,params[key]))
    fd.close()


    #
    #Run grompp
    #
    prfile_mdpfile='pr.mdp'
    prout='pr_out.mdp'
    pr_structure_file='pr.tpr'
    os.system('/usr/local/bin/grompp -f '+prfile_mdpfile+ ' -po ' +prout+ ' -c ' +after_em_grofile+ ' -o ' +pr_structure_file+ ' -p ' +topologyfile)

    #
    #Run mdrun
    #
    pr_trajectory_file='pr.trr'
    after_pr_grofile='pr.gro'
    after_pr_logfile='pr.log'
    after_pr_edrfile='pr.edr'
    os.system('/usr/local/bin/mdrun -s ' +pr_structure_file+ ' -o ' +pr_trajectory_file+ ' -c ' +after_pr_grofile+ ' -g ' +after_pr_logfile+ ' -e ' +after_pr_edrfile)
    #
    # Define the standard md.mdp file
    #
    params={'cpp':'/lib/cpp',
            'include':'-I../top',
            'constraints':'none',
            'integrator':'md',
            'dt':'0.001',
            'nsteps':'10000000', ###5000000
            'nstxout':'10000', ###Changed this from 100 (collect data every 0.1ps) to 10000 (collect data every 10ps)
            'lincs_iter':'2',
            'coulombtype':'PME',
            'tau_p':'1',
            'tc-grps':'Protein OTHER',
            'tau_t':'0.1 0.1',
            'gen_temp':'300',
            'ref_t':'300 300',
            'energygrps':'Protein OTHER',
            'Tcoupl':'berendsen',
            'ref_p':'1',
            'gen_vel':'yes',
            'compressibility':'4.5e-5'}
    fd=open('md.mdp','w')
    for key in params.keys():
        fd.write('%s\t\t=  %s\n' %(key,params[key]))
    fd.close()


    #
    #Run grompp
    #
    mdfile_mdpfile='md.mdp'
    mdout='md_out.mdp'
    md_structure_file='md.tpr'
    os.system('/usr/local/bin/grompp -f '+mdfile_mdpfile+ ' -po ' +mdout+ ' -c ' +after_pr_grofile+ ' -o ' +md_structure_file+ ' -p ' +topologyfile)

    #
    #Run mdrun
    #
    md_trajectory_file='md.trr'
    after_md_grofile='md.gro'
    after_md_logfile='md.log'
    after_md_edrfile='md.edr'
    os.system('/usr/local/bin/mdrun -s ' +md_structure_file+ ' -o ' +md_trajectory_file+ ' -c ' +after_md_grofile+ ' -g ' +after_md_logfile+ ' -e ' +after_md_edrfile)
    

    #
    #Make a pdb file after the run
    #

    md_trajectory_file='md.trr'
    md_structure_file='md.tpr'
    pdb_file='after_md.pdb'
    os.system('/usr/local/bin/trjconv -f ' +md_trajectory_file+ ' -o ' +pdb_file+ ' -s ' +md_structure_file)  




    
    return
    

if __name__=="__main__":
    main()

    
