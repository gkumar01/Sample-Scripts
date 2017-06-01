#!/bin/env python

#
# A program to compare PKA conformations
#

import math, sys, string, os
sys.path.append('/usr/local/dislin/python')

import Protool
import Numeric
#import dislin_driver_matrices

def generate_matrix(filename):
    """Generate the full pairwise distance matrix"""
    #
    # First check if this matrix has already been calculated
    #
    matrix_name=filename+'_matrix'
    import os
    if os.path.isfile(matrix_name):
        print 'Loading matrix from file'
        try:
            fd=open(matrix_name)
            import pickle
            matrix=pickle.load(fd)
            fd.close()
            return matrix
        except:
            print 'Loading failed - recalculating matrix'
    #
    # No it was not calculated..
    #
    X1=Protool.structureIO()
    X1.readpdb(filename)
    residues=X1.residues.keys()
    residues.sort()
    #
    # Get the atoms to calculate for
    #
    atoms_tmp=X1.atoms.keys()
    atoms_tmp.sort()
    atoms_tmp=atoms_tmp[:2508] ###512

    include=['TPO']
    
    for residue in residues:
        resname=X1.resname(residue)
        if resname in include:
            atoms_tmp=X1.atoms.keys()
            atoms_tmp.sort()
            atoms_tmp=atoms_tmp[:2512]
    #
    # We exclude all backbone atoms except CA
    #
    exclude=['C','O','N','H',"O'", 'P', 'O1P', 'O2P', 'O3P']
    atoms=[]
    for atom in atoms_tmp:
        atom_name=X1.atname(atom)
        if atom_name in exclude:
            pass
        else:
            atoms.append(atom)
    #
    # Sort them
    #
    atoms.sort()
    #print len(atoms)


    distances=Numeric.zeros((len(atoms),len(atoms)),Numeric.Float)

    #
    # Start calculating the distances
    #
    counter=0
    row=0
    import sys
    print 'starting    ',
    for atom1 in atoms:
        counter=counter+1
        text='%%done %5.1f ' %(float(counter)/float(len(atoms))*100.0)
        back=(len(text)+1)*'\b'
        rt=back+text
        print rt,
        sys.stdout.flush()
        column=row
        for atom2 in atoms[row:]:
            distances[row][column]=X1.dist(atom1,atom2)
            distances[column][row]=distances[row][column]
            column=column+1
        #
        # Update the row counter
        #
        row=row+1
    print
    #
    # Save the matrix
    #
    print 'Saving matrix to file %s' %matrix_name
    fd=open(matrix_name,'w')
    import pickle
    pickle.dump(distances,fd)
    fd.close()
    #
    # All done
    #
    
    return distances

#
# -----
#

def main(file1,file2,file3):
    import Numeric
    print 'Generating matrix1'
    
    matrix1=generate_matrix(file1) ###1bkx.corrected.pdb     1bkx.small.pdb
    print 'Generating matrix2'
    matrix2=generate_matrix(file2) ###1atp.corrected.pdb     1atp.small.pdb
    #
    # Calculate matrix A
    #
    matrix_shape1 = Numeric.shape(matrix1)
    print "matrix_shape1", matrix_shape1
    matrix_shape2 = Numeric.shape(matrix2)
    print "matrix_shape2", matrix_shape2

    matrix_A = abs(matrix1-matrix2)
    matrix_DM_A=Numeric.array(matrix_A)
    #
    # Generate matrix 3
    #

    print 'Generating matrix3'
    matrix3=generate_matrix(file3) ###1j3h.corrected.pdb      1j3h.small.pdb
    #
    # Calculate matrix B
    #
    matrix_B = abs(matrix3-matrix2)
    #
    # What happens here?
    #
    matrix_DM_B=Numeric.array(matrix_B)
    #
    #Because the last value in matrix_A and matrix_B is 0, and to make matrix_C I have to divide A/B =C than I had to change the >0.01 value in m_B to -999
    #

    for row in range(len(matrix_B)): 
        for column in range(len(matrix_B)): 
            #print matrix_B[row][column]
            if abs(matrix_DM_B[row][column])<0.01:
                matrix_DM_B[row][column]=-999
                #print row,column,'is 1'
                #print '"1.00" occurs %d times in matrix_B' % \


    #
    #Only take disances from matrix_A and matrix_B which are bigger than 0.5A and InputFile[1] is DM_A.txt and InputFile[2] is DM_B.txt
    #

    import Numeric
    new_matrix_B=Numeric.zeros((len(matrix_DM_A),len(matrix_DM_A)),Numeric.Float)
    matrix_C = Numeric.zeros((len(matrix_DM_A),len(matrix_DM_A)),Numeric.Float) ###Numeric.Float gives all the values numers wiht "extra digits" instead of only rounding numbers!s



    for row in range(len(matrix_DM_A)):
         for column in range(len(matrix_DM_A[row])):
             if (matrix_DM_B[row][column] >= 0.5):
                 matrix_C[row][column] = matrix_DM_A[row][column]/matrix_DM_B[row][column]
                 new_matrix_B[row][column]=matrix_DM_B[row][column]
             else:
                 matrix_C[row][column]=-999
                 new_matrix_B[row][column]=-999
    
    #print "matrix_C", matrix_C   ######

             #if abs(matrix_C[row][column])<0.01:
             #    print matrix_DM_A[row][column],matrix_DM_B[row][column],matrix_C[row][column]
    #print "matrix_C", matrix_C
    #print Numeric.shape(matrix_C) ###How this matrix looks like; (x,y)

   

    return matrix_C, new_matrix_B ###if you have a return command before the end of your schript you want to run it will not calculate the rest of your script!!




if __name__=="__main__":
    import sys
    main(sys.argv[1],sys.argv[2],sys.argv[3])

        
