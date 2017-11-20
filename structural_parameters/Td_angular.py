# -*- coding: utf-8 -*-
"""
Made by Ernesto Martinez-Baez 11-09-17

Code to calculate Td angular structure parameter on a defined core structure
from an xyz file.
"""
import sys
import numpy as np
import math

with open(sys.argv[2], 'w') as f_out:
    
    f_out.write(str('##Tetrahedral parameter \n##Td_param \
                     angles [1 2 3 4 5 6] \n'))


    f_xyz = open(sys.argv[1], 'r')
    flines_xyz = f_xyz.readlines()
    center_element = "Al"
    vertex_elements = "O"
    N = 282
#    raw_data = np.loadtxt('new_2xyz.xyz', dtype=np.ndarray)
    t = 0
    for k in range(0,len(flines_xyz),N+2):
        
        data_lines = flines_xyz[k:k+N+2]
        data = np.loadtxt(data_lines, dtype=np.ndarray, skiprows=2)
        coord_data = []
        
        def center_atom(a): 
            center_xyz = []  
            stripped_xyz = []
            for i in xrange(len(a)):
              stripped_xyz = a[i].split()
              if stripped_xyz[0] == center_element :
                center_xyz = stripped_xyz[0:4]
                center_xyz.insert(len(center_xyz),i-2)
            return(center_xyz)    
        
        
        
        def distance_matrix(a):
            D = []
        
            dxx = []
            dyy = []
            dzz = []
        
            j=0
            for i in xrange(len(a)):
              if a[i][0] == vertex_elements :
        
                dxx.append(math.pow(math.fabs(
                float(center_atom(data_lines)[1])) - float(float(a[i][1])),2))
                dyy.append(math.pow(math.fabs(
                float(float(center_atom(data_lines)[2])) - float(float(a[i][2]))),2))
                dzz.append(math.pow(math.fabs(
                float(float(center_atom(data_lines)[3])) - float(float(a[i][3]))),2))
                D.append([math.pow(float(dxx[j]) + float(dyy[j]) + float(dzz[j]),0.5),i])        
                
                j += 1
              else:
                  continue
            D_array = np.asarray(D, dtype=float)  
            b = D_array[np.lexsort((D_array[:,1],D_array[:,0]))]   # https://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.lexsort.html
            #Will sort based on the D_array[:,0] which are the first elements of each 2D array and if equal sorting will look at 2nd elements sorting
            return(b)    
        
        def Td_parameter(a):
            cos_ij = 0
            q = 0
            Td_q = 0
            angles = []
            for i in range(0,3):
                for j in range(i+1,4):
                    cx2 = math.pow(float(a[int(distance_matrix(a)[i,1]),1])-
                                  float(a[int(distance_matrix(a)[j,1]),1]),2)
                    cy2 = math.pow(float(a[int(distance_matrix(a)[i,1]),2])-
                                  float(a[int(distance_matrix(a)[j,1]),2]),2)
                    cz2 = math.pow(float(a[int(distance_matrix(a)[i,1]),3])-
                                  float(a[int(distance_matrix(a)[j,1]),3]),2)
                    c2 = cx2 + cy2 +cz2
                    
                    cos_ij = float((math.pow(float(distance_matrix(a)[i,0]),2)+math.pow(float(distance_matrix(a)[j,0]),2)
                             -c2))/float((2*float(distance_matrix(a)[i,0])*float(distance_matrix(a)[j,0])))
                    
                    angles.append(math.degrees(math.acos(cos_ij)))
                    #print (cos_ij, angle)
                    
                    q += math.pow((1.0/3.0 + cos_ij),2)
        
            Td_q = 1.0-float(3*q)/float(8)
             
            return(Td_q, angles)
        
        output = Td_parameter(data)
#        print output
        
        f_out.write(str(t) + '   ' +
                    str(output[0]) + '   ' + 
                    str(output[1][0])[:7] + '   ' +
                    str(output[1][1])[:7] + '   ' +
                    str(output[1][2])[:7] + '   ' +
                    str(output[1][3])[:7] + '   ' +
                    str(output[1][4])[:7] + '   ' +
                    str(output[1][5])[:7] + '   ' + '\n')
        t = t + 1
