#26 novembre 2017
#Create a freq and opacity file for a specific set of points from a fortran data file

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import glob
import multiprocessing as mp

h_Js = 6.626070040e-34
h_eVs = 4.135667662e-15
list_element = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn']


#Open fortran file and return a clean list

def cut_fortran_to_list(line, delimiter, debug = 0):
    line_list = []
    temp_str = []
    column_number = 0
    char_location = 0

    for char_number in xrange(len(line)):
        if debug == 1: print "char_number : " + str(char_number) + " column_number :  " + str(column_number) + " Final line index : " + str(char_per_column[column_number] + char_location)
        if debug == 1: print "Line number : " + str(line_number)

        if char_number >= delimiter[column_number] + char_location:
            try:   
                temp_flt = double(''.join(temp_str))
                line_list.append(temp_flt)

            except ValueError:
                line_list.append(''.join(temp_str))

            column_number += 1
            char_location = char_number
            if debug == 1: print "Switching column : " + str(temp_flt) + " :Number of char : " + str(len(temp_str))

            temp_str = []
            temp_flt = 0

        temp_str.append(line[char_number])

    return line_list

def Rydberg_to_freq(value):
	#2 Ryd = 27.211396 eV -- Source : http://greif.geo.berkeley.edu/~driver/conversions.html
	return (double(value) * 13.605698) / h_eVs

def MBarn_to_cm2(value):
	#megabarn	Mb = 10-18 cm2 -- Source : https://en.wikipedia.org/wiki/Barn_(unit)
	return double(value) * 1e-18


def open_Topbase_data (filename):
	#filename = 'Data_xsection.txt'
	data_file = open(filename, 'r')

	data_delimiter = [len('  9.900000E-01'), len(' 6.476E+00')]
	header_delimiter = [len('       1'), len('   1'), len('   1'), len('   200'), len('    1'), len('  -1.00000E+00'), len('     101')]

	data_list = []
	lines = data_file.readlines()

	for line_index in xrange(len(lines)):
		#To skip the header
		if lines[line_index][1] != '=' and lines[line_index][-2] != 'P' and lines[line_index][0] != '<':

			if len(lines[line_index]) > int(sum(data_delimiter)) + 2:
				data_list.append([])
				temp = cut_fortran_to_list(lines[line_index], header_delimiter)
				data_list[-1].append(temp)
				data_list[-1].append([[],[]])

			elif len(lines[line_index]) < int(sum(header_delimiter) + 2):
				temp = cut_fortran_to_list(lines[line_index], data_delimiter)
				
				data_list[-1][1][0].append(Rydberg_to_freq(temp[0]))

				data_list[-1][1][1].append(MBarn_to_cm2(temp[1]))

	data_file.close()
	return data_list


def openFortranData(name, column_size):
	data = genfromtxt(name, delimiter = column_size)
	return data


def linear_interpolation (x_1, x_2, y_1, y_2, x_3):
    #x_3 est le point qui se trouve entre les deux points et dont on cherche le y
    a = (y_2 - y_1) / (x_2 - x_1)
    b = y_1 - (a * x_1)
    
    return (a*x_3) + b
    
    
def createFreq_list(original_data_list, weight_data):

	data_list = []
	freq_list = []

	
	for j in range(len(weight_data)):
		freq_list.append(weight_data[j][0])


	for index in xrange(len(original_data_list)):
		data_list.append([original_data_list[index][0], [[],[]]])

		data_list[-1][1][0] = freq_list
		for i in xrange(len(freq_list)):
			data_list[-1][1][1].append(-1)

		for element_index in xrange(len(original_data_list[index][1][0]) - 1):
			try:
				x1 = original_data_list[index][1][0][element_index]
				x2 = original_data_list[index][1][0][element_index + 1]
				y1 = original_data_list[index][1][1][element_index]
				y2 = original_data_list[index][1][1][element_index + 1]
				
				for freq_index in xrange(len(freq_list)):
					if freq_list[freq_index] >= x1 and freq_list[freq_index] < x2:
						x3 = freq_list[freq_index]
						temp = linear_interpolation(x1, x2, y1, y2, x3)

						#data_list[-1][1][0].append(x3)
						data_list[-1][1][1][freq_index] = temp

					if freq_list[freq_index] < x1 or freq_list[freq_index] > x2:
						pass

			except IndexError:
				print "IndexError in linear_interpolation loop"
				break

	return data_list


#Write all the big list to files

def write_xsection_data(workfile, one_array, header, columns = ['x','y'], filetype = ".txt", opencommand = 'w'):
    #For column based arrays of one or multiple arrays. Array of column arrays -> All x in one column and all y in the next
       
	index_line = 0
	index_column = 0
	file_num = 0    
    
	#One array of multiple columns
	if len(np.shape(one_array)) == 2 or len(np.shape(one_array)) == 1:
		try :
			if workfile[-4:] == ".txt" or workfile[-4:] == filetype:
				f = open(workfile, opencommand)       
			if workfile[-4:] != ".txt" or workfile[-4:] != filetype:
				f = open(workfile + filetype, opencommand)
		except IOError:
			if workfile[-4:] == ".txt" or workfile[-4:] == filetype:
				f = open(workfile,'w')       
			if workfile[-4:] != ".txt" or workfile[-4:] != filetype:
				f = open(workfile + filetype,'w')


		f.write('{:6d}'.format(int(header[0])))
		f.write('\t')
		f.write('{:3d}'.format(int(header[1])))
		f.write('\t')
		f.write('{:3d}'.format(int(header[2])))
		f.write('\t')
		f.write('{:6d}'.format(int(header[3])))
		f.write('\t')
		f.write('{:6d}'.format(int(header[4])))
		f.write('\t')
		f.write('%1.10e'%double(header[5]))
		f.write('\t')
		f.write('{:6d}'.format(int(header[6])))
		f.write('\n')

		#for i in header[1:]:
			#f.write('\t')
			#f.write(str(int(i)))
		#f.write('\n')

		f.write(columns[0])
		for column_name in columns[1:]:
			f.write("\t")
			f.write(column_name)
		f.write("\n")

		while index_line < len(one_array[0]):
			while index_column < len(one_array):
				temp_str = '%1.15e'%one_array[index_column][index_line]
				f.write(temp_str.replace('e', 'D'))
				#f.write("{:15.10E}".format(Decimal(str(one_array[index_column][index_line]))))

				if index_column + 1 < len(one_array):
					f.write("\t")

				index_column += 1

			index_column = 0
			f.write("\n")
			index_line += 1

		index_line = 0
		f.close() 
		return
	return



def write_Topbase_data(one_element_data_list):

	for element in one_element_data_list:
		if len(element[1][0]) > 0:
			element_name = list_element[int(element[0][1]) - 1]

			if os.path.exists('xsection_data/' + str(element_name)):
				filename = 'xsection_data/' + str(element_name) + '/Z' + str(int(element[0][1])) + 'E' + str(int(element[0][2])) + '.txt'
				#filename = 'TB_data/TB_xsection_' + str(element_name) + '_i=' + str(element[0][0]) + '_Ne=' + str(element[0][2]) + '.txt'
				print "Writing file : " + filename
				write_xsection_data(filename, element[1], element[0], columns = ['Freq', 'xsec'], opencommand = 'a')

			else:
				os.mkdir('xsection_data/' + str(element_name))
				filename = 'xsection_data/' + str(element_name) + '/Z' + str(int(element[0][1])) + 'E' + str(int(element[0][2])) + '.txt'
				#filename = 'TB_data/TB_xsection_' + str(element_name) + '_i=' + str(element[0][0]) + '_Ne=' + str(element[0][2]) + '.txt'
				print "Writing file : " + filename
				write_xsection_data(filename, element[1], element[0], columns = ['Freq (Hz)', 'xsection (cm2)'], opencommand = 'a')

	return

def work_on_one_file(filename):
	print "Working on : " + str(filename)
	original_data_list = open_Topbase_data(filename)
	freq_weight_file = openFortranData("newgrid.txt", [len("    0.5995849160E+17"), len("    0.5878283490E+15"), len("    0.5000000000E+02")])
	freq_data_list = createFreq_list(original_data_list, freq_weight_file)
	return freq_data_list

if os.path.exists('xsection_data/'):
	pass
else:
	os.mkdir('xsection_data/')


list_fichier = glob.glob('TB_data/Topbase_*')
os.system('rm -r ./xsection_data/*')
work_pool = mp.Pool(processes = 8)
result_list = []

result_list = work_pool.map_async(work_on_one_file, list_fichier)

"""
for filename in list_fichier:
	print "Working on : " + str(filename)
	original_data_list = open_Topbase_data(filename)
	freq_weight_file = openFortranData("newgrid.txt", [len("    0.5995849160E+17"), len("    0.5878283490E+15"), len("    0.5000000000E+02")])
	freq_data_list = createFreq_list(original_data_list, freq_weight_file)
	write_Topbase_data(freq_data_list)
	"""


for freq_data_list in result_list.get():
	write_Topbase_data(freq_data_list)



work_pool.close()
work_pool.join()
