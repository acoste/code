# -*- coding:Utf-8 -*-

#Author : Arthur COSTE
#Creation Date : July 12th, 2012
#Update Date : July 16th, 2012

#Status : Working
# Version : (RC) 1.1

#Purpose : this script is processing DWI headers to remove certain gradient directions 
#	   We have to modify the dimension, the pointed DWI and the gradient directions set
#	   this is usefull while performing Leave One Out to process headers
#	   To realize Leave Multiple Out we can use that script mutliples times

# Input : Original DWI header, New DWI header path, new DWI file name, gradient direction to be removed
#	  Size of header, number of gradient directions, line of dwi file linked in header

# Output : New DWI header

# usage: text_file_processing.py [-h] [-v] [-i I] [-o O] [-dwi DWI] [-grad GRAD] [-sh SH] [-gn GN] [-dl DL] [-dwil DWIL]
# default values to use if LOO : sh = 20, gn = 25, dl = 7, dwil = 17

#-----------------------------------------------------------------------------------------------------------------------

import argparse
import os

#parsing command line arguments
parser = argparse.ArgumentParser(epilog="Default Values for nhdr file : sh=20, gn=25, dl=7, dwil=17")
parser.add_argument("-v","--verbose", help="output verbosity", action="store_true")
parser.add_argument("-i", type=str,help="Input hearder file path")
parser.add_argument("-o", type=str,help="Output modified header path")
parser.add_argument("-dwi", type=str, help="Associated DWI to the modified header")
parser.add_argument("-grad", type=int,help="Gradient Direction to remove")
parser.add_argument("-sh", type=int,help="Header number of lines")
parser.add_argument("-gn", type=int,help="Number of input gradient directions")
parser.add_argument("-dl", type=int,help="Header line containing dimension of file")
parser.add_argument("-dwil", type=int,help="Header line containing dwi pointed by header")
args = parser.parse_args()

#default values : number of grad dir = 25, dimension line = 7, header dimension = 20, dwil = 17

# parameters to set and handling errors
if args.i: source_file = args.i
else: print "ERROR : Input file not specified"
if args.o: new_file = args.o
else: print "ERROR : Output file not specified"
if args.dwi: dwi_file = args.dwi
else: print "ERROR : Associated DWI not specified"
if args.grad: remove_gd = args.grad 
else: print "ERROR : No gradient direction specified"
if args.sh: header_number_of_lines = args.sh 
else: print "ERROR : Header Size not specified"
if args.gn: number_of_gd = args.gn
else: print "ERROR : Number of Input Directions not specified"
if args.dl: dimension_line = args.dl
else: print "ERROR : Header Line containing dimension of file"
if args.dwil: dwi_line = args.dwil
else: print "ERROR : DWI line not specified"

# associated parameters
line_grad = header_number_of_lines + remove_gd
end_of_modifications = header_number_of_lines + number_of_gd
old_dimension = number_of_gd+1
new_dimension = number_of_gd
numbering_discontinuity = header_number_of_lines+10
dwi_line_2 = dwi_line-14
encoding_line = dwi_line-4
old_extension = "gzip"
patho, new_extension = os.path.splitext(args.dwi)
#processing

if args.i and args.o and args.grad and args.sh and args.gn and args.dl and args.dwil:

	if remove_gd==0:
		print "ERROR : cannot remove Baseline Image"
	if remove_gd<0:
		print "ERROR : inconsistent designation of Image"
	if remove_gd>0:	

		#open files
		source_file = open(args.i, 'r')
		new_file = open(args.o,'a')
	
		for i, line in enumerate (source_file):

			# modify file dimension
			if i == dwi_line_2:
				if args.verbose: print "new dwi associated to new header: ", dwi_file
				base_content = "content: exists("
				old_content = base_content + str(args.i.replace(".nhdr",".raw")) +str(",0)")
				new_content = base_content + str(args.dwi) +str(",0)")
				if args.verbose: print old_content
				if args.verbose: print new_content
				new_file.write(line.replace(old_content,new_content,1))

			if i == dimension_line:
				if args.verbose: print "new dimension of set is now: ", new_dimension
				new_file.write(line.replace(str(old_dimension),str(new_dimension),1))

			if i == encoding_line:
				if args.verbose: print "new file extension is ", new_extension
				base_encoding = "encoding: "
				old_encoding = base_encoding + str(old_extension) 
				new_encoding = base_encoding + str(new_extension.replace(".","",1)) 
				if args.verbose: print old_encoding
				if args.verbose: print new_encoding
				new_file.write(line.replace(old_encoding,new_encoding,1))

			if i == dwi_line:
				if args.verbose: print "new dwi associated to new header: ", dwi_file
				base_dim = "data file: "
				old_str = base_dim + str(args.i.replace(".nhdr",""))+str(".raw.gz")
				new_str = base_dim + str(args.dwi) 
				new_file.write(line.replace(old_str,new_str,1))

			if i == line_grad:				# we do not copy the associated information
				print "removing gradient :", remove_gd
				pass

			if i < line_grad and  i != dimension_line and i != dwi_line and i != dwi_line_2 and i != encoding_line:		#while inferior at the chosen we don't change anything
				new_file.write(line)

			if i > line_grad and i<=end_of_modifications:	#then we process grad dirs with a small trick for numbers		

				if i > dimension_line and i < numbering_discontinuity:
					base_string_1 = base_string_2 = "_000"
				if i == numbering_discontinuity:
					base_string_1 = "_00"
					base_string_2 = "_000"

				if i > numbering_discontinuity:
					base_string_1=base_string_2= "_00"

				current_gd = base_string_1 + str(i-header_number_of_lines)
				new_gd = base_string_2 + str(i-header_number_of_lines-1)
				new_file.write(line.replace(current_gd,new_gd,1))
				if args.verbose: print "current gradient dir number : ", current_gd
				if args.verbose: print "new gradient dir number : ",new_gd		

			if i > end_of_modifications:			# all fields after the gradient are left identical
				new_file.write(line)
			i = i + 1

		#close files
		source_file.close()
		new_file.close()

else: print "ERROR : MISSING ARGUMENTS"
