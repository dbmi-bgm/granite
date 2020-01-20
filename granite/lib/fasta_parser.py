#!/usr/bin/env python


######################################################################################
##	fasta_parser.py
#
#	Author: Michele Berselli
#		University of Padova
#		berselli.michele@gmail.com
#
##	LICENSE:
#		Copyright (C) 2018  Michele Berselli
#
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation.
#
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################


import sys, argparse, os
from bitarray import bitarray


########################################
####		CLASS FastaHandler     ####
########################################
class FastaHandler(object):
	''' a class to implement an object to handle multiple fasta sequences '''

	########################################
	####		CLASS FastaSequence    ####
	########################################
	class FastaSequence(object):
		''' a class to implement an object representing a fasta sequence '''

		## FUNCTIONS ##
		def __init__(self, header, sequence):
			self.__header = header
			self.__sequence = sequence
		#end def __init__

		def get_header(self):
			return self.__header
		#end def get_header

		def get_sequence(self):
			if type(self.__sequence) == bitarray:
				return self.__bitarray_to_seq()
			else:
				return self.__sequence
			#end if
		#end def get_seq()

		def get_bitarray(self):
			if type(self.__sequence) != bitarray:
				raise ValueError('Sequence is not encoded as a bitarray, use get_sequence() instead!')
			else:
				return self.__sequence
			#end if
		#end def get_bitarray

		def __bitarray_to_seq(self):
			s, i, max_i = '', 0, len(self.__sequence) - 2
			encode = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
			while i <= max_i:
				s += encode[str(int(self.__sequence[i])) + str(int(self.__sequence[i + 1]))]
				i += 2
			#end while
			return s
		#end def __bitarray_to_seq

		def write_sequence(self, fo, max_char_per_line=0):
			fo.write('>' + self.__header + '\n')
			i, sequence = 0, self.get_sequence()
			if not max_char_per_line:
				fo.write(sequence + '\n')
			else:
				max_i = len(sequence) - max_char_per_line
				while i <= max_i:
					fo.write(sequence[i:i + max_char_per_line] + '\n')
					i += max_char_per_line
				#end while
				if sequence[i:]: fo.write(sequence[i:] + '\n')
			#end if
		#def write_sequence

	#end class FastaSequence

	## FUNCTIONS ##
	def __init__(self):
		self.__sequences = []
	#end def __init__

	def parse(self, inputfile):
		''' adds fasta sequences as FastaSequence objects '''
		header, sequence = None, []
		fi = self.__open(inputfile)
		for line in fi:
			if line.startswith('>'):
				if header: self.__sequences.append(self.FastaSequence(header, ''.join(sequence)))
				header, sequence = line.rstrip()[1:], []
			else:
				sequence.append(line.rstrip())
			#end if
		#end for
		if header: self.__sequences.append(self.FastaSequence(header, ''.join(sequence)))
		fi.close()
	#end def parse

	def parse_binary(self, inputfile):
		''' adds fasta sequences binary encoded as FastaSequence objects, 
		works only for upper or lower canonical bases A, C, T, G '''
		encode = {
				'A': '00', 'a': '00',
				'C': '01', 'c': '01',
				'G': '10', 'g': '10',
				'T': '11', 't': '11'
				}
		header, sequence = None, bitarray()
		fi = self.__open(inputfile)
		for line in fi:
			if line.startswith('>'):
				if header: self.__sequences.append(self.FastaSequence(header, sequence))
				header, sequence = line.rstrip()[1:], bitarray()
			else:
				try:
					for c in line.rstrip(): sequence.extend(encode[c])
				except:
					self.__sequences = []
					raise ValueError('Non-canonical base found in sequence ' + header + '!')
				#end try
			#end if
		#end for
		if header: self.__sequences.append(self.FastaSequence(header, sequence))
		fi.close()
	#end def parse_binary

	def parse_generator(self, inputfile):
		''' create a generator of the sequences, 
		each iteration return a (header, sequence) tuple '''
		header, sequence = None, []
		fi = self.__open(inputfile)
		for line in fi:
			if line.startswith('>'):
				if header: yield (header, ''.join(sequence))
				header, sequence = line.rstrip()[1:], []
			else:
				sequence.append(line.rstrip())
			#end if
		#end for
		if header: yield (header, ''.join(sequence))
		fi.close()
	#end def parse_generator

	def __open(self, inputfile):
		''' check the existance of inputfile and opens it if exists '''
		if os.path.isfile(inputfile):
			return open(inputfile, 'r')
		else:
			raise ValueError('Input file is missing!')
		#end if
	#end def __open

	def iter_sequences(self):
		''' create a generator of the stored sequences,
		each iteration return a FastaSequence object '''
		return iter(self.__sequences)
	#end def iter_sequences

	def get_sequences(self):
		''' return the list of the stored FastaSequence objects '''
		return self.__sequences
	#end def get_sequences

	def write_sequences(self, fo, max_char_per_line=0):
		for s in self.__sequences:
			s.write_sequence(fo, max_char_per_line)
		#end for
	#end def write_sequences

#end class FastaHandler

