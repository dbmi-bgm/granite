#!/usr/bin/env python

#################################################################
#
#	novoCaller 2
#		contact: berselli.michele@gmail.com
#
#################################################################


#################################################################
#
#	LIBRARIES
#
#################################################################
import pysam
import numpy as np
import sys, os
import argparse
import re


#################################################################
#
#	FUNCTIONS
#
#################################################################
#################################################################
#	GT_ordering_alternate (original)
#################################################################
def GT_ordering_alternate(ALT_count):
	''' '''
	combos=(ALT_count+1)*(ALT_count+2)//2
	ordering=np.empty([combos,2])
	count=0
	for a1 in range(0,ALT_count+1):
		for a2 in range(a1,ALT_count+1):
			ordering[count,0]=a1
			ordering[count,1]=a2
			count=count+1
		#end for
	#end for
	return ordering
#end def

#################################################################
#	row_gen (original)
#################################################################
def row_gen(GT1,GT2,alt_count,mut_rate):
	''' '''
	N=alt_count
	combos=(N+1)*(N+2)//2
	row=np.zeros(combos)
	count=0
	for a1 in range(N+1):
		for a2 in range(N+1):
			for a3 in range(N+1):
				for a4 in range(N+1):
					P=1.0
					if a1==GT1[0]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					#end if
					if a2==GT1[1]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					#end if
					if a3==GT2[0]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					#end if
					if a4==GT2[1]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					#end if
					count+=1

					for b1 in [a1,a2]:
						for b2 in [a3,a4]:
							gt_work=np.sort([b1,b2])
							index=(2*N+3-gt_work[0])*gt_work[0]//2+gt_work[1]-gt_work[0]
							row[index]=row[index]+0.25*P
						#end for
					#end for
				#end for
			#end for
		#end for
	#end for

	return row
#end def

#################################################################
#	table_gen (original)
#################################################################
def table_gen(alt_count,mut_rate):
	''' '''
	N=alt_count
	II_prev=-1

	combos=(N+1)*(N+2)//2
	table=np.zeros([combos**2,combos])
	for a1 in range(N+1):
		for a2 in range(a1,N+1):
			for a3 in range(N+1):
				for a4 in range(a3,N+1):
					GT1=[a1,a2]
					GT2=[a3,a4]
					I1=(2*N+3-GT1[0])*GT1[0]//2+GT1[1]-GT1[0]
					I2=(2*N+3-GT2[0])*GT2[0]//2+GT2[1]-GT2[0]
					II=I1*combos+I2

					if II<=II_prev:
						sys.exit("error in II calc!!!\n")
					#end if

					row=row_gen(GT1,GT2,alt_count,mut_rate)
					table[II,:]=row
				#end for
			#end for
		#end for
	#end for

	return table
#end def

#################################################################
#	GT_likelihood_wrt_allele_calc (original)
#################################################################
def GT_likelihood_wrt_allele_calc(ALT_count):
	''' '''
	ordering=GT_ordering_alternate(ALT_count)

	combos=(ALT_count+1)*(ALT_count+2)//2
	table=np.zeros([combos,ALT_count+1])*1.

	for i in range(combos):
		a1=int(ordering[i,0]) # casting to int to be used as indexes
		a2=int(ordering[i,1])
		table[i,a1]=table[i,a1]+0.5
		table[i,a2]=table[i,a2]+0.5
	#end for

	return table
#end def

#################################################################
#	encode_chrom
#################################################################
def encode_chrom(chrom):
	''' convert chr IDs to numbers to be used by other functions '''
	encode = {'M': 0, 'MT': 0, 'X': 23, 'Y': 24}
	chrom_repl = chrom.replace('chr', '')
	if chrom_repl in encode:
		return encode[chrom_repl]
	else:
		try:
			return int(chrom_repl)
		except Exception:
			return 25
		#end try
	#end if

#end def

#################################################################
#	buffering_bams
#################################################################
def buffering_bams(bams_infofile):
	''' return a list containing reading buffers to bam files,
	return also a list with the corrisponding IDs associated
	to the bam files '''
	bamfiles, IDs  = [], []
	with open(bams_infofile) as fi:
		for line in fi:
			line_strip = line.rstrip()
			if line_strip:
				try:
					ID, filepath = line_strip.split('\t') #ID	path/to/file
					bamfile = pysam.AlignmentFile(filepath, "rb" )
					IDs.append(ID)
					bamfiles.append(bamfile)
				except Exception:
					sys.exit('ERROR in parsing the bams info file, expected two columns\n')
				#end try
			#end if
		#end for
	#end with

	return bamfiles, IDs
#end def

#################################################################
#	get_ADs (original)
#################################################################
def get_ADs(samfile,chrom,position_actual,REF,MQ_thresh,BQ_thresh):
	''' '''
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	try:
		if chrom==0:
			CC="M"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup(CC, position, position+1)
	#end try

	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip :
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper():
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
							#end if
						else:
							if pileupread.alignment.is_reverse:
								ADr[1]=ADr[1]+1
							else:
								ADf[1]=ADf[1]+1
							#end if
						#end if
					#end if
				#end if
			#end for
		#end if
	#end for

	return ADf,ADr
#end def

#################################################################
#	get_ADs_deletion (original)
#################################################################
def get_ADs_deletion(samfile,chrom,position_actual,REF_0,MQ_thresh,BQ_thresh):
	''' '''
	REF=REF_0
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	try:
		if chrom==0:
			CC="M"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup(CC, position, position+1)
	#end try

	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						indel_val=pileupread.indel
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val>=0:
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
							#end if
						else:
							if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val<0:
								if pileupread.alignment.is_reverse:
									ADr[1]=ADr[1]+1
								else:
									ADf[1]=ADf[1]+1
								#end if
							#end if
						#end if
					#end if
				#end if
			#end for
		#end if
	#end for

	return ADf,ADr
#end def

#################################################################
#	get_ADs_insertion (original)
#################################################################
def get_ADs_insertion(samfile,chrom,position_actual,REF_0,MQ_thresh,BQ_thresh):
	''' '''
	REF=REF_0
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	try:
		if chrom==0:
			CC="M"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		#end if
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		#end if
		if chrom==23:
			CC="X"
		#end if
		if chrom==24:
			CC="Y"
		#end if
		SP=samfile.pileup(CC, position, position+1)
	#end try

	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						indel_val=pileupread.indel
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val<=0:
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
							#end if
						else:
							if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val>0:
								if pileupread.alignment.is_reverse:
									ADr[1]=ADr[1]+1
								else:
									ADf[1]=ADf[1]+1
								#end if
							#end if
						#end if
					#end if
				#end if
			#end for
		#end if
	#end for

	return ADf,ADr
#end def

#################################################################
#	get_ADs_combined (original)
#################################################################
def get_ADs_combined(samfile,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh):
	''' '''
	split_ALT=ALT.split(",")
	if len(split_ALT)>1:
		ADf,ADr = get_ADs(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	#end if
	if len(REF)>1:
		ADf,ADr = get_ADs_deletion(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	#end if
	if len(ALT)>1:
		ADf,ADr = get_ADs_insertion(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	#end if
	ADf,ADr = get_ADs(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)

	return ADf,ADr
#end def

#################################################################
#	get_all_ADs_combined (original)
#################################################################
def get_all_ADs_combined(unrelated_samfiles,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh):
	''' '''
	ADfs=[]
	ADrs=[]
	for samfile in unrelated_samfiles:
		ADf,ADr = get_ADs_combined(samfile,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh)
		ADfs.append(ADf)
		ADrs.append(ADr)
	#end for

	ADfs=np.array(ADfs)
	ADrs=np.array(ADrs)

	return ADfs,ADrs
#end def

#################################################################
#	M1_L_calc_aux (original)
#################################################################
def M1_L_calc_aux(rho,k):
	''' '''
	ALT_count=1
	M1_L_k=np.zeros(ALT_count+1)
	default=np.log((1.-rho)/ALT_count)
	M1_L_k=M1_L_k+default
	M1_L_k[k]=np.log(rho)

	return M1_L_k
#end def

#################################################################
#	M1_L_calc (original)
#################################################################
def M1_L_calc(AD,rho):
	''' '''
	ALT_count=1
	if AD.size-1 != ALT_count:
		sys.exit("ERROR in M1_L_calc\n")
	#end if
	M1_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L_calc_aux(rho,k)
		M1_L.append(M1_L_k)
	#end for

	return M1_L
#end def

#################################################################
#	M2_L_calc_aux (original)
#################################################################
def M2_L_calc_aux(M1_L_k,GT_likelihood_wrt_allele_L):
	''' '''
	ALT_count=1
	if M1_L_k.size-1 != ALT_count:
		sys.exit("ERROR in M2_L_calc_aux\n")
	#end if
	combos=(ALT_count+1)*(ALT_count+2)//2
	temp_table=GT_likelihood_wrt_allele_L+np.tile(M1_L_k.reshape([1,ALT_count+1]),[combos,1])
	M2_L_k=np.zeros(combos)
	for i in range(combos):
		row=temp_table[i,:]
		row_max=np.max(row)
		row=row-row_max
		M2_L_k[i]=np.log(np.sum(np.exp(row)))+row_max
	#end for

	return M2_L_k
#end def

#################################################################
#	M2_L_calc (original)
#################################################################
def M2_L_calc(M1_L,GT_likelihood_wrt_allele_L):
	''' '''
	ALT_count=1
	if (M1_L[0]).size-1 != ALT_count:
		sys.exit("ERROR in M2_L_calc\n")
	#end if
	M2_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L[k]
		M2_L_k=M2_L_calc_aux(M1_L_k,GT_likelihood_wrt_allele_L)
		M2_L.append(M2_L_k)
	#end for

	return M2_L
#end def

#################################################################
#	GT_marg_L_calc (original)
#################################################################
def GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L):
	''' '''
	GT_marg_L=prior_L
	ALT_count=1
	if ADf.size-1 != ALT_count or ADr.size-1 != ALT_count:
		sys.exit("ERROR in GT_marg_L_calc\n")
	#end if
	for k in range(ALT_count+1):
		M2_L_k=M2_L_f[k]
		GT_marg_L=GT_marg_L+ADf[k]*M2_L_k
	#end for
	for k in range(ALT_count+1):
		M2_L_k=M2_L_r[k]
		GT_marg_L=GT_marg_L+ADr[k]*M2_L_k
	#end for

	return GT_marg_L
#end def

#################################################################
#	M3_L_calc_aux (original)
#################################################################
def M3_L_calc_aux(GT_marg_L,M2_L_k):
	''' '''
	M3_L_k=GT_marg_L-M2_L_k

	return M3_L_k
#end def

#################################################################
#	M3_L_calc (original)
#################################################################
def M3_L_calc(GT_marg_L,M2_L):
	''' '''
	ALT_count=1
	if len(M2_L)-1 != ALT_count:
		sys.exit("ERROR in M3_L_calc\n")
	#end if
	M3_L=[]
	for k in range(ALT_count+1):
		M2_L_k=M2_L[k]
		M3_L_k=M3_L_calc_aux(GT_marg_L,M2_L_k)
		M3_L.append(M3_L_k)
	#end for

	return M3_L
#end def

#################################################################
#	M4_L_calc_aux (original)
#################################################################
def M4_L_calc_aux(M3_L_k,GT_likelihood_wrt_allele_L):
	''' '''
	ALT_count=1
	if (GT_likelihood_wrt_allele_L.shape)[1]-1 != 1:
		sys.exit("ERROR in M4_L_calc_aux\n")
	#end if
	combos=(ALT_count+1)*(ALT_count+2)//2
	temp_table=GT_likelihood_wrt_allele_L+np.tile(M3_L_k.reshape([combos,1]),[1,ALT_count+1])
	M4_L_k=np.zeros(ALT_count+1)
	for i in range(ALT_count+1):
		column=temp_table[:,i]
		column_max=np.max(column)
		column=column-column_max
		M4_L_k[i]=np.log(np.sum(np.exp(column)))+column_max
	#end for

	return M4_L_k
#end def

#################################################################
#	M4_L_calc (original)
#################################################################
def M4_L_calc(M3_L,GT_likelihood_wrt_allele_L):
	''' '''
	ALT_count=1
	if (GT_likelihood_wrt_allele_L.shape)[1]-1 != ALT_count:
		sys.exit("ERROR in M4_L_calc\n")
	#end if
	M4_L=[]
	for k in range(ALT_count+1):
		M3_L_k=M3_L[k]
		M4_L_k=M4_L_calc_aux(M3_L_k,GT_likelihood_wrt_allele_L)
		M4_L.append(M4_L_k)
	#end for

	return M4_L
#end def

#################################################################
#	A_marg_L_calc (original)
#################################################################
def A_marg_L_calc(M1_L,M4_L):
	''' '''
	ALT_count=1
	if len(M1_L)-1 != ALT_count:
		sys.exit("ERROR in A_marg_L_calc\n")
	#end if
	A_marg_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L[k]
		M4_L_k=M4_L[k]
		A_marg_L_k=M1_L_k+M4_L_k
		A_marg_L.append(A_marg_L_k)
	#end for

	return A_marg_L
#end def

#################################################################
#	T_term_calc_for_rho (original)
#################################################################
def T_term_calc_for_rho(A_marg_L,AD):
	''' '''
	if len(A_marg_L)!=AD.size:
		sys.exit("ERROR in T_term_calc\n")
	#end if
	ALT_count=AD.size-1
	T1_term=0.
	T2_term=0.
	for k in range(ALT_count+1):
		A_marg_L_k=A_marg_L[k]
		A_marg_temp=np.exp(A_marg_L_k-np.max(A_marg_L_k))
		A_marg=A_marg_temp/np.sum(A_marg_temp)
		T1_term=T1_term+A_marg[k]*AD[k]
		T2_term=T2_term+(1.-A_marg[k])*AD[k]
	#end for

	return T1_term,T2_term
#end def

#################################################################
#	GT_marg_L_to_GT_marg (original)
#################################################################
def GT_marg_L_to_GT_marg(GT_marg_L):
	''' '''
	M=np.max(GT_marg_L)
	GT_marg_L=GT_marg_L-M
	GT_marg=np.exp(GT_marg_L)
	S=np.sum(GT_marg)
	GT_marg=GT_marg/S
	joint_probty_term=np.log(S)+M

	return GT_marg,joint_probty_term
#end def

#################################################################
#	EM_step (original)
#################################################################
def EM_step(ADf_list,ADr_list,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D_original,allele_freq):
	''' '''
	D=np.zeros(3)
	D[0]=D_original[0]
	D[1]=D_original[1]
	D[2]=D_original[2]

	if allele_freq<=0.:
		AF=0.
	else:
		AF=allele_freq
	#end if

	f0=(1.-AF)**2.
	f2=AF**2.
	f1=1.-f0-f2
	D=np.array([f0,f1,f2])*1000.+2.

	T1_f=a-1.
	T2_f=b-1.
	T1_r=a-1.
	T2_r=b-1.
	T_for_prior=D-1.
	joint_probty=             (a-1.)*np.log(rho_f_old)+(b-1.)*np.log(1.-rho_f_old)
	joint_probty=joint_probty+(a-1.)*np.log(rho_r_old)+(b-1.)*np.log(1.-rho_r_old)
	for i in range(3):
		joint_probty=joint_probty+(D[i]-1)*prior_L_old[i]
	#end for

	if len(ADf_list)!=len(ADr_list):
		sys.exit("ERROR1 in EM_step\n")
	#end if

	for i in range(len(ADf_list)):
		ADf=ADf_list[i]
		ADr=ADr_list[i]
		M1_L_f = M1_L_calc(ADf,rho_f_old)
		M1_L_r = M1_L_calc(ADr,rho_r_old)
		M2_L_f = M2_L_calc(M1_L_f,GT_likelihood_wrt_allele_L)
		M2_L_r = M2_L_calc(M1_L_r,GT_likelihood_wrt_allele_L)
		GT_marg_L = GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L_old)
		M3_L_f = M3_L_calc(GT_marg_L,M2_L_f)
		M3_L_r = M3_L_calc(GT_marg_L,M2_L_r)
		M4_L_f = M4_L_calc(M3_L_f,GT_likelihood_wrt_allele_L)
		M4_L_r = M4_L_calc(M3_L_r,GT_likelihood_wrt_allele_L)
		A_marg_L_f = A_marg_L_calc(M1_L_f,M4_L_f)
		A_marg_L_r = A_marg_L_calc(M1_L_r,M4_L_r)

		T1_term_f,T2_term_f = T_term_calc_for_rho(A_marg_L_f,ADf)
		T1_term_r,T2_term_r = T_term_calc_for_rho(A_marg_L_r,ADr)

		T1_f = T1_f + T1_term_f
		T2_f = T2_f + T2_term_f
		T1_r = T1_r + T1_term_r
		T2_r = T2_r + T2_term_r

		GT_marg,joint_probty_term = GT_marg_L_to_GT_marg(GT_marg_L)
		joint_probty = joint_probty + joint_probty_term

		T_for_prior = T_for_prior + GT_marg
	#end for

	rho_f_new=1./(1.+T2_f/T1_f)
	rho_r_new=1./(1.+T2_r/T1_r)
	prior_new=T_for_prior/np.sum(T_for_prior)
	prior_L_new=np.log(prior_new)

	return rho_f_new,rho_r_new,prior_L_new,joint_probty
#end def

#################################################################
#	EM_full (original)
#################################################################
def EM_full(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq):
	''' '''
	joint_probty_s=[]
	joint_probty_new=np.nan
	for i in range(3):
		joint_probty_old=joint_probty_new
		rho_f_new,rho_r_new,prior_L_new,joint_probty_new = EM_step(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)
		rho_f_old=rho_f_new
		rho_r_old=rho_r_new
		prior_L_old=prior_L_new
		joint_probty_s.append(joint_probty_new)
	#end for
	while np.abs(joint_probty_old-joint_probty_new)>10**-7:
		joint_probty_old=joint_probty_new
		rho_f_new,rho_r_new,prior_L_new,joint_probty_new = EM_step(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)
		rho_f_old=rho_f_new
		rho_r_old=rho_r_new
		prior_L_old=prior_L_new
		joint_probty_s.append(joint_probty_new)
	#end while

	return rho_f_new,rho_r_new,prior_L_new,joint_probty_s
#end def

#################################################################
#	GTL_L_calc (original)
#################################################################
def GTL_L_calc(ADf,ADr,rho_f,rho_r,GT_likelihood_wrt_allele_L):
	''' '''
	M1_L_f = M1_L_calc(ADf,rho_f)
	M1_L_r = M1_L_calc(ADr,rho_r)
	M2_L_f = M2_L_calc(M1_L_f,GT_likelihood_wrt_allele_L)
	M2_L_r = M2_L_calc(M1_L_r,GT_likelihood_wrt_allele_L)
	prior_L = np.zeros(3)
	GTL_L = GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L)
	GTL_L=GTL_L-np.max(GTL_L)

	return GTL_L
#end def

#################################################################
#	posterior_probty_calc_exact (original)
#################################################################
def posterior_probty_calc_exact(prior_L,table_L,C_GL_L,M_GL_L,D_GL_L):
	''' '''
	combos=3
	work_column=np.empty(combos**2)
	for I1 in range(combos):
		for I2 in range(combos):
			II=I1*combos+I2
			work_column[II]=prior_L[I1]+prior_L[I2]+M_GL_L[I1]+D_GL_L[I2]
		#end for
	#end for

	work_table=table_L+np.tile(C_GL_L,[combos**2,1])+np.tile(np.reshape(work_column,[combos**2,1]),[1,combos])
	work_table=work_table-np.max(work_table)
	work_table=np.exp(work_table)
	work_table=work_table/np.sum(work_table)
	PP=np.max(np.array([work_table[0][1],work_table[0][2]]))

	return PP,work_table
#end def

#################################################################
#	denovo_P_calc (original)
#################################################################
def denovo_P_calc(ADfs,ADrs,rho_f,rho_r,GT_likelihood_wrt_allele_L,table_L,prior_L):
	''' '''
	M_GL_L=GTL_L_calc(ADfs[0],ADrs[0],rho_f,rho_r,GT_likelihood_wrt_allele_L)
	D_GL_L=GTL_L_calc(ADfs[1],ADrs[1],rho_f,rho_r,GT_likelihood_wrt_allele_L)
	C_GL_L=GTL_L_calc(ADfs[2],ADrs[2],rho_f,rho_r,GT_likelihood_wrt_allele_L) # child is the last one
	PP,work_table = posterior_probty_calc_exact(prior_L,table_L,C_GL_L,M_GL_L,D_GL_L)

	return PP,work_table
#end def

#################################################################
#	PP_calc (original)
#################################################################
def PP_calc(trio_samfiles,unrelated_samfiles,chrom,pos,REF,ALT,allele_freq,MQ_thresh,BQ_thresh):
	''' '''
	ADfs,ADrs = get_all_ADs_combined(unrelated_samfiles,chrom,pos,REF,ALT,MQ_thresh,BQ_thresh)
	ADfs_U=ADfs
	ADrs_U=ADrs
	rho_f_old=0.8
	rho_r_old=0.8
	prior_old=np.array([1./3,1./3,1./3])
	prior_old=prior_old/np.sum(prior_old)
	prior_L_old=np.log(prior_old)
	GT_likelihood_wrt_allele = GT_likelihood_wrt_allele_calc(1)
	GT_likelihood_wrt_allele_L=np.log(GT_likelihood_wrt_allele)
	a=2.
	b=2.
	D=np.array([2.,2,2])

	rho_f_new,rho_r_new,prior_L_new,joint_probty_s = EM_full(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)

	AF_unrel=0.
	for i in range(ADfs.shape[0]):
		temp1=GTL_L_calc(ADfs[i],ADrs[i],rho_f_new,rho_r_new,GT_likelihood_wrt_allele_L)
		temp=temp1+prior_L_new
		temp=temp-np.max(temp)
		temp=np.exp(temp)
		temp=temp/np.sum(temp)

		AF_unrel=AF_unrel+temp[1]+temp[2]*2.
	#end for

	AF_unrel=AF_unrel/2./ADfs.shape[0]

	ADfs,ADrs = get_all_ADs_combined(trio_samfiles,chrom,pos,REF,ALT,MQ_thresh,BQ_thresh)

	table=table_gen(1,1e-8)
	table_L=np.log(table)
	PP,work_table = denovo_P_calc(ADfs,ADrs,rho_f_new,rho_r_new,GT_likelihood_wrt_allele_L,table_L,prior_L_new)

	return PP,ADfs,ADrs,ADfs_U,ADrs_U,rho_f_new,rho_r_new,prior_L_new,AF_unrel
#end def

#################################################################
#	ALT_read_check_in_parents
#################################################################
def ALT_read_check_in_parents(ADfs, ADrs, thr=3):
	''' check if parents have alternate reads over threshold '''
	if len(ADfs) != 3 or len(ADrs) != 3:
		sys.exit("ERROR in retrieving stranded AD counts, missing information for trio\n")
	#end if
	alt_count = ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1]

	if alt_count > thr:
		return False
	else:
		return True
	#end if

#end def

#################################################################
# RUNNER (main function)
#################################################################
def runner(args):
	''' parse the input vcf file and calls the functions to run the analysis '''

	# Variables
	is_allele_freq_thr = True if args['allelefreqthr'] else False
	allele_freq_thr = float(args['allelefreqthr']) if is_allele_freq_thr else 1.
	MQ_thr, BQ_thr = -100., -100.
	replace_RSTR, replace_novoCaller = False, False

	# Buffers
	output_writer = open(args['outputfile'], 'w')
	input_reader = open(args['inputfile'], 'r')

	# Data structures
	variants_passed = []
	trio_column_idxes = {}

	# Opening bam files and getting bam associated IDs
	sys.stderr.write('Reading unrelated and trio bamfiles...\n')
	sys.stderr.flush()

	unrelated_bamfiles, IDs_unrelated = buffering_bams(args['unrelatedbams'])
	trio_bamfiles, IDs_trio = buffering_bams(args['triobams']) # [parent, parent, child]

	# Checking bam info files for trio is complete
	if len(trio_bamfiles) != 3:
		sys.exit('ERROR in bams info file for trio, missing information for some family member\n')
	#end if

	# Reading input file
	count, add_RSTR, added_RSTR = 0, False, False
	add_novoCaller, added_novoCaller = False, False
	line_strip = input_reader.readline().rstrip()
	while line_strip:
		if line_strip.startswith('#'): # reading a header line
			if line_strip.startswith('##'): # header definition line
				if line_strip.split('=')[0] == '##FORMAT':
					add_RSTR = True
				elif line_strip.split('=')[0] == '##INFO':
					add_novoCaller = True
				#end if

				# Adding new format and novoCaller tag definition to header
				if add_RSTR and not added_RSTR:
					added_RSTR = True
					RSTR_tag = '##FORMAT=<ID=RSTR,Number=4,Type=Integer,Description="\
								Reference and alternate allele read counts by strand (Rf,Af,Rr,Ar)">\n'
					output_writer.write(RSTR_tag)
				elif add_novoCaller and not added_novoCaller:
					added_novoCaller = True
					novoCaller_tag = '##INFO=<ID=novoCaller,Number=2,Type=Float,Description="\
									Statistics from novoCaller 2. Format:\'Post_prob|AF_unrel\'">\n'
					output_writer.write(novoCaller_tag)
				#end if

				# Checking if fields are already in the vcf and need to be replaced
				#	if not write header line
				if line_strip.startswith("##FORMAT=<ID=RSTR"):
					replace_RSTR = True
				elif line_strip.startswith("##INFO=<ID=novoCaller")):
					replace_novoCaller = True
				else:
					output_writer.write(line_strip + '\n') # writing header line
				#end if
			elif line_strip.startswith('#CHROM'): # header columns line
				output_writer.write(line_strip) # writing header line

				# Saving the idxes of columns for trio IDs in vcf
				for i, ID in enumerate(line_strip.split('\t')[9:]):
					if ID in IDs_trio:
						trio_column_idxes.setdefault(ID, i)
					#end if
				#end for

				# Checking genotypes trio are complete
				if len(trio_column_idxes) != 3:
					sys.exit('ERROR in vcf structure, missing genotype information for trio\n')
				#end if

				# Writing unrelated samples IDs in the order associate bam files are given
				for ID in IDs_unrelated:
					output_writer.write('\t' + ID)
				#end for
				output_writer.write('\n')
			#end if
		else:
			count += 1
			sys.stderr.write('Analyzing variant: ' + str(count) + '\n')
			sys.stderr.flush()

			# Splitting line and parsing basic information
			try:
				CHROM, POS, _, REF, ALT, _, _, INFO = line_strip.split('\t')[:8]
				chrom_int = encode_chrom(CHROM)
				pos_int = int(POS)
			except Exception:
				sys.exit('ERROR in vcf structure, missing essential columns\n')
			#end try

			# Getting allele frequency
			is_novoAF, allele_freq = False, 0.
			for ID in INFO.split(";"):
				if ID.startswith('novoAF='):
					is_novoAF = True
					try:
						allele_freq = float(ID.split('=')[1])
					except Exception: # novo_AF field is not a float as expected
						if is_allele_freq_thr: # error if allele frequency threshold was provided
							sys.exit('ERROR in input parsing, novoAF field is in the wrong format \
									but looks like you expect to filter based on allele frequency')
						else: # set to 0. if no filtering by allele frequency is requested
							allele_freq = 0.
						#end if
					#end try
					break
				#end if
			#end for

			# Check if novo_AF field missing, error if allele frequency threshold was provided
			if not is_novoAF and is_allele_freq_thr:
				sys.exit('ERROR in input parsing, novoAF field missing but looks like you expect \
						to filter based on allele frequency')
			#end if

			# Calculate statistics
			if allele_freq <= allele_freq_thr: # hard filter on gnomAD allele frequency to substitute the one on unrelated AF
				PP, ADfs, ADrs, ADfs_U, ADrs_U, _, _, _, AF_unrel = PP_calc(trio_bamfiles, unrelated_bamfiles, chrom_int, pos_int, REF, ALT, allele_freq, MQ_thr, BQ_thr)
				#if AF_unrel < 0.01 and ALT_read_checker_in_parents(ADfs, ADrs):
				if ALT_read_check_in_parents(ADfs, ADrs):
					variant = [PP, line_strip, ADfs, ADrs, ADfs_U, ADrs_U, AF_unrel]
					variants_passed.append(variant)
				#end if
			#end if
		#end if

		line_strip=input_reader.readline().rstrip()
	#end while

	sys.stderr.write('\n...Writing results for ' + str(count) + ' variants\n')
	sys.stderr.flush()

	# Writing results
	for variant in sorted(variants_passed, key=lambda x: x[0], reverse=True):
		PP, line_strip, ADfs, ADrs, ADfs_U, ADrs_U, AF_unrel = variant

		# Formatting the printed output
		output_writer.write('\t'.join(line_strip.split('\t')[:7]) + '\t')

		# INFO
		infos = line_strip.split('\t')[7]
		if replace_novoCaller:
			# Removing novoCaller field
			infos_to_write = re.sub(r'novoCaller=.*?;', '', infos)
		else:
			infos_to_write = infos
		#end if

		# Adding novoCaller to INFO
		if infos_to_write[-1] == ';':
			output_writer.write(infos_to_write + 'novoCaller={0}|{1};\t'.format(PP, AF_unrel))
		else:
			output_writer.write(infos_to_write + ';novoCaller={0}|{1};\t'.format(PP, AF_unrel))
		#end if

		# FORMAT
		formats = line_strip.split('\t')[8]
		if replace_RSTR:
			idx_RSTR, formats_to_write = 0, ''
			# Removing RSTR field
			for i, tag in enumerate(formats.split(':')):
				if tag == 'RSTR':
					idx_RSTR = i
				else:
					formats_to_write += tag + ':'
				#end if
			#end for
			genotypes = []
			for gtp in line_strip.split('\t')[9:]:
				gtp_as_list = gtp.split(':')
				del gtp_as_list[idx_RSTR]
				genotypes.append(':'.join(gtp_as_list))
			#end for
		else:
			formats_to_write = formats + ':'
			genotypes = line_strip.split('\t')[9:]
		#end if

		# Adding tag RSTR
		output_writer.write(formats_to_write + 'RSTR\t')

		# Genotypes
		empty_genotype = './.' + ':' * genotypes[0].count(':')

		# Modifying genotypes with AD by strand info
		for i in range(len(genotypes)):
			genotypes[i] += ':'
		#end for
		for i in range(3):
			genotypes[trio_column_idxes[IDs_trio[i]]] += '{0},{1},{2},{3}'.format(int(ADfs[i][0]), int(ADfs[i][1]), int(ADrs[i][0]), int(ADrs[i][1]))
		#end for
		for i in range(len(ADfs_U)):
			genotypes.append(empty_genotype + ':{0},{1},{2},{3}'.format(int(ADfs_U[i][0]), int(ADfs_U[i][1]), int(ADrs_U[i][0]), int(ADrs_U[i][1])))
		#end for

		output_writer.write('\t'.join(genotypes) + '\n')
	# end for

	# Closing files buffers
	input_reader.close()
	output_writer.close()
	for buffer in unrelated_bamfiles:
		buffer.close()
	#end for
	for buffer in trio_bamfiles:
		buffer.close()
	#end for
#end def


#################################################################
#
# MAIN
#
#################################################################
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Bayesian de novo variant caller')

	parser.add_argument('-i', '--inputfile', help='input vcf file, must contain novoAF=<float> in INFO field to filter by allele frequency', required=True)
	parser.add_argument('-o', '--outputfile', help='output file to write results, vcf format', required=True)
	parser.add_argument('-u', '--unrelatedbams', help='tsv file containing ID<TAB>Path/to/file for unrelated bam files \
													used to train the model', required=True)
	parser.add_argument('-t', '--triobams', help='tsv file containing ID<TAB>Path/to/file for family bam files, \
												the proband must be listed as last', required=True)
	parser.add_argument('-a', '--allelefreqthr', help='threshold to filter by population allele frequency', required=False)

	args = vars(parser.parse_args())

	runner(args)

#end if
