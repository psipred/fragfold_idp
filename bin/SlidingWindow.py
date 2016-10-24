from optparse import OptionParser
import numpy as np
import random
import itertools
from itertools import combinations
import pickle
import glob
import re


def ReadFile(readf, splitsgn = None):
# wdid (what does it do?): reads a files and returns it as a split list of lines
# input: filename, *optional: split sign for lines (eg. ',' for csv)
# output: list of lines split by ' ' or split sign (if specified)

	finput = open(readf, 'r')
	lines = finput.readlines()
	finput.close()

	if splitsgn:
		plik = []
		for line in lines:
			plik.append( line.replace(splitsgn,' ').split() )
	else:
		plik = lines

	return plik

def TranslateAA(seq):
# wdid: translates three-letter amino-acid code into single-letter code
# input: list or string with 3-letter sequence
# output: string (no spaces) with single-letter amino-acid code

	if type(seq) is list:
		seq = str(' '.join(seq))
	elif type(seq) is str and ',' in seq:
		seq = seq.replace(',', ' ')

	d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
		'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
		'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
		'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	seq = seq.upper().replace(' ','')

	if len(seq) % 3 != 0:
		raise ValueError('Input length should be a multiple of three')

	y = ''
	for i in range(len(seq)/3):
		y += d[seq[3*i:3*i+3]]
	return y

def Centroid(X):
	return sum(X)/len(X)

def GenerateCB(CA,N,C):
# wdid: generates xyz for CB
# input: 3 xyz vectors (lists) - CA, N, C
# output: vector (list) with CB xyz coordinates

# constants
	cacb_dist = 1.538
	sin_teth = 0.7912
	cos_teth = 0.6115

	nca = [CA[0]-N[0], CA[1]-N[1], CA[2]-N[2]]
	cca = [CA[0]-C[0], CA[1]-C[1], CA[2]-C[2]]
	xx = [nca[0]+cca[0], nca[1]+cca[1], nca[2]+cca[2]]
	yy = [nca[1]*cca[2]-nca[2]*cca[1], nca[2]*cca[0]-nca[0]*cca[2], nca[0]*cca[1]-nca[1]*cca[0]]
	sx = cacb_dist * cos_teth / (xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2])**(0.5)
	sy = cacb_dist * sin_teth / (yy[0]*yy[0]+yy[1]*yy[1]+yy[2]*yy[2])**(0.5)
	bx = [ sx*x for x in xx ]
	by = [ sy*y for y in yy ]

	return [ CA[0]+bx[0]+by[0], CA[1]+bx[1]+by[1], CA[2]+bx[2]+by[2] ]

def ReadPDB(filename, chain_mode = 'CA'):
# wdid: reads a PDB file (either main-chain or CA-only) and returns it
# input: filename incl. path, *optional: specify if main-chain [MC] is to be read, main-chain+CB [CB] or [CA] only [default]
# output: numpy array/s of residues and conformations - the order of output N, CA, CB, C, O; order of info in array - [conformation,amino-acid,xyz]

# pre-processing

	if chain_mode != 'MC' and chain_mode != 'CA' and chain_mode != 'CB':
		raise ValueError('Chain mode in ReadPDB() incorrectly specified! Specify either CA, CB or MC')

# read file
	f = open(filename, 'r')
	lines = f.readlines()
	f.close()

#setup variables
	conf = 0
	aa = 0
	max_aa = 0
	at_count = 0

# pre-allocate numpy array/s
# max. 500 residues in *MAX_N_CONF* conformations
	max_n_conf = 500
	CA_xyz = np.zeros((max_n_conf, 500, 3))
	if chain_mode == 'MC' or chain_mode == 'CB':
		C_xyz = np.zeros((max_n_conf, 500, 3))
		N_xyz = np.zeros((max_n_conf, 500, 3))
		O_xyz = np.zeros((max_n_conf, 500, 3))
	if chain_mode == 'CB':
		CB_xyz = np.zeros((max_n_conf, 500, 3))


# read PDB into numpy array
	for line in lines:
		if line.startswith('ATOM') and 'CA' in line:
			CA_xyz[conf,aa,:] = re.sub('-', ' -', line[29:54]).split()
			at_count += 1
		elif line.startswith('ATOM') and ' C ' in line[0:20] and ('MC' in chain_mode or 'CB' in chain_mode):
			C_xyz[conf,aa,:] = re.sub('-', ' -', line[29:54]).split()
			at_count += 1
		elif line.startswith('ATOM') and ' N ' in line[0:20] and ('MC' in chain_mode or 'CB' in chain_mode):
			N_xyz[conf,aa,:] = re.sub('-', ' -', line[29:54]).split()
			at_count += 1
		elif line.startswith('ATOM') and ' O ' in line[0:20] and ('MC' in chain_mode or 'CB' in chain_mode):
			O_xyz[conf,aa,:] = re.sub('-', ' -', line[29:54]).split()
			at_count += 1
		elif line.startswith('ATOM') and 'CB' in chain_mode and 'CB' in line:
			CB_xyz[conf,aa,:] = re.sub('-', ' -', line[29:54]).split()
			at_count += 1
		if at_count == 4 and 'GLY' in line and 'CB' in chain_mode:
			CB_xyz[conf,aa,:] = GenerateCB(CA_xyz[conf,aa,:],N_xyz[conf,aa,:],C_xyz[conf,aa,:])
			at_count += 1

		if at_count == 1 and 'CA' in chain_mode:
			aa += 1
			at_count = 0
		elif at_count == 4 and 'MC' in chain_mode:
			aa += 1
			at_count = 0
		elif at_count == 5 and 'CB' in chain_mode:
			aa += 1
			at_count = 0

		if line.startswith('END\n') or line.startswith('TER') or line.startswith('END '):
			if aa > 0:
				conf += 1
			if aa > max_aa:
				max_aa = aa
			aa = 0

# produce output
	CA_xyz = CA_xyz[:conf,:max_aa,:]
	if chain_mode == 'MC' or chain_mode == 'CB':
		C_xyz = C_xyz[:conf,:max_aa,:]
		N_xyz = N_xyz[:conf,:max_aa,:]
		O_xyz = O_xyz[:conf,:max_aa,:]
	if chain_mode == 'CB':
		CB_xyz = CB_xyz[:conf,:max_aa,:]

# return
	if conf == 1:
		CA_xyz = CA_xyz[0]
		if chain_mode == 'CB' or chain_mode == 'MC':
			C_xyz = C_xyz[0]
			N_xyz = N_xyz[0]
			O_xyz = O_xyz[0]
		if chain_mode == 'CB': CB_xyz = CB_xyz[0]
	if chain_mode == 'CA':
		return CA_xyz
	elif chain_mode == 'MC':
		return N_xyz, CA_xyz, C_xyz, O_xyz
	elif chain_mode == 'CB':
		return N_xyz, CA_xyz, CB_xyz, C_xyz, O_xyz


def SavePDB(xyz_CA, residues, target):

	import datetime.date
	formatPDB = "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n"

	temp_f = 10.0
	x_ca = []
	y_ca = []
	z_ca = []
#	x_cb = []
#	y_cb = []
#	z_cb = []

	for coords_CA in xyz_CA:	# split xyzCA!!!
		x_ca.append(coords_CA.split()[0])
		y_ca.append(coords_CA.split()[1])
		z_ca.append(coords_CA.split()[2])
#		x_cb.append(coords_CB[0])
#		y_cb.append(coords_CB[1])
#		z_cb.append(coords_CB[2])

	res_name = []
	for res in residues:
		res_name.append(res.replace('\n',''))
#		res_name.append(res.replace('\n',''))

	f = open(target, 'w')

	at_no = 1
	res_no = 0
	at_name = ' CA ' # !!!!!!!!

	if header_info:
		f.write('HEADER    %s') % header_info

	f.write('EXPDTA    GENERATED USING PYTHON BY TOMASZ KOSCIOLEK ON %s') % datetime.date.today().strftime('%d/%m/%Y')

	for line in zip(res_name,x_ca,y_ca,z_ca):

		res_no += 1

		f.write(str(formatPDB) % (at_no, at_name, line[0], res_no, float(line[1]), float(line[2]), float(line[3]), temp_f))
		at_no += 1

	f.write('END')
	f.close()



def RMSD_dumb(arrA, arrB, per_res = False):

### find number of atoms
	if len(arrA.shape) > 1:
		n_at = arrA.shape[-2]
		n_dim = arrA.shape[-1]
	else:
		n_at = 1
		n_dim = arrA.shape[0]

	rmsd = []
	for A, B in zip(arrA, arrB):
		rmsd.append([ (A[i]-B[i])**2.0 for i in range(n_dim) ])

	if per_res:
		return [ np.sqrt(sum(r)) for r in rmsd ]	##### CHECK HERE #####
	else:
		return np.sqrt(sum([sum(r) for r in zip(*rmsd)])/n_at)

def Kabsch(arrA, arrB, output = False, per_res = False):

	C = np.dot(np.transpose(arrA), arrB)	# covariance matrix
	v, s, w = np.linalg.svd(C)

	# if the two proteins are reflections of each other
	d = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
	if d:
		s[-1] = -s[-1]
		v[:,-1] = -v[:,-1]

	# rotation matrix
	U = np.dot(v, w)

	# rotate arrA
	arrA = np.dot(arrA, U)

	if output:
		return arrA, RMSD_dumb(arrA, arrB)
	else:
		if per_res:
			return RMSD_dumb(arrA,arrB, per_res)
		else:
			return RMSD_dumb(arrA,arrB)

def RMSD_fit(P, Q, per_res = False, thr_val = 1e-3, max_iter = 20, mean_struc = False):
    """ Varies the distance between P and Q, and optimizes rotation for each step
    until a minimum is found.

    ADAPTED FROM: https://github.com/charnley/rmsd
    """
### pre-processing
    assert(P.shape == Q.shape)
### ##############

    no_steps = 0
    step_size = P.max(0)
    threshold = step_size*thr_val
    rmsd_best = Kabsch(P, Q)
    P_prev = P					# added instead of just *P*
    while True:
	P = P_prev				#
	no_steps += 1
	for i in range(3):
            temp = np.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = Kabsch(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P_prev[:,i] += step_size[i]	#
            else:
                rmsd_new = Kabsch(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P_prev[:,i] -= step_size[i]	#
                else:
                    step_size[i] /= 2
        if (step_size<threshold).all() or no_steps > max_iter:
            break

    if mean_struc:
        return RMSD_dumb(P, Q, True), np.array((P+Q)/2.0)
    else:
     	if per_res:
		return Kabsch(P, Q, per_res=True)
     	else:
		return rmsd_best

def multi_RMSD(*args, **opts):
	''' ALL PAIR-WISE COMBINATIONS '''
# wdid: calculates RMSD (total or per-residue) on a set of models
# input: a set of xyz arrays of coordinates (eg. prot1_xyz, prot2_xyz...) OR a multidimensional array with various conformations (eg. prot[conf,xyz]). NUMPY arrays required
# output: a vector with per-residue RMSD values OR mean total RMSD

#### pre-processing
	per_res = False
	max_iter = False
	thr_val = False
	if opts:
		if opts['per_res']: per_res = True
		if opts['max_iter']: max_iter = opts['max_iter']
		if opts['thr_val']: thr_val = opts['thr_val']

# check input
	if len(args) > 1:
		n_conf = len(args)
		prot1 = args[0]
		if len(prot1.shape) > 1:
			n_aa = prot1.shape[-2]
		else:
			n_aa = 1
		proteins = np.zeros((n_conf,n_aa,3))

		for n in range(n_conf):
			prot = args[n] - Centroid(args[n])
			proteins[n,:,:] = prot
		del prot1
		del prot
	elif len(args) == 1:
		prot1 = args[0]
		n_conf = prot1.shape[0]
		n_aa = prot1.shape[1]
		proteins = np.zeros((n_conf,n_aa,3))
		# do a centroid flip for each protein
		for n in range(n_conf):
			prot = args[0][n] - Centroid(args[0][n])
			proteins[n,:,:] = prot
		del prot1
		del prot
	else:
		raise ValueError('No valid input provided in PerResRMSD')
#############

# calculate the list of all pairs
	combinations_list = list(combinations(range(n_conf), 2))
	n_combinations = len(combinations_list)


# calculate RMSD (tot OR per-res) for all pairs

	pairwise_RMSD = [None] * n_combinations
	comb_no = 0
	for comb in combinations_list:
		pairwise_RMSD[comb_no] = RMSD_fit(proteins[comb[0],:,:], proteins[comb[1],:,:], per_res, thr_val, max_iter)
		comb_no += 1

	pairwise_RMSD = np.array(pairwise_RMSD)**2

# calculate mean value accross all pairs
	mean_RMSD = np.zeros(pairwise_RMSD.shape[-1])
	for pair in pairwise_RMSD:
		mean_RMSD += pair

	if per_res:
		return np.sqrt(mean_RMSD / n_combinations)
	else:
		return np.sqrt(mean_RMSD[0] / n_combinations)


def Sliding_Window(*args, **opts):

#### pre-processing
# prepare proteins np.array, and find out the length of the protein and number of conformations to analyse
	for arg in args:
		if type(arg) is not (np.ndarray or numpy.ndarray):
			raise TypeError('Numpy array required in PerResRMSD function!')
		if 'first_shape' in locals():
			if first_shape != arg.shape:
				raise ValueError('Shape mismatch in PerResRMSD')
		first_shape = arg.shape
	sw = 1
	if opts:
		if opts['sw']:
			sw = opts['sw']

# check input
	if len(args) > 1:
		n_conf = len(args)
		prot1 = args[0]
		if len(prot1.shape) > 1:
			n_aa = prot1.shape[-2]
		else:
			n_aa = 1
		proteins = np.zeros((n_conf,n_aa,3))

		for n in range(n_conf):
			proteins[n,:,:] = args[n]
		del prot1
	elif len(args) == 1:
		proteins = args[0]
		n_conf = proteins.shape[0]
		n_aa = proteins.shape[1]
	else:
		raise ValueError('No valid input provided in PerResRMSD')
##############

# main thing
	SW_RMSD = {key: [] for key in range(n_aa)}
	tmp = np.zeros(sw)

	for aa in range(n_aa): #-sw+1):
		begin_sw = int(aa-((sw-1)/2.0))
		if begin_sw < 0: begin_sw = 0
		end_sw = int(aa+((sw-1)/2.0))
		if end_sw > n_aa: end_sw = n_aa

 		tmp = multi_RMSD(proteins[:,begin_sw:end_sw,:], per_res = True, max_iter=3, thr_val = 1.0)
#		tmp = multi_RMSD(proteins[:,aa:(aa+sw),:], per_res = True, max_iter=3, thr_val = 1.0)
#		print begin_sw, end_sw, tmp

		for t in range(len(tmp)):
#			SW_RMSD[aa+t].append(tmp[t])
			SW_RMSD[begin_sw+t].append(tmp[t])

	out = []
	for key in SW_RMSD:
		out.append( sum(SW_RMSD[key])/len(SW_RMSD[key]) )

	return np.array(out)


def DumpVars(v_name, n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e):
# name

# NMR_n_ens
# NMR_PerRes_RMSD
# NMR_tot_RMSD
# NMR_InterRes_RMSD

# FF_n_ens
# FF_ens_list
# FF_PerRes_RMSD
# FF_tot_RMSD
# FF_InterRes_RMSD
# FF_tot_E
	with open(v_name, 'w') as f:
		pickle.dump([n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e], f)

def LoadVars(v_name):

	with open(v_name) as f:
		n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e = pickle.load(f)
	return n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e


def GenerateEnsStats(nmr_inp, ff_inp, FF_ens_size, FF_n_ens, SW_size, name_pickle):

	##########################
	##	input files	##
	##    1) nmr_inp	##
	##    2) ff_inp		##
	##########################

	#############
	# NMR input #
	#############

	## READ PDB FILES
	NMR_PDB_ens = ReadPDB(nmr_inp) # mode CA

	## FIND NUMBER OF CONFORMATIONS IN THE FILE
	NMR_n_conf = NMR_PDB_ens.shape[0]
	NMR_n_ens = 1
	prot_len = NMR_PDB_ens.shape[1]

	## GENERATE *PER-RESIDUE SLIDING WINDOW RMSD*
	## STORE ARRAY (VARIABLE)
	NMR_PerRes_RMSD = Sliding_Window(NMR_PDB_ens, sw = int(SW_size))

	## INTER-RESIDUE RMSD FULL ALIGNMENT - **NO** SLIDING WINDOW!
	NMR_InterRes_RMSD = np.zeros((NMR_n_ens, prot_len))
	for ens in range(NMR_n_ens):
		NMR_InterRes_RMSD[ens,:] = multi_RMSD(NMR_PDB_ens, per_res = True, max_iter = 6, thr_val = 0.2)

	## CALCULATE RMSD FOR EACH PAIR OF MODELS IN AN ENSEMBLE
	NMR_tot_RMSD = np.zeros((NMR_n_ens, NMR_n_conf*(NMR_n_conf-1)/2))
	for ens in range(NMR_n_ens):
		model_iter = 0
		for confA in range(NMR_n_conf):
		 for confB in range(NMR_n_conf):
		  if confA > confB:
		   NMR_tot_RMSD[ens,model_iter] = multi_RMSD(NMR_PDB_ens[confA,:,:], NMR_PDB_ens[confB,:,:], per_res = False, max_iter = 3, thr_val = 1.0)
		   model_iter += 1

	#############
	# FF  input #
	#############

	## READ *ENS* FILES
	if '*' in ff_inp:
		FF_ens = []
		for inp_file in glob.glob(ff_inp):
			FF_ens.append(ReadPDB(inp_file))
		FF_ens = np.array(FF_ens)
	else:
		FF_ens = ReadPDB(ff_inp) # in CA mode

	## GENERATE 100 RANDOM ENSEMBLES OF 10 STRUCTURES
	n_conf = FF_ens.shape[0]
	prot_len = FF_ens.shape[1]

	conf_list = [i for i in range(n_conf)]
	FF_ens_list = [random.sample(conf_list, FF_ens_size) for i in xrange(FF_n_ens)]

	## CALCULATE *PER-RESIDUE SLIDING-WINDOW RMSD* FOR EACH ENSEMBLE IN FF_ens_list
	FF_PerRes_RMSD = np.zeros((FF_n_ens, prot_len))			# np.array
	for ens in enumerate(FF_ens_list):
		FF_PerRes_RMSD[ens[0],:] = Sliding_Window(FF_ens[ens[1]], sw = int(SW_size))

	## INTER-RESIDUE RMSD FULL ALIGNMENT - **NO** SLIDING WINDOW!
	FF_InterRes_RMSD = np.zeros((FF_n_ens, prot_len))			# np.array
	for ens in enumerate(FF_ens_list):
		FF_InterRes_RMSD[ens[0],:] = multi_RMSD(FF_ens[ens[1]], per_res = True, max_iter = 6, thr_val = 0.2)


	## TOTAL RMSD FOR EACH PAIR OF MODELS IN AN ENSEMBLE
	FF_tot_RMSD = np.zeros((FF_n_ens, FF_ens_size*(FF_ens_size-1)/2))	# np.array
	for ens in range(FF_n_ens):
		model_iterator = 0
		for FF_conf_A, FF_conf_B in itertools.product(FF_ens_list[ens],FF_ens_list[ens]):
			if FF_conf_A != FF_conf_B and (FF_conf_A > FF_conf_B):
				FF_tot_RMSD[ens,model_iterator] = multi_RMSD(FF_ens[FF_conf_A], FF_ens[FF_conf_B], per_res=False, max_iter=3, thr_val = 1.0)
				model_iterator += 1

	## FF ENERGY FOR EACH MODEL IN AN ENSEMBLE
	FF_tot_E = np.zeros((FF_n_ens, FF_ens_size))				# np.array
	f = open(ff_inp, 'r')
	lines = f.readlines()
	f.close()
	FF_E_list = np.zeros(n_conf)
	model_iterator = 0
	for line in lines:
		if line.startswith('HEADER') and 'Potential' in line:
			FF_E_list[model_iterator] = float(line.split()[3])
			model_iterator += 1
	for ens in range(FF_n_ens):
		model_iterator = 0
		for model_id in FF_ens_list[ens]:
			FF_tot_E[ens, model_iterator] = FF_E_list[model_id]
			model_iterator += 1

	if name_pickle:
		# n_n, n_tot, n_pr, n_ir, ff_n, ff_list, ff_tot, ff_pr, ff_ir, ff_e
		DumpVars(name_pickle, NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_ens_list, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E)

	return CalculateEnsStats(NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E)

def CalculateEnsStats(NMR_n_ens, NMR_tot_RMSD, NMR_PerRes_RMSD, NMR_InterRes_RMSD, FF_n_ens, FF_tot_RMSD, FF_PerRes_RMSD, FF_InterRes_RMSD, FF_tot_E):
	#################
	#### RESULTS ####
	#################

	### allocate numpy arrays
	# Rs
	res_Rs = np.zeros(FF_n_ens)

	# FF RMSD
	res_FF_uRMSD = np.zeros(FF_n_ens)
	res_FF_sdRMSD = np.zeros(FF_n_ens)
	res_FF_uIrRMSD = np.zeros(FF_n_ens)
	res_FF_medRMSD = np.zeros(FF_n_ens)
	res_FF_deltaRMSD = np.zeros(FF_n_ens)

	# FF E
	res_FF_uE = np.zeros(FF_n_ens)
	res_FF_sdE = np.zeros(FF_n_ens)
	res_FF_deltaE = np.zeros(FF_n_ens)
	res_FF_medE = np.zeros(FF_n_ens)

	# NMR RMSD
	res_NMR_uRMSD = np.zeros(NMR_n_ens)
	res_NMR_sdRMSD = np.zeros(NMR_n_ens)
	res_NMR_uIrRMSD = np.zeros(NMR_n_ens)
	res_NMR_medRMSD = np.zeros(NMR_n_ens)
	res_NMR_deltaRMSD = np.zeros(NMR_n_ens)


	### Rs
	for ens in range(FF_n_ens):
		res_Rs[ens] = Rs(NMR_PerRes_RMSD, FF_PerRes_RMSD[ens])

#### RMSD ####

########
## FF ##
########
	### mean RMSD
	for ens in range(FF_n_ens):
		res_FF_uRMSD[ens] = np.mean(FF_tot_RMSD[ens])
	### RMSD SD
		res_FF_sdRMSD[ens] = np.std(FF_tot_RMSD[ens]) 	#### I HAVE TO CONSIDER ONLY ACTUAL *UNIQUE* VALUES
	### RMSD spread (deltaRMSD)
		res_FF_deltaRMSD[ens] = max(FF_tot_RMSD[ens])-min(FF_tot_RMSD[ens])
	### median RMSD
		res_FF_medRMSD[ens] = np.median(FF_tot_RMSD[ens])
	### mean inter-residue RMSD
		res_FF_uIrRMSD[ens] = np.mean(FF_InterRes_RMSD[ens])

#########
## NMR ##
#########
	### mean RMSD
	for ens in range(NMR_n_ens):
		res_NMR_uRMSD[ens] = np.mean(NMR_tot_RMSD[ens])
	### RMSD SD
		res_NMR_sdRMSD[ens] = np.std(NMR_tot_RMSD[ens])
	### mean inter-residue RMSD
		res_NMR_uIrRMSD[ens] = np.mean(NMR_InterRes_RMSD[ens])
	### RMSD spread
		res_NMR_deltaRMSD[ens] = max(NMR_tot_RMSD[ens])-min(NMR_tot_RMSD[ens])
	### median RMSD
		res_NMR_medRMSD[ens] = np.median(NMR_tot_RMSD[ens])

##############

#### Energy ####
	### mean E
	for ens in range(FF_n_ens):
		res_FF_uE[ens] = np.mean(FF_tot_E[ens])
	### E SD
		res_FF_sdE[ens] = np.std(FF_tot_E[ens])
	### energy spread (deltaE)
		res_FF_deltaE[ens] = max(FF_tot_E[ens])-min(FF_tot_E[ens])
	### median E
		res_FF_medE[ens] = np.median(FF_tot_E[ens])
################


### RETURN ###
	return res_Rs, res_FF_uRMSD, res_FF_sdRMSD, res_FF_uIrRMSD, res_FF_medRMSD, res_FF_deltaRMSD, res_FF_uE, res_FF_sdE, res_FF_deltaE, res_FF_medE, res_NMR_uRMSD, res_NMR_sdRMSD, res_NMR_uIrRMSD, res_NMR_medRMSD, res_NMR_deltaRMSD
##############
