import numpy as np

def generatePWMs(n,length_min,length_max=None):
	if length_max is None:
		length_max = length_min +1

	matrices = []
	for i in range(n):
		L = np.random.randint(length_min,length_max)

		# noise = 
		matrix= np.random.uniform(-4,2,(L,4))
		matrix = np.exp(matrix)
		matrix/= np.sum(matrix,1)[:,None]
		matrices.append(matrix)
	return matrices

def generateSequences(n,matrices,length_min,length_max = None):
	if length_max is None:
		length_max = length_min +1

	# cwms = []
	# for pwm in matrices:
	# 	cwms = np.cumsum(pwm,1)

	sequences = []
	origins=[]
	for s in range(n):
		L = np.random.randint(length_min,length_max)

		baseSeq = np.random.choice(4,L)
		
		nRandomSeqProb = 0.05
		if np.random.rand() > nRandomSeqProb:
			
			mat_idx = np.random.randint(0,len(matrices))
			origins.append(mat_idx)
			pwm_length,_ = matrices[mat_idx].shape

			start_idx = np.random.randint(0,L-pwm_length)
			for j in range(0,pwm_length):
				base = np.random.choice(4,p=matrices[mat_idx][j])

				baseSeq[start_idx+j] = base
		else:
			origins.append(-1)
		sequences.append("".join(['A','T','C','G'][q] for q in baseSeq))
	return sequences,origins


def generateFakeData(nMatrix,nSequence):


	pwms = generatePWMs(nMatrix,8,12)
	sequences,origins = generateSequences(nSequence,pwms,25,30)

	for i,pwm in enumerate(pwms):
		path = f"pwms/synthetic_{i}.pfm"
		np.savetxt(path,pwm.T,fmt='%.3f')

	# print(s)
	f = open("synthetic_reads.txt","w")
	f.write("\n".join(sequences)+"\n")
	f.close()
	f = open("synthetic_origins.txt","w")
	f.write("".join([f"{i} {origins[i]}\n" for i in range(len(origins))]))
generateFakeData(15,10000)