import numpy as np

forward_map = ['A','C','G','T']
complement_map = ['T','G','C','A']

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

def fancyEncode(sequence):
	bitfield = 0
	for base in sequence:
		bitfield = (bitfield << 2) + base
	return bitfield
def fancyDecode(value,length):
	out = ""
	for i in range(length-1,-1,-1):
		pos = i * 2
		base = (value >> pos) &	3
		out += forward_map[base]
	return out
	
def rc_encoding(sequence):
	value = fancyEncode(sequence)
	length = len(sequence)

	bitfield = 0
	for i in range(0,length):
		pos = i*2
		base = (value >> pos) & 3
		complement = base ^ 3
		bitfield = (bitfield << 2) + complement
	return bitfield

def validateComplement(sequence,baseSeq,direction):
	if direction > 0:
		#validates that the encoder + decoder work
		enc = fancyEncode(baseSeq)
		dec = fancyDecode(enc,len(baseSeq))
		if (sequence != dec):
			raise ValueError("Sequence encoder-decoder does not function")
	else:
		enc = rc_encoding(baseSeq)
		dec = fancyDecode(enc,len(baseSeq))
		
		if (sequence != dec):
			raise ValueError("Reverse complementing does not function")
		



def generateSequences(n,matrices,length_min,length_max = None):
	if length_max is None:
		length_max = length_min +1

	sequences = []
	origins=[]
	directions = []
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

		if np.random.rand() > 0.5:
			map = forward_map
			seq = baseSeq
			dir = 1
		else:
			map = complement_map
			seq = np.flip(baseSeq)
			dir = -1

		sequence = "".join(map[q] for q in seq)
		validateComplement(sequence,baseSeq,dir)
		sequences.append(sequence)
		directions.append(dir)
	return sequences,origins,directions


def generateFakeData(nMatrix,nSequence):


	pwms = generatePWMs(nMatrix,8,12)
	sequences,origins,directions = generateSequences(nSequence,pwms,25,30)

	for i,pwm in enumerate(pwms):
		path = f"pwms/synthetic_{i}.pfm"
		np.savetxt(path,pwm.T,fmt='%.3f')

	# print(s)
	f = open("synthetic_reads.txt","w")
	f.write("\n".join(sequences)+"\n")
	f.close()
	f = open("synthetic_origins.txt","w")

	dirmap = ["-","","+"]
	f.write("".join([f"{i} {origins[i]} {dirmap[1+directions[i]]}\n" for i in range(len(origins))]))
generateFakeData(15,10000)