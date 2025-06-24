#!/usr/bin/env python
from numpy import *
import torch
import argparse
import subprocess
import re
import gzip
import csv

def fastq_iter(handle):
    while True:
        lines = [handle.readline() for _ in range(4)]
        if not lines[0]:
            break
        yield lines

def open_maybe_gz(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

print("Defining Functions")

alphabet = "acgt"
def decode(x): # convert list of ints to string
    s = "".join([alphabet[xx] for xx in x if alphabet[xx]])
    return s
def encode(st): # convert a string into a list of ints
    x = [alphabet.index(ss) for ss in st if ss in alphabet]
    return x
def seqtomer(seq) : # return list of trimers in a seq
    ans = [int(16*seq[k] + 4*seq[k+1] + seq[k+2]) for k in range(len(seq)-2)]
    return ans
def mertobin(mer) : # trimer list to occupancy uint64 bitmap
    ibin = 0
    for m in mer :
        ibin |= (1 << m)
    return ibin

def makeerrors(seq,srate,irate,drate) : #error mode applied to a sequence seq
    # note: modifies (also returns) seq
    n = len(seq)
    # substitutions
    ns = random.binomial(n,srate*1.3333) # 3/4 of substitutions are "effective"
    ndx = random.randint(low=0,high=n,size=ns)
    vals = random.randint(low=0,high=4,size=ns)
    seq[ndx] = vals
    # deletions
    nd = random.binomial(n,drate)
    ndx = random.randint(low=0,high=n,size=nd)
    seq = delete(seq,ndx)    
    # insertions (at specified rate into smaller seq)
    ni = random.binomial(len(seq),irate)
    ndx = random.randint(low=0,high=len(seq)+1,size=ni)
    vals = random.randint(low=0,high=4,size=ni)
    seq = insert(seq,ndx,vals)
    # pad or truncate to original length
    nn = len(seq)
    if nn > n :
        seq = seq[:n]
    elif nn < n :
        seq = concatenate((seq,random.randint(low=0,high=4,size=n-nn)))
    return seq

m0 = uint64(0x5555555555555555)  # binary: 0101...
m1 = uint64(0x3333333333333333)  # binary: 00110011..
m2 = uint64(0x0f0f0f0f0f0f0f0f)  # binary:  4 zeros,  4 ones ...
m3 = uint64(0x00ff00ff00ff00ff)  # binary:  8 zeros,  8 ones ...
m4 = uint64(0x0000ffff0000ffff)  # binary: 16 zeros, 16 ones ...
m5 = uint64(0x00000000ffffffff)  # binary: 32 zeros, 32 ones
def popcount(x):
# https://github.com/google/jax/blob/6c8fc1b031275c85b02cb819c6caa5afa002fa1d/jax/lax_reference.py#L121-L150
    x = (x & m0) + ((x >>  1) & m0)  # put count of each  2 bits into those  2 bits
    x = (x & m1) + ((x >>  2) & m1)  # put count of each  4 bits into those  4 bits
    x = (x & m2) + ((x >>  4) & m2)  # put count of each  8 bits into those  8 bits
    x = (x & m3) + ((x >>  8) & m3)  # put count of each 16 bits into those 16 bits
    x = (x & m4) + ((x >> 16) & m4)  # put count of each 32 bits into those 32 bits
    x = (x & m5) + ((x >> 32) & m5)  # put count of each 64 bits into those 64 bits
    return x

try :
    from popcll_torch import popcll
    #https://github.com/iamgroot42/popcll_torch , thanks to Anshuman Suri (as9rw@virginia.edu)!
    mypopcount = popcll
except :
    mypopcount = popcount
    print("popcll_torch not found, so using slower popcount")

def find_runs(x):
    #Find runs of consecutive items in an array.
    # credit: https://gist.github.com/alimanfoo/c5977e87111abe8127453b21204c1065
    n = x.shape[0]
    # find run starts
    loc_run_start = empty(n, dtype=bool)
    loc_run_start[0] = True
    not_equal(x[:-1], x[1:], out=loc_run_start[1:])
    run_starts = nonzero(loc_run_start)[0]
    # find run values
    run_values = x[loc_run_start]
    # find run lengths
    run_lengths = diff(append(run_starts, n))
    return run_values, run_starts, run_lengths

def chemfilter(seq, homomax=3, atmax=22, cgmax=22) :
    # returns whether seq satisfies chemistry constraints
    bc = bincount(seq,minlength=4)
    if (bc[0]+bc[3] > atmax) or (bc[1]+bc[2] > cgmax) : return False
    _,_,run_lengths = find_runs(seq)
    if max(run_lengths) > homomax : return False
    return True

def allcoses(mer, tcosvecs) : # correlate a mer against all the cosine templates
    mmer = torch.LongTensor(mer).to(device)
    ncos = tcosvecs.shape[0]
    cosvec = torch.zeros(ncos, 64, dtype=torch.float, device=device)
    for k in range(ncos) :
        source = tcoses[k, torch.arange(len(mmer), dtype=torch.long, device=device)]
        cosvec[k,:].index_add_(0,mmer,source)
    return torch.sum(torch.unsqueeze(cosvec,dim=1)*tcosvecs,dim=2)

def prank(arr, descending=False) : # returns rank of each element in torch array
    argsrt = torch.argsort(arr, descending=descending)
    rank = torch.zeros(arr.shape, dtype=torch.float, device=device)
    rank[argsrt] = torch.arange(len(argsrt),dtype=torch.float,device=device)
    return rank    

class ApproximateLevenshtein :
    def __init__(s, M, N, Q, zsub, zins, zdel, zskew):
        torch.set_grad_enabled(False) # just in case not done elsewhere!
        s.M = M # length of seq1
        s.N = N # length of each seq2
        s.Q = Q # number of seq2s
        (s.zsub, s.zins, s.zdel, s.zskew) = (zsub, zins, zdel, zskew)
        s.tab = torch.zeros(N+1,Q, device=device)
        
    def __call__(s,seq1,seq2) :
        assert (len(seq1) == s.M) and (seq2.shape[1] == s.N) and (seq2.shape[0] == s.Q)
        s.tab[:,:] = (s.zskew * torch.arange(s.N+1., device=device)).unsqueeze(1) # force broadcast
        for i in range(1,s.M+1) :
            diag = s.tab[:-1,:] + torch.where(seq1[i-1] == seq2.t(), 0., s.zsub) # diagonal move
            s.tab[0,:] += s.zskew
            s.tab[1:,:] += s.zdel # down move
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], diag) # or diag if better
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], s.tab[:-1,:] + s.zins) # right move
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], s.tab[:-1,:] + s.zins) # repeat (>= 0 times) as you can afford
           # N.B.: M >= N gives better approx than N > M, so change arg order accordingly
        return s.tab[s.N,:]

class ParallelLevenshtein :
    def __init__(s, M, N, Q, zsub, zins, zdel, zskew):
        torch.set_grad_enabled(False) # just in case not done elsewhere!
        s.M = M # length of seq1
        s.N = N # length of each seq2
        s.Q = Q # number of seq2s
        (s.zsub, s.zins, s.zdel, s.zskew) = (zsub, zins, zdel, zskew)
        MN1 = M + N + 1
        s.blue = torch.zeros(Q, MN1, MN1, device=device)
        s.bluef = s.blue.view(Q, MN1 * MN1)
        s.ndxr = torch.zeros(M*N, dtype=int, device=device) # index of mer matches array into flat blue
        for m in torch.arange(M,device=device) :
            for n in torch.arange(N,device=device) :
                s.ndxr[n + N*m] = (3*M+2*N+2) + (M+N)*m + (M+N+2)*n
        s.lls = torch.zeros(MN1+1,dtype=torch.int,device=device)       
        s.rrs = torch.zeros(MN1+1,dtype=torch.int,device=device)       
        for i in range(2,MN1+1) :
            s.lls[i] = abs(M - i + 1) + 1
            s.rrs[i] = (M+N-1) - abs(- i + 1 + N )

    def __call__(s, seq1, sseq2): # single seq1, tensor of sseq2s
        assert (len(seq1) == s.M) and (sseq2.shape[1] == s.N) and (sseq2.shape[0] == s.Q)    
        (M1,N1,MN,MN1,MN2) = (s.M + 1, s.N + 1, s.M + s.N, s.M + s.N + 1, s.M + s.N + 2)
        abmatch = (seq1.view(1,s.M,1) != sseq2.view(s.Q,1,s.N)).type(torch.float) * s.zsub
        s.bluef[:,s.ndxr] = abmatch.view(s.Q,s.M*s.N)
        s.bluef[:,torch.arange(s.M,MN2*N1,MN2)] = (s.zskew*torch.arange(N1,device=device)).unsqueeze(0)
        s.bluef[:,torch.arange(s.M,MN*M1,MN)] = (s.zskew*torch.arange(M1,device=device)).unsqueeze(0)
        for k in torch.arange(2,MN1,device=device) :
            ll = s.lls[k]
            rr = s.rrs[k]
            slis = torch.arange(ll,rr+1,2,device=device)
            s.blue[:,k,slis] = torch.minimum(
                s.blue[:,k,slis] + s.blue[:,k-2,slis],
                torch.minimum(
                    s.blue[:,k-1,slis-1] + s.zdel,
                    s.blue[:,k-1,slis+1] + s.zins
                )
            )
        return s.blue[:,-1,s.N]

print("all functions now defined.")



print("all functions now defined.")

#Parsing arguments

parser = argparse.ArgumentParser(prog='ST', description="""Decode ST barcodes""")
parser.add_argument ("-N", dest="Num_bar", type=int, default=42000, help="Number of spots/barcodes")
parser.add_argument ("-M", dest="Len_bar", type=int, default=36, help="Length barcodes")
parser.add_argument ("-B", dest="Bar_file", default="./barcodes", help="File with barcodes, generated with the script 1_BarcodeGen.py")
parser.add_argument("-R1", "--input_R1", dest="R1_path", help="Path to the R1 fastq files")
parser.add_argument("-R2", "--input_R2", dest="R2_path", help="Path to the R2 fastq files")
parser.add_argument("--Ntriage", dest="Ntriage", type=int, help="Number of barcodes passed from triage to Levenshtein")
parser.add_argument("--Nthresh", dest="Nthresh", type=int, help="user set to Levenshtein score greater than which is called an erasure")
parser.add_argument("--out-prefix", dest="out_prefix", help="Prefix for out file")
#new arguments for paralellization
parser.add_argument ("-C", dest="Cuda", type=int, default=0, help="Cuda number")
parser.add_argument ("--Seg_id", dest="segmentid", help="For example 1/2, 2/2")
args = parser.parse_args()

print("arguments parsed")

# extracting argument values
N=args.Num_bar
M=args.Len_bar
Ncos=4
filename = args.Bar_file
read1= args.R1_path
read2= args.R2_path
Ntriage = args.Ntriage
Nthresh = args.Nthresh
Outfile = args.out_prefix
cudano = args.Cuda
segmentid = args.segmentid


# Read barcodes from CSV file
codes = []
with open(filename, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        # Extract barcode from 'cell' column and convert to lowercase, add newline for compatibility
        codes.append(row['cell'].lower() + '\n')

assert N == len(codes)
assert M == len(codes[0][:-1])
allseqs = []
alltrimers = []
allbitmaps = zeros(N, dtype=uint64)
cosvecs = torch.zeros((Ncos,N,64),dtype=torch.float)
coses = zeros((Ncos,M)) 
for k in range(Ncos) :
    coses[k,:] = cos(pi*arange(M)*(k+1.)/(M-1.))
tcoses = torch.tensor(coses, dtype=torch.float)
for i,code in enumerate(codes) :
    seq = encode(code[:-1])
    allseqs.append(seq)
    mer = seqtomer(seq)
    mmer = torch.LongTensor(mer)
    alltrimers.append(mer)
    allbitmaps[i] = mertobin(mer)
    for k in range(Ncos):
        source = tcoses[k,arange(M-2)]
        cosvecs[k,i,:].index_add_(0,mmer,source)
print("finished making code auxilliary tables, now pickling")        
pickledict = {"N" : N, "M" : M, "allseqs" : allseqs,
    "alltrimers" : alltrimers, "allbitmaps" : allbitmaps, "coses" : coses, "cosvecs" : cosvecs}


# DECODING Step 1. Move the code (assumed loaded above) and reads (from file) to the GPU.

print("moving to CUDA device")
# verify cuda and set device
if torch.cuda.is_available() :
    device = torch.device(f"cuda:{cudano}")
    cudaname = torch.cuda.get_device_name()
else :
    raise RuntimeError("Required GPU not found! Exiting.")
print(f"{segmentid}: Using {device} device {cudaname}.")

# Use subprocess to execute the command
out = subprocess.Popen(['zcat', read1], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
# Use wc -l to count lines directly without storing the intermediate content
val = subprocess.Popen(['wc', '-l'], stdin=out.stdout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].strip()
#checking errors
if not val.isdigit() : raise RuntimeError(f"failed to count lines in {read1}")
totlines = int(val)
#parse segment ID
num,denom = re.match('([0-9]+)/([0-9]+)',segmentid).groups()
startline = int(int(totlines*(float(num)-1.)/float(denom))/4)
endline = int(int(totlines*float(num)/float(denom))/4)

#Copy tensors to gpu
torch.set_grad_enabled(False)
tallseqs = torch.tensor(array(allseqs), device=device)
talltrimers = torch.tensor(array(alltrimers), device=device)
tallbitmaps = torch.tensor(allbitmaps.astype(int64), dtype=torch.int64, device=device) # 
tcoses = torch.tensor(coses, dtype=torch.float, device=device)
tcosvecs = cosvecs.to(device)

# DECODING Step 2 (main step). Primary and secondary triage, followed by Levenshtein

torch.set_grad_enabled(False)

mydist = ApproximateLevenshtein(M,M,Ntriage, 1.,1.,1.,1.)
#mydist = ParallelLevenshtein(M,M,Ntriage, 1.,1.,1.,1.)

Ncos = tcosvecs.shape[0]
dists = torch.zeros(Ncos+1, N, dtype=torch.float, device=device) # will hold distances for each read
allrank = torch.zeros(Ncos+1 ,N, dtype=torch.float, device=device)

#Define a demux function and then run for each matepair
def bar_demux (seq):
  # primary and secondary triage
  mer = seqtomer(seq)
  foo = int64(uint64(mertobin(mer)))
  seqbin = torch.tensor(foo,dtype=torch.int64,device=device)
  xored = torch.bitwise_xor(seqbin,tallbitmaps)
  dists[0,:] = 64. - mypopcount(xored)
  cosvec = torch.zeros(Ncos, 64, dtype=torch.float, device=device)
  for k in range(Ncos) :
      cosvec[k,mer] =  tcoses[k, torch.arange(len(mer), dtype=torch.long, device=device)]
  dists[1:,:] = torch.sum(torch.unsqueeze(cosvec,dim=1)*tcosvecs,dim=2) # all cosine distances
  for k in range(Ncos+1) :
      allrank[k,:] = prank(dists[k,:], descending=True) # rank them all
  offset = 1.
  fm = torch.prod(offset+allrank,dim=0)
  fmargsort = torch.argsort(fm) #store in a torch array the index of each barcode in the original list of barcodes
  # Levenshtein distance
  tseq1 = torch.tensor(seq,device=device) #take the read
  tseq2 = tallseqs[fmargsort[:Ntriage],:] #take the top Ntriage barcodes
  ans = mydist(tseq1,tseq2) #array with all the distances of each read with all the barcodes that passed the triage
  ia = torch.argmin(ans) # index into fmargsort of best, closest distance
  ibest = fmargsort[ia] # index of best codeword in codes
  ibest = (ibest if ans[ia] <= Nthresh else -1)
  if ibest == -1:
    return ("erasure")
  else:
    return codes[ibest]

print("Starting Demux")

Read1_out = Outfile + "_" + num + "of" + denom + "_R1.fastq.gz"
Read2_out = Outfile + "_" + num + "of" + denom + "_R2.fastq.gz"

with open_maybe_gz(read1) as r1_handle, open_maybe_gz(read2) as r2_handle, \
     gzip.open(Read1_out, "wt") as w1_handle, gzip.open(Read2_out, "wt") as w2_handle:
    i = 0
    for r1_lines, r2_lines in zip(fastq_iter(r1_handle), fastq_iter(r2_handle)):
        # 1) segment slicing
        if i < startline:
            i += 1
            continue
        if i >= endline:
            break

        # 2) normalize and strip Ns
        read = r1_lines[1].strip().lower().replace("n", "")

        # 3) grab exactly M bases (pad on the left with 'a' if too short)
        if len(read) < M:
            raw = read[-M:]
            raw = raw.rjust(M, "t")
        else:
            raw = read[:M]

        # 4) encode & demux
        seq = encode(raw)
        barcode = bar_demux(seq)

        # 5) if we got a barcode, trim, rename, write
        if barcode != "erasure":
            # trim the first M bases off R1
            trimmed_r1_seq = r1_lines[1][M:]
            trimmed_r1_qual = r1_lines[3][M:]
            # rebuild the header: keep everything after the first '_'
            prefix, rest = r2_lines[0].strip().split("_", 1)
            newname = f"{prefix}_{barcode.strip().upper()}_{rest}\n"
            # Write to output
            w1_handle.write(newname)
            w1_handle.write(trimmed_r1_seq)
            w1_handle.write(r1_lines[2])
            w1_handle.write(trimmed_r1_qual)
            w2_handle.write(newname)
            w2_handle.write(r2_lines[1])
            w2_handle.write(r2_lines[2])
            w2_handle.write(r2_lines[3])
        i += 1

print(f"Demux completed. Output written to {Read1_out} and {Read2_out}.")