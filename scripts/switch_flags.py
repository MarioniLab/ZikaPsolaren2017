import pysam
import sys

if len(sys.argv)!=4:
    print sys.argv[0], "<IN> <OUT> <first/second>"
    sys.exit(1)

curfile=sys.argv[1]
readsin=pysam.Samfile(curfile, "rb")
readsout=pysam.Samfile(sys.argv[2], "wb", template=readsin)
addfirst=(sys.argv[3]=="first")

# Defining a function to pull out 5' clip length.
def get_5clip(xread):
    if xread.is_unmapped:
        return -1
    cig=xread.cigar
    offset=0
    if xread.is_reverse:
        if cig[-1][0]==5:
            offset=cig[-1][1]
    else:
        if cig[0][0]==5:
	    offset=cig[0][1]
   
    # STAR only hard clips the chimeric segment and soft clips the main read.
    # This hack ensures that the chimeric segment is reported as the 5' end
    # if both the main and chimeric alignments have no 5' hard clipping.
    if xread.is_supplementary: 
        if offset==0:
            offset=-1    
    return offset

# Setting variables.
okay=True
currentname=""
laststore=[]
index_5=0
clip_5=0

# Looping across the file.
while okay:
    try:
        nextread=readsin.next()
    except StopIteration:
        okay=False
        
    if laststore and (not okay or nextread.qname!=currentname):
        # Need to pick the primary alignment as the 5'-most read (i.e., least hard-clipping on the 5' end).    
        for i, oldread in enumerate(laststore):
            oldread.is_paired=True
            oldread.is_supplementary=(i!=index_5) 
            if addfirst:
                oldread.is_read1=True
            else:
                oldread.is_read2=True
            readsout.write(oldread)
        laststore=[]

    if not okay or nextread.is_secondary: # Skipping secondary alignments produced by STAR.
        continue

    new_clip_5=get_5clip(nextread)
    if laststore:
        laststore.append(nextread)
        if new_clip_5 < clip_5: 
            clip_5=new_clip_5
            index_5=len(laststore)-1
    else:             
        laststore=[nextread]
        currentname=nextread.qname       
        index_5=0
        clip_5=new_clip_5

# Closing file handles.
readsin.close()
readsout.close()

