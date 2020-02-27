import numpy

nres = 144
out = open('mat_yDNA1_hDNA1.align', 'w')
for i in range(nres):
    line = f"DA yDNA1 {i+1} DA hDNA1 {i+1} 9 \n"
    #print(line)
    out.write(line)
out.close()
