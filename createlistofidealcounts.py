import sys

if len(sys.argv) != 3:
    print "You must supply input and output filenames"
    sys.exit()

inputfilename = sys.argv[1]
outputfilename = sys.argv[2]

with open(inputfilename) as f:
    list0 = []
    for line in f:
        v = line.lstrip()
        [c,h] = v.split(' ')
        list0.append([int(h,16),int(c)]);

    list0.sort(key=lambda item: (item[0], item[1]))

with open(outputfilename,'w') as f:
    for i in range(0,len(list0)):
        f.write(format(list0[i][0],'x') + "," + str(list0[i][1]) + "\n");

