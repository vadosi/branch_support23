'''
Input is an alignment, output is a set of trees built by fastme with different orders of taxa
'''

import os
from sys import argv, exit
from random import shuffle
from subprocess import Popen, PIPE

# FUNCTIONS

def readmatrix(fn):
  '''readmatrix: input -- name of file with distance matrix in "new PHYLIP format" (supported by FastME)
  returns distance matrix as a dictionary whose keys are frozen sets of pairs of sequence names'''
  fs = frozenset
  matrix = list()
  names = list()
  result = dict()
  inflow = open(fn, "r")
  try:
    num = int(inflow.readline())
  except:
    print ("File", fn, "has wrong format")
    return None
  raw = 0
  while raw < num:
    line = inflow.readline()
    s = line.strip().split()
    if len(s) - 1 != num:
      print ("File", fn, "has wrong format")
      return None
    names += [s[0]]
    matrix += [list()]
    column = 0
    for x in s[1:]:
      try:
        matrix[raw] += [float(x)]
      except:
        print ("File", fn, "has wrong format")
        return None
      column += 1
    raw += 1
  # while raw < num
  inflow.close()
  for i in range(len(names)):
    for j in range(len(names)):
      key = fs((names[i], names[j]))
      if key not in result.keys():
        result[key] = matrix[i][j]
      else:
        if result[key] != matrix[i][j]:
          print ("Matrix is not symmetric:", matrix[i][j], "is not equal to", matrix[j][i])
          return None
  return result

def listnames(matrix):
  '''listnames: input -- distance matrix, returns list of names'''
  resset = set()
  for pair in matrix.keys():
    resset |= pair
  return list(resset)

def writematrix(names, matrix):
  '''writenames: input -- list of names and distance matrix, returns string -- matrix in PHYLIP format,
  with order according to the list'''

  result = format(len(names), "5d") + "\n"
  for x in names:
    result += x
    for y in names:
      result += format(matrix[frozenset((x, y))], "10.6f")
    result += "\n"
  return result

# MAIN

try:
  infile = argv[1]
except:
  print ("Usage:\npython", argv[0], "alignment n")
  print ("Example:\npython", argv[0], "alignment.phy 100")
  exit(1)

try:
  num = int(argv[2])
except:
  num = 100

if "." in infile:
  matrixname = ".".join(infile.split(".")[:-1]) + ".dist"
else:
  matrixname = infile + ".dist"

command = ["fastme", "-i", infile, "-O", matrixname, "-p", "-c"]
fastme = Popen(command, stdout = PIPE, stderr = PIPE)
(log, err) = fastme.communicate()
if fastme.returncode != 0:
  print ("FastME error")
  print (err)
  exit(1)
os.remove(infile + "_fastme_stat.txt")

matrix = readmatrix(matrixname)
names = listnames(matrix)

result = matrixname + '.nwk'

if os.access(result, os.F_OK):
  os.remove(result)
for _ in range(num):
  shuffle(names)
  tmpflow = open("tmp", "w")
  tmpflow.write(writematrix(names, matrix))
  tmpflow.close()
  command = ["fastme", "-i", "tmp", "-o", "temp.tre", "-m", "B"]
  fastme = Popen(command, stdout = PIPE, stderr = PIPE)
  (log, err) = fastme.communicate()
  os.remove("tmp")
  os.remove("tmp_fastme_stat.txt")
  os.system("cat temp.tre >> " + result)
  os.remove("temp.tre")
print(result, '  DONE')