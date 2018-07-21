import sys
GUIDE_LENGTH = 20
def revComp(seq):
  dic = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  return ''.join([dic[ele] for ele in list(reversed(list(seq)))])

def findGuide(chromfile, guidefile):
  with open(guidefile, "a+") as outfile:
    with open(chromfile, "r") as infile:
      seq = ''.join(''.join(infile.readlines()).split('\n'))
      for idx, base in enumerate(seq):
        if seq[idx+GUIDE_LENGTH+1:idx+GUIDE_LENGTH+3] == "GG":
          writeString = seq[idx:idx+GUIDE_LENGTH]
          if all(writeString[i].isupper() for i in range(GUIDE_LENGTH)):
            outfile.write(writeString + "\n")
        if seq[idx:idx+2] == "CC":
          writeString = seq[idx+3:idx+3+GUIDE_LENGTH]
          if all(writeString[i].isupper() for i in range(GUIDE_LENGTH)):
            outfile.write(revComp(writeString) + "\n")

def main():
  chromfile = sys.argv[1]
  chromnumb = sys.argv[2]
  guidefile = "../GUIDES/guideschr" + str(chromnumb) + ".txt"
  findGuide(chromfile, guidefile)

if __name__ == "__main__":
  main()
