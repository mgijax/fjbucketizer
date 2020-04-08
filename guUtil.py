import re
def extractID(col9,patt):
    m = re.search(patt, col9, re.I)
    if m:
        id = m.group('id')
    else:
        id = col9
    return id

if __name__ == "__main__":
    import sys
    for line in sys.stdin:
        tokens = line.split("\t")
        if len(tokens) != 9:
            continue
        c9 = tokens[8]
        print(extractID(c9), "\t", c9)

