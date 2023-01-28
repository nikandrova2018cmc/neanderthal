from testScript import *
def main():
    start = input("Choose either you want to generate everytime a new ancestry branch or use the current one for tesitng:\n" \
            "1-Generate new branch(Y)\n" \
            "2-Don't generate(N)")
    if start.capitalize() == "Y":
        begin = Neanderthal()
        vector = [[], [], [], []]
        x = []
        while begin.s < 30:
            x.append(begin.s)
            for k in range(4):
                vector[k].append(begin.slow(begin.s)[k])
            begin.s += 5
        print(vector)
        for j in range(3):
            print(vector[j])
            print(x)
            begin.graph(x,vector[j])
        vector[3] = sorted(vector[3],key=lambda x:x[0])
        print(vector[3])
        X = [x[0] for x in vector[3]]
        Y = [x[1] for x in vector[3]]
        print(X)

    elif start.capitalize() == "N":
        begin = Neanderthal()
        vector = [[], [], [], []]
        x = []
        while begin.f < 30:
            x.append(begin.f)
            for k in range(4):
                vector[k].append(begin.fast(begin.f)[k])
            begin.f += 5
        print(vector)
        for j in range(3):
            print(vector[j])
            print(x)
            begin.graph(x,vector[j])
        vector[3] = sorted(vector[3],key=lambda x:x[0])
        print(vector[3])
        X = [x[0] for x in vector[3]]
        Y = [x[1] for x in vector[3]]
        print(X)
    else:
        begin = Neanderthal()
        vector = [[], [], [], []]
        x = []
        while begin.l < 30:
            x.append(begin.l)
            for k in range(4):
                vector[k].append(begin.length_dna(begin.l)[k])
            begin.l += 5
        print(vector)
        for j in range(3):
            print(vector[j])
            print(x)
            begin.graph(x,vector[j])
        vector[3] = sorted(vector[3],key=lambda x:x[0])
        print(vector[3])
        X = [x[0] for x in vector[3]]
        Y = [x[1] for x in vector[3]]
        print(X)
if __name__ == "__main__":
    main()
