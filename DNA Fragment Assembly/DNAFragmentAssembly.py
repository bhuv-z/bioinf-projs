# DNA Fragment Assembly

import os

def load_seqs(filename):
    
    if not os.path.exists(filename):
        print("File does not exist")
        exit(0)
        
    with open(filename) as fin:
        return [str(f) for f in fin.read().strip().splitlines()]


def assemble(frags):
    temp = frags.copy()
    
    for i in temp:
        for j in temp:
            if i!= j and i in j:
                frags.remove(i)
    del temp
    
    pair_score = {}
    
    for i in frags:
        for j in frags:
            if i == j:
                continue
            if i in j:
                frags.remove(i)
                print(f"Removing {i} since the entire fragment is a part of {j}")
                continue
            predecessor = None
            for k in range(0,len(i)):
                if j.startswith(i[k:]):
                    predecessor = i[k:]
                    break
            if predecessor:
                pair_score[(i, j)] = [predecessor, len(predecessor)]
                print(f"Overlap between {i} and {j} is {predecessor} of length {len(predecessor)}")
            else:
                pair_score[(i, j)] = ["", 0]
    
    pair_score_sorted = {y[0]:y[1] for y in sorted(pair_score.items(), key=lambda x:x[1][1], reverse=True)}
    
    print("\nSorted: ")
    for (x,y) in pair_score_sorted.keys():
        print(f"Overlap from {x} to {y} is {pair_score_sorted[(x,y)][1]}")
    
    path = []
    path_info = []
    edge_tracker = []
    while True:
        for pair in pair_score_sorted.keys():
            if path == []:
                path.append(pair[0])
                path_info.append(pair_score_sorted[pair])
                path.append(pair[1])
                continue
            if pair[0] == path[len(path)-1]:            
                possible_pairs = {key[0]:key[1] for key in pair_score_sorted.items() if key[0][0]==pair[0] and key not in edge_tracker}
                for key in possible_pairs.keys():
                    if key[1] not in path:
                        path.append(key[1])
                        path_info.append(possible_pairs[key])
                        break
        if len(path) == len(frags):
            break
    
    print("\nPath Assembly")
    print("-------------")
    assembled = ""
    indent = 0
    out = ""
    for i in range(0, len(frags)):
        out = " "*indent+" ".join(path[i])
        print(out)
        if i < len(frags)-1:
            indent = len(out)
            if path_info[i][1] != 0:
                indent -= path_info[i][1]*2
                assembled += path[i][:0-path_info[i][1]]
            else:
                assembled += path[i]
            indent += 1
            
        else:
            assembled += path[i]
    
    score = str(sum([i[1] for i in path_info]))
    
    print("-"*len(out))
    print("\nHamiltonian Path Score: "+score)
    print("\nAssembled Path: "+" ".join(assembled)+"\n")
            
        
                        
    
if __name__ == "__main__":
    filename = input("Enter Fragment Filename: ")
    frags = load_seqs(filename)
    
    for i in range(0, len(frags)):
        print(f"Fragment{i} is {frags[i]}")    
    
    assemble(frags)
    
    
    
    