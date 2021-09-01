import os
import pandas as pd
import io


class tw_pair_dict(dict):
    def __len__(self):
	        return dict.__len__(self) / 2
    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        dict.__setitem__(self, (key[1], key[0]), value)
# (a, b) : value

class Matrix:
    def __init__(self):
        self.seq = ''      
        self.matrix = tw_pair_dict()
        
    # load distance matrix
    def load_matrix(self, filename):
        if not os.path.isfile(filename):
                print("File does not exist")
                return None  
        
        sub_map = tw_pair_dict()
        with open(filename, "r") as fin:
            lines = fin.readlines()
            num_ent = lines[0]
            text = "".join(lines[1:])
            self.seq = lines[1].strip().split("\t")
            df = pd.read_csv(io.StringIO(text.strip()), sep='\t')
            for col in df.columns.values:
                for j, row in df.iterrows():
                    if (col, df.columns.values[j]) not in sub_map.keys() and (df.columns.values[j], col) not in sub_map.keys():
                        sub_map[(col, df.columns.values[j])] = row[col]
        self.matrix = sub_map
        return sub_map

    def print_matrix(self, matrix):
            string = '   '+'\t'.join(list(self.seq))
            for i in range(0, len(self.seq)):
                string += "\n" + self.seq[i] + ': ' + '\t'.join([str(matrix[(j, self.seq[i])]) for j in self.seq])
            
            print(string)


class UPGMA:
    def __init__(self, matrix:Matrix):
        self.codes = matrix.seq         # matrix codes
        self.length = len(self.codes)   # initial length of matrix
        self.matrix = matrix.matrix     # distance matrix
        
    def run(self):
        # og_mat = self.matrix
        clusters = {}   # cluster cluster pairs
        newick = {}
        code_order = {} # og index of codes
        
        for i in range(0, self.length):
            code_order[self.codes[i]] = i
            newick[self.codes[i]] = self.codes[i]
            clusters[self.codes[i]] = {}
            clusters[self.codes[i]][0] = self.codes[i]
        
        max_dist = 999999
        num_clusters = self.length
        
        while num_clusters > 1:
            # calculates distance and finds smallest    
            smallest = max_dist
            smallestI = 0
            smallestJ = 0
            
            cluster_list = list(clusters.keys())
            
            for i in range(0, num_clusters-1):
                for j in range(i+1, num_clusters):
                    print(f"\nCluster {i}: {cluster_list[i]} Cluster {j}: {cluster_list[j]}")

                    temp_dist = 0
                    # calcualate distance betqween clusters s the unweighted average of all distances
                    for k in range(0, len(cluster_list[i])):
                        for m in range(0, len(cluster_list[j])):
                            # print(clusters[cluster_list[i]][k])
                            # print(clusters[cluster_list[j]][m])
                            temp_dist += self.matrix[(clusters[cluster_list[i]][k],clusters[cluster_list[j]][m])]
                    temp_dist = temp_dist/(len(cluster_list[i])*len(cluster_list[j]))
                    
                    if temp_dist < smallest:
                        smallest, smallestI, smallestJ = temp_dist, i, j
            
            clusterI = cluster_list[smallestI]
            clusterJ = cluster_list[smallestJ]
            merge = clusterI+clusterJ
            print(f"Merging Clusters: {clusterI} and {clusterJ} with distance {smallest}")
            i = 0
            clusters[merge] = {}
            for j in range(0, len(clusterI)):
                clusters[merge][i] = clusters[clusterI][j]
                i+=1
            
            for j in range(0, len(clusterJ)):
                clusters[merge][i] = clusters[clusterJ][j]
                i+=1
            
            newick[merge] = f"({newick[clusterI]},{newick[clusterJ]})"
            
            del clusters[clusterI]
            del clusters[clusterJ]
            del newick[clusterI]
            del newick[clusterJ]
            num_clusters-=1
            print(f"Newick format: {newick}")
        
        cluster_list = list(clusters.keys())
        print(f"Newick format: {newick[cluster_list[0]]}")


class NJ:
    def __init__(self, matrix:Matrix):
        self.codes = matrix.seq         # matrix codes
        self.length = len(self.codes)   # initial length of matrix
        self.matrix = matrix.matrix     # distance matrix

    def run(self):
        newick = {}
        clusters = {}
        for i in range(0, self.length):
            newick[self.codes[i]] = self.codes[i]
            clusters[self.codes[i]] = {}
            for j in range(0, self.length):
                clusters[self.codes[i]][self.codes[j]] = float(self.matrix[(self.codes[i], self.codes[j])])
            
            
        max_dist = 999999
        num_clusters = self.length
        while num_clusters > 2:
            cluster_list = list(clusters.keys())
            r_vals = [0]*num_clusters
            for i in range(0, num_clusters):
                temp_val = 0
                for j in range(0, num_clusters):
                    # if(cluster_list[i] not in clusters.keys() or cluster_list[j] not in clusters[cluster_list[i]].keys()):
                    #     continue
                    temp_val += clusters[cluster_list[i]][cluster_list[j]]
                r_vals[i] = temp_val/(num_clusters-2)
        
            print("R Values")
            for i in range(0, num_clusters):
                print(f"{cluster_list[i]} {r_vals[i]} ", end=' ', flush=True)
            # print("\n",end=' ', flush=False)
            
            smallest = max_dist
            smallestI = 0
            smallestJ = 0
            
            print("\nTD Matrix: ")
            for i in range(0, num_clusters-1):
                for j in range(i+1, num_clusters):
                    t_d = clusters[cluster_list[i]][cluster_list[j]] - r_vals[i] - r_vals[j]
                    
                    print(f"TD({cluster_list[i]},{cluster_list[j]}) = {t_d},  ", end=' ',flush=True)
                    if t_d < smallest:
                        smallest, smallestI, smallestJ = t_d, i, j
                print("\n")    
            
            # merge and insert clusters into cluster dictionary
            clusterI = cluster_list[smallestI]
            clusterJ = cluster_list[smallestJ]
            merge = clusterI+clusterJ
            
            # get branch lengths
            branch_1 = (clusters[clusterJ][clusterI] + r_vals[smallestI] - r_vals[smallestJ]) / 2            
            branch_2 = (clusters[clusterJ][clusterI] + r_vals[smallestJ] - r_vals[smallestI]) / 2
        
            print(f"Merging Clusters {clusterI} and {clusterJ}")
            print(f"\tDistance between {clusterI} and ancestral node = {branch_1}")
            print(f"\tDistance between {clusterJ} and ancestral node = {branch_2}")
            
            for i in range(0, num_clusters):
                if cluster_list[i] != clusterI and cluster_list[i] != clusterJ:
                    dist_1 = clusters[cluster_list[i]][clusterI]
                    dist_2 = clusters[cluster_list[i]][clusterJ]
                    # assign new distance to pair
                    if merge not in clusters.keys():
                        clusters[merge] = {}
                        clusters[merge][merge] = 0
                    clusters[merge][cluster_list[i]] = (dist_1 + dist_2 - clusters[clusterI][clusterJ]) / 2
                    clusters[cluster_list[i]][merge] = clusters[merge][cluster_list[i]]
            
            newick[merge] = f"({newick[clusterI]},{newick[clusterJ]})"
            
            # delete stuff no longer neeed
            del clusters[clusterI]
            del clusters[clusterJ]
            del newick[clusterI]
            del newick[clusterJ]
            
            # new list of cluster keys
            cluster_list = list(clusters.keys())
            num_clusters -=1
        
        print(f"Distance between remaining clusters: {clusters[cluster_list[0]][cluster_list[1]]}{clusters[cluster_list[1]][cluster_list[0]]}")

        out = " "
        for j in range(0, num_clusters):
            out += newick[cluster_list[j]] + " "
        print(out)
            
            
if __name__=="__main__":
    matrix = Matrix()
    fn = input("Enter Distance Matrix Name: ")
    if os.path.exists(fn):
        dm = matrix.load_matrix(fn)
        matrix.print_matrix(dm)
    else:
        print("File does not exist: ")
        exit(0)
    
    algos = {
        "1": UPGMA,
        "2": NJ
    }
    
    algo = input("Select Algorithm:\n\t1. UPGMA\n\t2. Neighbor Joining\nEnter 1 or 2:")
    if algo in algos.keys():
        a = algos[algo](matrix)
        a.run()
    else:
        print("Invalid Input")
        exit(0)