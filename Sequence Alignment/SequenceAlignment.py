import pandas as pd
import io
import os
import random

class SeqAlign:
    
    def __init__(self, seq1=None, seq2=None, submatrix=None):
        
        while not seq1:
            text = input("Input 1st sequence file name: ")
            seq1 = self.load_seq(text)
        
        while not seq2:
            seq2 = self.load_seq(input("Input 2nd sequence file name: "))
                
        self.seq1 = " " + seq1.strip() # offset for matrix
        self.seq2 = " " + seq2.strip() # not offsetting 2nd
        
        self.seq_type = None # input("\nAre the sequences Nucleotide or Peptide: ")
        
        self.sub_map = submatrix
        # self.alignment = input("\nEnter Alignment Type (Global, Local, Semi-Global): ")
        # self.gap = input("\nEnter gap penalty: ")
        self.get_subscore = lambda i,j :  self.sub_map[(i,j)] if (i,j) in self.sub_map.keys() else self.sub_map[(j,i)]
        
    def load_seq(self, filename):   
        if not os.path.isfile(filename):
               print("File does not exist")
               return None      
        with open(filename, "r") as fin:
            lines = fin.readlines()
            if lines[0].startswith(">"):
                seq = "".join(lines[1:])
            else:
                seq = "".join(lines)
        return seq
    
    def load_matrix(self, filename):
        if not os.path.isfile(filename):
               print("File does not exist")
               return None  
        
        is_pam = input("Is it a non comma-delimited PAM matrix [y/n]: ")
        
        sub_map = {}
        
        if is_pam.lower() == "n":
            with open(filename, "r") as fin:
                lines = fin.readlines()
                if lines[0].startswith(">"):
                    text = "".join(lines[1:])
                else:
                    text = "".join(lines)
                df = pd.read_csv(io.StringIO(text.strip()), sep=',')
                for col in df.columns.values:
                    for j, row in df.iterrows():
                        if (col, df.columns.values[j]) not in sub_map.keys() and (df.columns.values[j], col) not in sub_map.keys():
                            sub_map[(col, df.columns.values[j])] = row[col]
                    
        else:
            with open(filename, "r") as fin:
                text = fin.read()
                text = text.strip("\n").strip("\t")
                lines = text.splitlines()
                labels = list("".join(lines[0].split(" ")))
                table = []
                for i in range (1, len(lines)):
                    table.append(lines[i].split(" ")[1:])
                
                for i in range(0, len(labels)):
                    for j in range (0, len(table)):
                        if (labels[i], labels[j]) not in sub_map.keys() and (labels[j], labels[i]) not in sub_map.keys():
                            sub_map[(labels[i], labels[j])] = float(table[i][j])
                
        return sub_map 
                
    
    def glbl(self, gap:int):
           
        seq1 = self.seq1
        seq2 = self.seq2
        
        sequence_matrix = []
        
        # initialize matrix
        # +1 column and +1 row than seq sizes
        for y in range(0, len(seq2)):
            if y == 0:
                sequence_matrix.append([(0+gap*x, []) for x in range(0,len(seq1))])
            else:
                sequence_matrix.append([(0+gap*y, [])]+[None]*(len(seq1)-1)) 
        # pprint.pprint(sequence_matrix)
        val1 = 0
        val2 = 0
        val3 = 0
        # self.matrix_print(sequence_matrix)
        # d : 0, h : 1, v : 2
        # i : seq2 : row ; j : seq1 : col
        for i in range(1, len(seq2)):
            for j in range(1, len(seq1)):
                # sub or match
                # diag = 0
                diag = sequence_matrix[i-1][j-1][0]
                val1 = (self.get_subscore(seq2[i], seq1[j]) + diag, ["d"])
                
                # deletion
                #hrz = 0
                hrz = sequence_matrix[i][j-1][0]
                val2 = (hrz+gap, ["h"])
                    
                # insertion
                #vrt = 0
                vrt = sequence_matrix[i-1][j][0]
                val3 = (vrt+gap, ["v"]) 

                if val1[0] == val2[0] == val3[0]:
                    sequence_matrix[i][j] = (val1[0], ["d","h","v"])
                elif val1[0] == val2[0]:
                    temp = (val1[0], ["d", "h"])
                    sequence_matrix[i][j] = max([temp, val3], key=lambda x:x[0])
                elif val3[0] == val2[0]:
                    temp = (val3[0], ["v", "h"])
                    sequence_matrix[i][j] = max([temp, val1], key=lambda x:x[0])
                elif val1[0] == val3[0]:
                    temp = (val3[0], ["d", "v"])
                    sequence_matrix[i][j] = max([temp, val2], key=lambda x:x[0])
                else:
                    sequence_matrix[i][j] = max([val1, val2, val3], key=lambda x:x[0])

        print("Global OPT Matrix")
        self.print_matrix(sequence_matrix)
        
        #optimal_alignment = []        
        
        i = len(seq2)-1
        j = len(seq1)-1
        
        return self.backtrack_dfs(sequence_matrix, i, j, gap) # optimal_alignment, self.output_alignment(optimal_alignment)
          
    def local(self, gap):
        
        seq1 = self.seq1
        seq2 = self.seq2
        
        sequence_matrix = []
        
        # initialize matrix
        for y in range(0, len(seq2)):
            if y == 0:
                sequence_matrix.append([(0, [])]*len(seq1))
            else:
                sequence_matrix.append([(0, [])]+[None]*(len(seq1)-1))
                
        val1 = 0
        val2 = 0
        val3 = 0
        max_cell = [0, [None,None]]
        for i in range(1, len(seq2)):
            for j in range(1, len(seq1)):
                # sub or match
                # diag = 0
                diag = sequence_matrix[i-1][j-1][0]
                val1 = (self.get_subscore(seq2[i], seq1[j]) + diag, ["d"])
                
                # deletion
                #hrz = 0
                hrz = sequence_matrix[i][j-1][0]
                val2 = (hrz+gap, ["h"])
                    
                # insertion
                #vrt = 0
                vrt = sequence_matrix[i-1][j][0]
                val3 = (vrt+gap, ["v"]) 

                val_0 = (0, [])
                if val1[0] == val2[0] == val3[0]:
                    sequence_matrix[i][j] = max([val_0, (val1[0], ["d","h","v"])])
                elif val1[0] == val2[0]:
                    temp = (val1[0], ["d", "h"])
                    sequence_matrix[i][j] = max([val_0, temp, val3], key=lambda x:x[0])
                elif val3[0] == val2[0]:
                    temp = (val3[0], ["v", "h"])
                    sequence_matrix[i][j] = max([val_0, temp, val1], key=lambda x:x[0])
                elif val1[0] == val3[0]:
                    temp = (val3[0], ["d", "v"])
                    sequence_matrix[i][j] = max([val_0, temp, val2], key=lambda x:x[0])
                else:
                    sequence_matrix[i][j] = max([val_0, val1, val2, val3], key=lambda x:x[0])
                
                if sequence_matrix[i][j][0] > max_cell[0]:
                    max_cell[0] = sequence_matrix[i][j][0]
                    max_cell[1][0] = i
                    max_cell[1][1] = j
        
        print("Local OPT Matrix")
        self.print_matrix(sequence_matrix)
        
        return self.backtrack_dfs(sequence_matrix, max_cell[1][0], max_cell[1][1], gap)                

    def semi_glbl(self, gap):
        
        seq1 = self.seq1
        seq2 = self.seq2
        
        
        sequence_matrix = []
        for y in range(0, len(seq2)):
            if y == 0:
                sequence_matrix.append([(0, [])]*len(seq1))
            else:
                sequence_matrix.append([(0, [])]+[None]*(len(seq1)-1))
                
        val1 = 0
        val2 = 0
        val3 = 0
        for i in range(1, len(seq2)):
            for j in range(1, len(seq1)):
                # sub or match
                # diag = 0
                diag = sequence_matrix[i-1][j-1][0]
                val1 = (self.get_subscore(seq2[i], seq1[j]) + diag, ["d"])
                
                # deletion
                #hrz = 0
                hrz = sequence_matrix[i][j-1][0]
                val2 = (hrz+gap, ["h"])
                    
                # insertion
                #vrt = 0
                vrt = sequence_matrix[i-1][j][0]
                val3 = (vrt+gap, ["v"]) 

                if val1[0] == val2[0] == val3[0]:
                    sequence_matrix[i][j] = (val1[0], ["d","h","v"])
                elif val1[0] == val2[0]:
                    temp = (val1[0], ["d", "h"])
                    sequence_matrix[i][j] = max([temp, val3], key=lambda x:x[0])
                elif val3[0] == val2[0]:
                    temp = (val3[0], ["v", "h"])
                    sequence_matrix[i][j] = max([temp, val1], key=lambda x:x[0])
                elif val1[0] == val3[0]:
                    temp = (val3[0], ["d", "v"])
                    sequence_matrix[i][j] = max([temp, val2], key=lambda x:x[0])
                else:
                    sequence_matrix[i][j] = max([val1, val2, val3], key=lambda x:x[0])

        max_end = max([[sequence_matrix[a][len(seq1)-1][0], [a,len(seq1)-1]] for a in range(1, len(seq2))]+
                        [[sequence_matrix[len(seq2)-1][b][0], [len(seq2)-1,b]] for b in range(1, len(seq1)-1)], key=lambda x:x[0])

        print("Semi-Global OPT Matrix")
        self.print_matrix(sequence_matrix)

        return self.backtrack_dfs(sequence_matrix, max_end[1][0], max_end[1][1], gap)
               
    def backtrack_dfs(self, matrix, i, j, gap):
        in_i = i
        in_j = j
        
        if not in_i or not in_j: return (), "Matrix cannot be backtracked"
        
        paths = []
        previous_fork = []
        previous_fork_direction = None
        cell = matrix[i][j]
        curr_path = 0
        single_brach_tracker = [] 

        while True:
            
            if cell[1]:
                if "d" in cell[1]:
                    direction = "d"
                    i_move = i-1
                    j_move = j-1
                    score = self.get_subscore(self.seq2[i], self.seq1[j])
                    
                elif "h" in cell[1]:
                    direction = "h"
                    i_move = i
                    j_move = j-1
                    score = gap
                    
                elif "v" in cell[1]:
                    direction = "v"
                    i_move = i-1
                    j_move = j
                    score = gap
                    
                if not paths:        
                    paths.append({"score": score, "path": [[(i,j), direction]]})                   
                else:
                    paths[curr_path]["score"] += score
                    paths[curr_path]["path"].append([(i,j), direction])
                
                # update coords
                i = i_move
                j = j_move
                
                if len(cell[1]) > 1:
                    previous_fork = cell
                    previous_fork_direction = direction
                    single_brach_tracker.append(0)
                else:
                    single_brach_tracker.append(1)
                    
                cell = matrix[i_move][j_move]
            else:
                if 0 not in single_brach_tracker:
                    break
                single_brach_tracker = []
                previous_fork[1].remove(previous_fork_direction)
                i = in_i
                j = in_j
                cell = matrix[i][j]
                paths.append({"score": 0, "path": []})
                curr_path+=1
        
        optimal_alignment = max(paths, key=lambda x:x["score"])
        opt_path = optimal_alignment["path"]
        
        # for path in paths:
        #     if path["score"] == 14:
        #         path["path"].reverse()
        #         print(self.output_alignment(path["path"]))
        #         path["path"].reverse()

        
        opt_path.reverse()
        
        print(f'\nScore: {optimal_alignment["score"]}')    
        return optimal_alignment, self.output_alignment(opt_path)
        
    def output_alignment(self, alignment_list):
        s1 = []
        s2 = []

        
        for val in alignment_list:
            # horizonal is a gap on seq2, vertical is a gap on seq1
            if val[1] == "d":
                s2.append(self.seq2[val[0][0]])
                s1.append(self.seq1[val[0][1]])
            elif val[1] == "h":
                s1.append(self.seq1[val[0][1]])
                s2.append("-")
            else:
                s1.append("-")
                s2.append(self.seq2[val[0][0]])
        
        return " ".join(s1) + "\n" + " ".join(s2)
    

    def print_matrix(self, matrix):
        string = '\t'+'\t'.join(list(self.seq1))
        for i in range(0, len(matrix)):
            string += "\n" + self.seq2[i] + '\t' + '\t'.join([str(cell[0]) for cell in matrix[i]])
        
        print(string)
        
if __name__ == "__main__":
    
    while True:
        aligner = SeqAlign()
        
        submatrix = None
        while not submatrix:
            submatrix = aligner.load_matrix(input("\nEnter Subscore Matrix Filename: "))
                    
        aligner.sub_map = submatrix

        alignment_type = None
        while not alignment_type:
            alignment_type = input("\nEnter Alignment Type (Global, Local, Semi-Global): ")
            if alignment_type.lower() not in ["global", "local", "semi-global"]:
                alignment_type = None
                print("Invalid Alignment Type")
                
        gap = None
        while not gap:
            gap = input("\nEnter gap score: ")
            try:
                gap = float(gap)
            except ValueError:
                print("Invalid Gap Score")
                gap = None
        
        opt_str = None
        if alignment_type.lower() == "global":
            opt_list, opt_str = aligner.glbl(gap)
        
        elif alignment_type.lower() == "local":
            opt_list, opt_str = aligner.local(gap)
            
        elif alignment_type.lower() == "semi-global":
            opt_list, opt_str = aligner.semi_glbl(gap)
        
        print("\nOptimal Alignment:")
        print(opt_str)