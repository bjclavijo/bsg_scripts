import pybsg as bsg

class LDGLineConsensus(object):
    """docstring for LDGLineConsensus."""
    def __init__(self, ws, mldg, ldg, line):
        self.ws = ws
        self.mldg=mldg
        self.ldg=ldg
        self.line=line
        self.initial_line=line

    def add_connecting_nodes(self):
        new_line=[self.line[0]]
        added=0
        completed=0
        total=0
        for p in range(len(self.line)-1):
            total+=1
            node_from=self.line[p]
            node_to=self.line[p+1]
            #dist=[x.dist for x in self.ldg.links[abs(node_from)] if x.source==-node_from and x.dest==node_to][0]
            node_from_fw=[(x[0],x[1]) for x in self.mldg.fw_neighbours_by_distance(node_from,2)]
            for i in range(len(node_from_fw)):
                if node_from_fw[i][1] in self.line:#==node_to:
                    node_from_fw=node_from_fw[:i]
                    break
            node_to_bw=[(x[0],-x[1]) for x in self.mldg.fw_neighbours_by_distance(-node_to,2)]
            for i in range(len(node_to_bw)):
                if node_to_bw[i][1] in self.line:#==node_from:
                    node_to_bw=node_to_bw[:i]
                    break
            #print ("Finding common nodes between %d and %d, %d bp appart"%(node_from,node_to,dist))
            node_from_fw_ids=[x[1] for x in node_from_fw]
            node_to_bw_ids=[x[1] for x in node_to_bw]
            common=[x for x in node_from_fw_ids if x in node_to_bw_ids]
            if common:
                completed+=1
                #print (node_from_fw)
                #print (node_to_bw)
            for x in common:
                #print("common node %d: %d bp" % (x,len(self.ws.sg.nodes[abs(x)].sequence)))
                new_line.append(x)
                added+=1
            new_line.append(node_to)
        self.line=new_line
        #print ("%d nodes added in %d/%d gaps"% (added,completed,total))
        return added

    def produce_consensus(self,join_overlapping=True,walk_with_nanopore=True,fill_with_nanopore=True):
        #put first sequence read_first
        #iterate over gap and sequences
        ## if join overlapping, run join_overlapping heuristic
        ## if no gap_seq and walk_with_nanopore, run nanopore_walking
        ## if no gap_seq and fill_with_nanopore, run nanopore_filling
        ## if gap_seq, add gap seq, if not put Ns (at least 100)
        ## put sequence for next sequence in line
        pass
