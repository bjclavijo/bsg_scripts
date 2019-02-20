import pybsg as bsg
from statistics import median
import gzip

class ReadGraphNode(object):
    def __init__(self,id):
        self.id=id
        self.read_pos={}
        self.nexts=[]
        self.prevs=[]

class ReadGraph(object):
    def __init__(self,k):
        self.k=k
        self.sequence_to_nodes={}
        self.nodes=[]
        self.reads=[]
        self.nkcov={}

    def add_read(self,sequence,read_id):
        last_node=-1
        added=reused=0
        self.reads.append(read_id)
        for x in range(len(sequence)-self.k+1):
            seq=sequence[x:x+self.k]
            #print("checking for kmer: '%s'"%seq)
            if seq not in self.sequence_to_nodes.keys():
                #add new node with ++last_id
                added+=1
                self.sequence_to_nodes[seq]=len(self.nodes)
                self.nodes.append(ReadGraphNode(len(self.nodes)))
            else: reused+=1
            this_node=self.sequence_to_nodes[seq]
            self.nodes[this_node].read_pos[read_id]=x
            if last_node!=-1:
                self.nodes[last_node].prevs.append(this_node)
                self.nodes[this_node].nexts.append(last_node)
            last_node=this_node
        #print ("%d new nodes and %d re-uses"%(added,reused))

    def reads_pos_dif_to_median(self):
        for rid in self.reads:
            pmedians=[]
            for n in self.nodes:
                if rid in n.read_pos and len(n.read_pos)>3:
                    #print (n.read_pos)
                    pmedians.append(n.read_pos[rid]-median(n.read_pos.values()))
            print(rid,pmedians)
            hist(pmedians)
            show()

    def place_and_score(self,sequence):
        score=0
        max_score=0
        for x in range(len(sequence)-self.k+1):
            seq=sequence[x:x+self.k]
            max_score+=1
            if seq in self.sequence_to_nodes: score+=1
        #print ("score is %d/%d"%(score,max_score))
        return score/max_score

    def extend_paths(self,sg,origin,dest,max_dist,min_coverage=.85):
        '''Recursive function, returns [[solution1, solution2, ..., solutionN] , [unfinished1, unfinished2, ..., unfinishedN] ]'''
        #print ("Extending paths from %d to %d max_dist %d, min_coverage %.2f"%(origin,dest,max_dist,min_coverage))
        sols=[]
        extpaths=[[origin]]
        while extpaths:
            #print("\nexpansion round starting with extpaths=%s and sols=%s\n"%(extpaths,sols))
            new_extpaths=[]
            for p in extpaths:
                #print("Extending path %s"%p)
                links=sg.get_fw_links(p[-1])
                #print("links from last node %d: %s"%(p[-1],[l.dest for l in links]))
                for l in links:
                    #print ("%s + %d" %(p,l.dest))
                    if l.dest==dest:
                        #print ("  -> SOLVED!!!")
                        s=[x for x in p]
                        s.append(dest)
                        sols.append(s)
                    else:
                        ep=[x for x in p]
                        ep.append(l.dest)
                        epl=sum([len(sg.nodes[abs(n)].sequence) for n in ep[1:]])
                        #print ("  -> dist %d"%epl)
                        if epl>max_dist:
                            #print ("  -> TOO FAR ")
                            continue
                        pck=self.get_path_covered_kmers(ep[1:],sg)
                        pckp=pck[0]/pck[1]
                        #print ("  -> cov %d/%d = %.2f"%(pck[0],pck[1],pckp))
                        if pckp<min_coverage:
                            #print ("  -> TOO LOW ")
                            continue
                        new_extpaths.append(ep)
                    #print (new_extpaths)
            extpaths=new_extpaths
        return sols

    def get_path_covered_kmers(self,path,sg):
        cov=0
        total=0
        for n in path:
            x=self.get_node_covered_kmers(n,sg)
            cov+=x[0]
            total+=x[1]
        return [cov,total]

    def get_full_path_covered_kmers(self,path,sg):
        '''returns count of covered kmers for path. Uses internal cache'''
        p=bsg.SequenceGraphPath(sg)
        for x in path: p.nodes.append(x)
        sequence=p.get_sequence()

        score=0
        max_score=0
        for x in range(len(sequence)-self.k+1):
            seq=sequence[x:x+self.k]
            max_score+=1
            if seq in self.sequence_to_nodes: score+=1
        return [score,max_score]

    def get_node_covered_kmers(self,node_id,sg):
        '''returns count of covered kmers for node. Uses internal cache'''
        sequence=sg.nodes[abs(node_id)].sequence
        if node_id<1:
            n=bsg.Node(sequence)
            n.make_rc()
            sequence=n.sequence
        if node_id not in self.nkcov.keys():
            score=0
            max_score=0
            for x in range(len(sequence)-self.k+1):
                seq=sequence[x:x+self.k]
                max_score+=1
                if seq in self.sequence_to_nodes: score+=1
            self.nkcov[node_id]=[score,max_score]
        return self.nkcov[node_id]

    def walk_on_sg(self,sg,node_from,node_to):
        return self.extend_paths(sg,node_from,node_to, 20000)

class LDGLineConsensus(object):
    """docstring for LDGLineConsensus."""
    def __init__(self, ws, mldg, ldg, line, long_reads_file=""):
        self.ws = ws
        self.mldg=mldg
        self.ldg=ldg
        self.line=line
        self.initial_line=line
        self.long_reads_file=long_reads_file
        self.read_seqs={}
        if self.long_reads_file:
            with gzip.open(self.long_reads_file,'rt') as f:
                for l in f:
                    if len(l)==0: continue
                    if l[0]=='>':
                        name=l[1:]
                    else:
                        self.read_seqs[int(name)]=l.strip()
        #check the reads that are in the line node
        self.node_reads={}
        for x in line: self.node_reads[abs(x)]=set()

        for rid in range(len(ws.long_read_mappers[0].filtered_read_mappings)):
            for m in ws.long_read_mappers[0].filtered_read_mappings[rid]:
                n=abs(m.node)
                if n in self.node_reads.keys():
                    self.node_reads[n].add(rid)

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

    def walk_with_nanopore(self,node_from,node_to):
        print('Walking with nanopore between %d and %d'%(node_from,node_to))
        #get the reads, chop the appropriate regions
        chopped_seqs=[]
        for rid in self.node_reads[abs(node_from)].intersection(self.node_reads[abs(node_to)]):
            rc=False
            rnodes=[m.node for m in self.ws.long_read_mappers[0].filtered_read_mappings[rid]]
            if -node_from in rnodes and -node_to in rnodes: #read is reversed
                sstart=max([m.qEnd for m in self.ws.long_read_mappers[0].filtered_read_mappings[rid] if m.node==-node_from])
                send=min([m.qStart for m in self.ws.long_read_mappers[0].filtered_read_mappings[rid] if m.node==-node_to])
                rc=True
            elif node_from in rnodes and node_to in rnodes:
                sstart=max([m.qEnd for m in self.ws.long_read_mappers[0].filtered_read_mappings[rid] if m.node==node_from])
                send=min([m.qStart for m in self.ws.long_read_mappers[0].filtered_read_mappings[rid] if m.node==node_to])
            else: continue
            if send<sstart: continue
            sstart=max(sstart-199,0)
            send=min(send+199,len(self.read_seqs[rid]))
            if len(self.read_seqs)>1:
                seq=self.read_seqs[rid][sstart:send]
            else:
                raise NotImplementedError("can't yet get sequences from datastore, sorry!")
            if rc:
                tn=bsg.Node(seq)
                tn.make_rc()
                seq=tn.sequence
            chopped_seqs.append(seq)
        print('Chopped %d sequences, with sizes %s'% (len(chopped_seqs),[len(x) for x in chopped_seqs]))


        #Construct a directional DBG with the nanopore data
        rg=ReadGraph(15)
        for i in range(len(chopped_seqs)): rg.add_read(chopped_seqs[i],i)
        print('graph with %d nodes created, %d nodes with coverage>1'%(len(rg.nodes),len([x for x in rg.nodes if len(x.read_pos)>1])))
        #optional: try to delimit all nodes that are not anchors that are in direct connection to these nodes?

        # try to find a path in the graph that is coherent with the PB DDBG
        for min_id in [90,80,70,60,50]:
            rg_paths=rg.walk_on_sg(self.ws.sg,node_from,node_to)
            for rg_path in rg_paths:
                rg_path=rg_path[1:-1]
                ck=rg.get_full_path_covered_kmers(rg_path,self.ws.sg)
                print("Solution path %s, score %d / %d = %.3f" %(rg_path,ck[0],ck[1],ck[0]/ck[1]))
            if rg_paths:break
        if not rg_paths:
            print("Can't find solution with typical walking, FW:")
            for fn in self.ws.sg.get_fw_links(node_from):
                nseq=self.ws.sg.nodes[abs(fn.dest)].sequence
                if fn.dest<0:
                    tn=bsg.Node(nseq)
                    tn.make_rc()
                    nseq=tn.sequence
                print ("next neighbour is node %d (%d bp) with score %.2f "%(fn.dest,len(nseq),rg.place_and_score(nseq)))
                for fn2 in self.ws.sg.get_fw_links(fn.dest):
                    nseq2=self.ws.sg.nodes[abs(fn2.dest)].sequence
                    if fn2.dest<0:
                        tn=bsg.Node(nseq2)
                        tn.make_rc()
                        nseq2=tn.sequence
                    print ("   next next neighbour is node %d (%d bp) with score %.2f "%(fn2.dest,len(nseq2),rg.place_and_score(nseq2)))
                    for fn3 in self.ws.sg.get_fw_links(fn2.dest):
                        nseq3=self.ws.sg.nodes[abs(fn3.dest)].sequence
                        if fn3.dest<0:
                            tn=bsg.Node(nseq3)
                            tn.make_rc()
                            nseq3=tn.sequence
                        print ("      next next next neighbour is node %d (%d bp) with score %.2f "%(fn3.dest,len(nseq3),rg.place_and_score(nseq3)))
            print("Can't find solution with typical walking, BW:")
            for fn in self.ws.sg.get_fw_links(-node_to):
                nseq=self.ws.sg.nodes[abs(fn.dest)].sequence
                if fn.dest<0:
                    tn=bsg.Node(nseq)
                    tn.make_rc()
                    nseq=tn.sequence
                print ("next neighbour is node %d (%d bp) with score %.2f "%(fn.dest,len(nseq),rg.place_and_score(nseq)))
                for fn2 in self.ws.sg.get_fw_links(fn.dest):
                    nseq2=self.ws.sg.nodes[abs(fn2.dest)].sequence
                    if fn2.dest<0:
                        tn=bsg.Node(nseq2)
                        tn.make_rc()
                        nseq2=tn.sequence
                    print ("   next next neighbour is node %d (%d bp) with score %.2f "%(fn2.dest,len(nseq2),rg.place_and_score(nseq2)))
                    for fn3 in self.ws.sg.get_fw_links(fn2.dest):
                        nseq3=self.ws.sg.nodes[abs(fn3.dest)].sequence
                        if fn3.dest<0:
                            tn=bsg.Node(nseq3)
                            tn.make_rc()
                            nseq3=tn.sequence
                        print ("      next next next neighbour is node %d (%d bp) with score %.2f "%(fn3.dest,len(nseq3),rg.place_and_score(nseq3)))

    def produce_consensus(self,join_neighbours=True,join_overlapping=True,walk_with_nanopore=True,fill_with_nanopore=True):
        #put first sequence read_first

        unconnected_distances=[]
        n_stretches=0
        direct_connections=0
        short_indirect=0
        s=""
        #first node, no overlap possible to any previous
        n=self.line[0]
        if n>0:
            s+=self.ws.sg.nodes[n].sequence
        else:
            rnode=bsg.Node(self.ws.sg.nodes[-n].sequence)
            rnode.make_rc()
            s+=rnode.sequence

        #iterate over gap and sequences

        ## put sequence for next sequence in line
        self.joined_fw=[]
        self.dist_fw=[]
        for ni in range(1,len(self.line)):

            n=self.line[ni]
            pn=self.line[ni-1]
            try:
                ns1=[x for x in self.ws.linked_read_mappers[0].tag_neighbours[abs(pn)] if x.node==abs(n)][0].score
            except:
                ns1=0
            try:
                ns2=[x for x in self.ws.linked_read_mappers[0].tag_neighbours[abs(n)] if x.node==abs(pn)][0].score
            except:
                ns2=0
            connected=False
            mlc=[x.dist for x in self.mldg.get_fw_links(pn) if x.dest==n]
            print("\nJoining nodes %d (%d bp) and %d (%d bp), neighbour scores %.3f, %.3f" %
                (pn,len(self.ws.sg.nodes[abs(pn)].sequence),n,len(self.ws.sg.nodes[abs(n)].sequence),ns1,ns2))
            #HEURISTIC #1: direct connection on SG

            if join_neighbours and mlc and median(mlc)<200:
                for fl in self.ws.sg.get_fw_links(pn):
                    if fl.dest==n:
                        #print("nodes %d and %d are connected on sg, distance %d !" % (pn,n,fl.dist))
                        if fl.dist<0:
                            s=s[:fl.dist] #distance is negative, so overlap is removed
                            connected=True
                            direct_connections+=1
                            self.joined_fw.append(True)
                            self.dist_fw.append(fl.dist)
                            print ("Gap joined on SG direct neighbours!")

            #HEURISTIC #2: indirect overlapping connection on SG

            #Now an indirect short-range connection
            if join_overlapping and not connected and mlc:
                conn_paths=self.ws.sg.find_all_paths_between(pn,n,0,30)
                if len(conn_paths)==1:
                    conn_path=conn_paths[0]
                    #print ("found connection path between nodes %d and %d! " %(pn,n), [x for x in conn_path.nodes])
                    ps=conn_path.get_sequence_size_fast()
                    pd=ps-2*199#len(self.ws.sg.nodes[abs(pn)].sequence)-len(self.ws.sg.nodes[abs(n)].sequence)
                    #print ("path size %d, estimated distance %d, %d mldg conections with median %s" %
                    #(ps,pd,len(mlc),median(mlc)))
                    if pd>min(median(mlc)-100,median(mlc)*.8) and pd<max(median(mlc)+100,median(mlc)*1.2):
                        connected=True
                        if pd<0:
                            s=s[:pd]
                        if pd>0:
                            s+=conn_path.get_sequence()[199:-199]
                        short_indirect+=1
                        self.joined_fw.append(True)
                        self.dist_fw.append(pd)
                        print ("Gap joined on SG short indirect connection!")
            #HEURISTIC #3: path on SG supported by nanopore reads

            if walk_with_nanopore and not connected:
                a=self.walk_with_nanopore(pn,n)

            #HEURISTIC #4: path on SG supported by nanopore reads

            if fill_with_nanopore and not connected:
                pass

            if not connected:
                n_stretches+=1
                if mlc:
                    d=median(mlc)
                else:
                    d=100
                unconnected_distances.append(d)
                #s+=''.join(['N' for x in range(200)])
                s+=''.join(['N' for x in range(max(200,int(d)))])
                self.joined_fw.append(False)
                self.dist_fw.append(d)
            if n>0:
                s+=self.ws.sg.nodes[n].sequence
            else:
                rnode=bsg.Node(self.ws.sg.nodes[-n].sequence)
                rnode.make_rc()
                s+=rnode.sequence
        #print("consensus includes %d direct connections, %d short patches, and %d N stretches" %
        #(direct_connections,short_indirect,n_stretches))
        #unconnected_distances.sort()
        #print("Distances for unconnected neighbours: ", unconnected_distances)
        return s
