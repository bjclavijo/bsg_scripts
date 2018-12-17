#!/usr/bin/env python3
import pybsg as bsg
import pylab
from graphviz import Graph
import numpy as np
from statistics import mean, median

class LDGStats(object):
    """docstring for LDGStats."""
    def __init__(self, _ldg):
        self.ldg = _ldg

    def get_line_lenght(self,line):
        totalbp=len(self.ldg.sg.nodes[abs(line[0])].sequence)
        for i in range(len(line)-1):
            totalbp+=[x.dist for x in self.ldg.get_fw_links(line[i]) if x.dest==line[i+1]][0]
            totalbp+=len(self.ldg.sg.nodes[abs(line[i+1])].sequence)
        return totalbp

    def print_line_contiguity_stats(self,_min_nodes=2):
        lines=self.ldg.get_all_lines(_min_nodes)
        lines_lenghts=[self.get_line_lenght(x) for x in lines]
        lines_lenghts.sort(reverse=True)
        total=sum(lines_lenghts)
        print("Total bp in lines: %d" % total)
        t=0
        checkpoints=[0,10,20,30,40,50,60,80,90,100,200]
        ncp=0
        l=0
        for c in lines_lenghts:
            t+=c
            l+=1
            while t>=total*checkpoints[ncp]/100:
                print("N%d: %dbp, contig #%d"%(checkpoints[ncp],c,l))
                ncp+=1

class LRMatchesPlotter(object):
    """docstring for LRMatchesPlotter."""
    def __init__(self, lorm, fs):
        self.lorm=lorm
        self.fs=fs

    def plot_unfiltered(self,read_id,signed_nodes=False,highlight_nodes=[],filter_nodes=[]):
        m=[ x for x in self.lorm.get_raw_mappings_from_read(read_id) if not filter_nodes or abs(x.node) in filter_nodes]
        self.plot_matches(m,signed_nodes,highlight_nodes,title="RAW matches for read %d" % read_id)
        pylab.title
    def plot_unfiltered_in_node(self,node_id,signed_nodes=False,highlight_nodes=[]):
        reads=[]
        for x in self.lorm.mappings:
            if abs(x.node)==node_id and ((not reads) or x.read_id != reads[-1]):
                reads.append(x.read_id)
        for read_id in reads:
            matches=[x for x in self.lorm.mappings if x.read_id==read_id]
            if len(matches)>50: continue
            self.plot_matches(matches,signed_nodes,highlight_nodes=[node_id])


    def plot(self,read_id,signed_nodes=False,highlight_nodes=[]):
        self.plot_matches(self.lorm.filtered_read_mappings[read_id],signed_nodes,highlight_nodes,title="Filtered matches for read %d" % read_id)

    def plot_matches(self,matches,signed_nodes=False, highlight_nodes=[],light_highlight_nodes=[],title="",xlim=0):
        pylab.figure(figsize=self.fs)
        last_y=0
        nodes=[]
        if not matches: return
        tickpos=[]
        ticklabel=[]
        for m in matches:
            if abs(m.node) not in [abs(x) for x in nodes]:
                if signed_nodes: nodes.append(m.node)
                else: nodes.append(abs(m.node))
        read_last_pos=max([x.qEnd for x in matches])
        read_last_pos=max(xlim,read_last_pos)
        for n in nodes:
            nlength=len(self.lorm.getSequenceGraph().nodes[abs(n)].sequence)
            tickpos.append(last_y+nlength/2)
            ticklabel.append("%d (%d bp)"% (n,nlength))
            if abs(n) in highlight_nodes:
                pylab.axhspan(last_y, last_y+nlength, facecolor='y', alpha=0.25)
            if abs(n) in light_highlight_nodes:
                pylab.axhspan(last_y, last_y+nlength, facecolor='y', alpha=0.15)
            for m in matches:
                if abs(m.node)==abs(n):
                    xpos=(m.qStart,m.qEnd)
                    if m.node==n:
                        ypos=(m.nStart+last_y,m.nEnd+last_y)
                    else:
                        ypos=(nlength-m.nStart-15+last_y,nlength-m.nEnd-15+last_y)
                    pylab.plot(xpos,ypos,"x-")
                    s="%d%%" % (m.score*100.0/(m.nEnd-m.nStart))
                    pylab.text(mean(xpos),mean(ypos),s,ha="center", va="center",bbox=dict(boxstyle="round", color="white",alpha=.5))

            last_y+=nlength
            pylab.plot([0,read_last_pos],[last_y,last_y],"k:")
        if not title:title='Matches for read %d' % matches[0].read_id
        pylab.title(title)
        pylab.xlabel('Read position')
        pylab.ylabel('Node')
        pylab.xlim(left=0,right=read_last_pos)
        pylab.ylim(bottom=0,top=last_y)
        pylab.yticks(tickpos, ticklabel)
        pylab.show()

    def plot_all_matches_for_node(self,node,signed_nodes=False,highlight_all=[]):
        for frm in self.lorm.filtered_read_mappings:
            if [x for x in frm if abs(x.node)==abs(node)]:
                self.plot_matches(frm,signed_nodes,light_highlight_nodes=highlight_all,highlight_nodes=[abs(node)])

    def plot_all_matches_for_nodes(self,nodes,signed_nodes=False,highlight_all=[]):
        for frm in self.lorm.filtered_read_mappings:
            if all(n in [abs(x.node) for x in frm] for n in nodes):
                self.plot_matches(frm,signed_nodes,light_highlight_nodes=highlight_all,highlight_nodes=[abs(x) for x in nodes])


class LRMatchesAnalyser(object):
    """docstring for LRMatchesPlotter.
    Read classitication:
    0 - Not analysed
    1 - """
    def __init__(self, lorm, allowed_overlap=500):
        self.lorm=lorm
        self.allowed_overlap=allowed_overlap
        self.read_class_names={
        0: 'OK',
        1: 'No mappings',
        2: 'Significant overlaps',
        10: 'Not analysed',
        }

    def classify_reads(self):
        self.read_class = np.zeros((len(self.lorm.filtered_read_mappings),), dtype=np.int8)
        rid=-1
        for frm in self.lorm.filtered_read_mappings:
            rid+=1
            c=0
            if not frm:
                c=1
            else:
                for i in range(len(frm)-1):
                    if frm[i].qEnd-self.allowed_overlap>frm[i+1].qStart:
                        c=2
                        break
            self.read_class[rid]=c

    def report_classification(self):
        unique, counts = np.unique(self.read_class, return_counts=True)
        for uc in zip(unique, counts):
            print(self.read_class_names[uc[0]]," -> ",uc[1],"(%.2f%%)" % (int(uc[1])*100.0/len(self.read_class)))

class LDGPlotter(object):
    """docstring for LDGPlotter."""
    def __init__(self, ldg):
        self.ldg = ldg

    def plot_around_nodes(self,_nodes,radius=10,min_links=1):
        link_edges=[]
        initialnodes=[abs(x) for x in _nodes]
        s = Graph('structs', node_attr={'shape': 'plaintext'},engine='neato')
        s.attr(size='15,8',ratio='fill')
        nodes=initialnodes.copy()
        for x in range(radius):
            new_nodes=[]
            for n in nodes:
                for fnd in self.ldg.fw_neighbours_by_distance(n,min_links):
                    if abs(fnd[1]) not in nodes: new_nodes.append(abs(fnd[1]))
                for bnd in self.ldg.fw_neighbours_by_distance(-n,min_links):
                    if abs(bnd[1]) not in nodes: new_nodes.append(abs(bnd[1]))
            nodes+=new_nodes
            if x<radius: new_nodes=[]

        for n in nodes:
            if n in initialnodes: nodecolor='bgcolor="yellow"'
            elif n in new_nodes: nodecolor='bgcolor="cyan"'
            else: nodecolor=''
            s.node(str(n), '''<
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
          <TR>
            <TD PORT="+" bgcolor="darkgreen"> </TD>
            <TD %s>%d</TD>
            <TD PORT="-" bgcolor="red"> </TD>
          </TR>
        </TABLE>>''' % (nodecolor,n))
            for fnd in self.ldg.fw_neighbours_by_distance(n,min_links):
                if abs(fnd[1]) in nodes and abs(fnd[1])<=abs(n):
                    v1='%d:-'%n
                    if fnd[1]>0: v2='%d:+'%fnd[1]
                    else: v2='%d:-'% -fnd[1]
                    if (v1,v2) not in link_edges and (v2,v1) not in link_edges: link_edges.append((v1,v2))
            for bnd in self.ldg.fw_neighbours_by_distance(-n,min_links):
                if abs(bnd[1]) in nodes and abs(bnd[1])<=abs(n):
                    v1='%d:+'%n
                    if bnd[1]>0: v2='%d:+'%bnd[1]
                    else: v2='%d:-'% -bnd[1]
                    if (v1,v2) not in link_edges and (v2,v1) not in link_edges: link_edges.append((v1,v2))
        s.edges(link_edges)
        return s

    def plot_linked_space_around_node(self, node, clean_ldg=None):
        if clean_ldg is None: clean_ldg=self.ldg
        read_first_fconnection={}
        for l in clean_ldg.get_fw_links(node):
            if l.read_id not in read_first_fconnection or read_first_fconnection[l.read_id].dist<l.dist:
                read_first_fconnection[l.read_id]=l
        read_first_bconnection={}
        for l in clean_ldg.get_bw_links(node):
            if l.read_id not in read_first_bconnection or read_first_bconnection[l.read_id].dist<l.dist:
                read_first_bconnection[l.read_id]=l

        #create lists of distances for every node fw and bw and plot as histograms
        fnodes=set([x.dest for x in read_first_fconnection.values()])
        bnodes=set([x.dest for x in read_first_bconnection.values()])
        fdists={}
        bdists={}
        for l in self.ldg.get_fw_links(node):
            if l.dest in fnodes:
                if l.dest not in fdists:
                    fdists[l.dest]=[]
                fdists[l.dest].append(l.dist)

        for l in self.ldg.get_bw_links(node):
            if l.dest in bnodes:
                if l.dest not in bdists:
                    bdists[l.dest]=[]
                bdists[l.dest].append(l.dist)
        #
        #plot_dist_hists_from_dict(fdists)
        self.plot_dist_hists_from_dict_as_used_space(bdists)
        self.plot_dist_hists_from_dict_as_used_space(fdists)

    def plot_dist_hists_from_dict_as_used_space(self,ddict):
        BINSIZE=250
        nids=list(ddict.keys())
        nids.sort(key=lambda x:min([ddict[x]]))
        dists=[ddict[x] for x in nids]
        labels=["%d ( %d bp )" % (x, len(self.ldg.sg.nodes[abs(x)].sequence)) for x in nids]
        #compute how many bins and the limits:
        upper_limmit=(max([max(x) for x in dists])+5000)//BINSIZE*BINSIZE
        lower_limmit=(200+BINSIZE)//BINSIZE*BINSIZE
        bin_count=(upper_limmit-lower_limmit+1)//BINSIZE
        occ=[]
        for ni in nids:
            oc=[]
            nsize=len(self.ldg.sg.nodes[abs(ni)].sequence)
            for d in ddict[ni]:
                for x in range(0,nsize,BINSIZE):
                    oc.append(d+x+BINSIZE/2)
            occ.append(oc)
        pylab.figure(figsize=(10,3))
        if dists:
            pylab.hist(occ,histtype='barstacked',label=labels,bins=bin_count,range=(lower_limmit,upper_limmit))
        pylab.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        pylab.show()

def deselect_selfloops(u,ldg):
    for x in ldg.find_self_loops():
        u.selected_nodes[abs(x)]=False

def deselect_tips(u,ldg):
    removed=0
    for x in ldg.find_tips():
        removed+=1
        u.selected_nodes[abs(x)]=False
    return removed

def print_cardinalities(ldg):
    cardd={}
    for x in ldg.links:
        a=sum([1 for y in x if y.source<0])
        b=sum([1 for y in x if y.source>0])
        c=(a,b)
        if b<a: c=(b,a)
        if c not in cardd: cardd[c]=0
        cardd[c]+=1
    cc=[]
    for c in cardd.keys():
        cc.append((c,cardd[c]))
    cc.sort()
    for c in cc: print (c[0],c[1])

def print_lines_cardinalities(ldg,min_nodes=2):
    cardd={}
    for x in ldg.get_all_lines(min_nodes):
        a=len(ldg.get_fw_links(x[-1]))
        b=len(ldg.get_bw_links(x[0]))
        c=(a,b)
        if b<a: c=(b,a)
        if c not in cardd: cardd[c]=[0,0]
        cardd[c][0]+=1
        cardd[c][1]+=get_line_lenght(ldg,x)
    cc=[]
    for c in cardd.keys():
        cc.append((c,cardd[c]))
    cc.sort()
    for c in cc: print (c[0],c[1])

def reconnect_bubble_sides(ldg,mldg):
    from statistics import median
    ab=0
    ba=0
    both=0
    none=0
    cons=[]
    for x in ldg.find_bubbles(0,40000):
        #print(aldg2.links[abs(x[0])])
        #print(aldg2.links[abs(x[1])])
        #print(mldg.links[abs(x[0])])
        cab=mldg.are_connected(-x[0],x[1])
        cba=mldg.are_connected(-x[1],x[0])
        #print(cab,cba)
        c1=0
        c2=0
        if cba and not cab:
            ba+=1
            c1=-x[1]
            c2=x[0]
        elif cab and not cba:
            ab+=1
            c1=-x[0]
            c2=x[1]
        elif cab and cba: both+=1
        else: none+=1
        if c1 and c2:
            d=int(median([x.dist for x in mldg.links[abs(c1)] if x.source==c1 and x.dest==c2]))
            cons.append((c1,c2,d))
    print(ab,ba,both,none)
    for x in cons: ldg.add_link(x[0],x[1],x[2])

def describe_fw(node,min_links=3):
    for x in mldg.fw_neighbours_by_distance(node,min_links):
        print(x,end="")
        if u.selected_nodes[abs(x[1])]:  print(' ANCHOR',end='')
        for i in range(len(lines)):
            if abs(x[1]) in [abs(y) for y in lines[i]]:
                print(' found in line %d' % (i),end='')
        print()

def make_and_simplify_linkage():
    u=bsg.LinkageUntangler(ws)
    u.select_multi_linkage_linear_anchors(mldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    deselect_selfloops(u,ldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    for step in range(10):
        print("Tip removal round %d"%step)
        r=deselect_tips(u,ldg)
        print('%d tips removed'% r)
        if r==0: break
        ldg=u.make_nextselected_linkage(mldg)
        ldg.remove_transitive_links(10)
    reconnect_bubble_sides(ldg,mldg)
    ldg.remove_transitive_links(10)
    for i in range(len(ldg.links)):
        if (len([x for x in ldg.links[i] if x.source==i]) >1) or (len([x for x in ldg.links[i] if x.source==-i]) >1): u.selected_nodes[i]=False
    u.report_node_selection()
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    deselect_selfloops(u,ldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    for step in range(10):
        print("Tip removal round %d"%step)
        r=deselect_tips(u,ldg)
        print('%d tips removed'% r)
        if r==0: break
        ldg=u.make_nextselected_linkage(mldg)
        ldg.remove_transitive_links(10)
    reconnect_bubble_sides(ldg,mldg)
    ldg.remove_transitive_links(10)
    for i in range(len(ldg.links)):
        if (len([x for x in ldg.links[i] if x.source==i]) >1) or (len([x for x in ldg.links[i] if x.source==-i]) >1): u.selected_nodes[i]=False
    u.report_node_selection()
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    deselect_selfloops(u,ldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    for step in range(10):
        print("Tip removal round %d"%step)
        r=deselect_tips(u,ldg)
        print('%d tips removed'% r)
        if r==0: break
        ldg=u.make_nextselected_linkage(mldg)
        ldg.remove_transitive_links(10)
    reconnect_bubble_sides(ldg,mldg)
    ldg.remove_transitive_links(10)
    for i in range(len(ldg.links)):
        if (len([x for x in ldg.links[i] if x.source==i]) >1) or (len([x for x in ldg.links[i] if x.source==-i]) >1): u.selected_nodes[i]=False
    u.report_node_selection()
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    deselect_selfloops(u,ldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    for step in range(10):
        print("Tip removal round %d"%step)
        r=deselect_tips(u,ldg)
        print('%d tips removed'% r)
        if r==0: break
        ldg=u.make_nextselected_linkage(mldg)
        ldg.remove_transitive_links(10)
    reconnect_bubble_sides(ldg,mldg)
    ldg.remove_transitive_links(10)
    for i in range(len(ldg.links)):
        if (len([x for x in ldg.links[i] if x.source==i]) >1) or (len([x for x in ldg.links[i] if x.source==-i]) >1): u.selected_nodes[i]=False
    u.report_node_selection()
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    deselect_selfloops(u,ldg)
    ldg=u.make_nextselected_linkage(mldg)
    ldg.remove_transitive_links(10)
    for step in range(10):
        print("Tip removal round %d"%step)
        r=deselect_tips(u,ldg)
        print('%d tips removed'% r)
        if r==0: break
        ldg=u.make_nextselected_linkage(mldg)
        ldg.remove_transitive_links(10)
    reconnect_bubble_sides(ldg,mldg)
    ldg.remove_transitive_links(10)
    return ldg

def get_all_reads_between_as(ns,lorm,u,lines,do_print=True):
    read_ids=[]
    for i in lorm.reads_in_node[abs(ns[0])]:
        if not all(i in lorm.reads_in_node[abs(x)] for x in ns): continue

        nodes=[x.node for x in lorm.filtered_read_mappings[i]]
        if all(x in nodes for x in ns): read_ids.append(i)
        elif all(-x in nodes for x in ns): read_ids.append(-i)
        else: continue
        if not do_print: continue
        if read_ids[-1]>0:
            for x in lorm.filtered_read_mappings[i]:
                nm=""
                if u.selected_nodes[abs(x.node)]: nm+=" -=A=-"
                for li in range(len(lines)):
                    if x.node in lines[li] or -x.node in lines[li]: nm+=" line #%d"% li
                if x.node in ns: nm+=" <<--- Node in list"
                print(x,nm)

        if read_ids[-1]<0:
            for x in lorm.filtered_read_mappings[i][::-1]:
                nm=""
                if u.selected_nodes[abs(x.node)]: nm+=" -=A=-"
                for li in range(len(lines)):
                    if x.node in lines[li] or -x.node in lines[li]: nm+=" line #%d"% li
                if -x.node in ns: nm+=" <<--- Node in list"
                print(x,nm)
        print()
    return read_ids

def reads_in_lines(lorm,lines):
    ril=[set() for x in lines]
    for lid in range(len(lines)):
        for n in lines[lid]:
            for x in lorm.reads_in_node[abs(n)]: ril[lid].add(x)
    return ril

def reads_to_lines(lorm,reads_in_lines):
    rtl=[set() for x in range(len(lorm.filtered_read_mappings))]
    for lrsi in range(len(reads_in_lines)):
        for r in reads_in_lines[lrsi]:
            rtl[r].add(lrsi)
    return rtl

def kmer_to_number(seq):
    fkmer=0
    rkmer=0
    i1=0
    i2=len(seq)-1
    for x in seq:
        if x=='A':
            fkmer+=0<<2*i1
            rkmer+=3<<2*i2
        elif x=='C':
            fkmer+=1<<2*i1
            rkmer+=2<<2*i2
        elif x=='G':
            fkmer+=2<<2*i1
            rkmer+=1<<2*i2
        elif x=='T':
            fkmer+=3<<2*i1
            rkmer+=0<<2*i2
        else: return None
        i1+=1
        i2-=1
    if fkmer<rkmer: return fkmer
    else: return -rkmer

def sequence_to_kmer_numbers(seq,k,positive=True):
    numbers=[]
    for kstart in range(len(seq)-k):
        if positive:
            numbers.append(abs(kmer_to_number(seq[kstart:kstart+k])))
        else:
            numbers.append(kmer_to_number(seq[kstart:kstart+k]))
    return numbers;
class Nano10xDetailedAnalysis(object):
    """uses a LinkedReadMapper to select best haplotype solution in a LongRead."""
    def __init__(self, sg, lirm, lorm, k):
        super(Nano10xDetailedAnalysis, self).__init__()
        self.sg=sg
        self.lirm = lirm
        self.lorm = lorm
        self.k=k
        self.blrsg=bsg.BufferedSequenceGetter(lorm.datastore,10000000,1000000)

    def load_read_data(self,read_id):
        "sets mappings, nodes and nodesets, and read_seq"
        self.read_id=read_id
        self.mappings=[]
        self.nodes=set()
        max_map_end=0
        for i in range(len(self.lorm.mappings)):
            if self.lorm.mappings[i].read_id>read_id: break
            if self.lorm.mappings[i].read_id==read_id:
                self.mappings.append(self.lorm.mappings[i])
                self.nodes.add(abs(self.mappings[-1].node))
                max_map_end=max(max_map_end,self.mappings[-1].qEnd)
        print("creating all nodesets for read %d, with nodes: %s"%(read_id,self.nodes))
        self.nodesets=set()
        for node in self.nodes:
            nodeset=[n.node for n in self.lirm.tag_neighbours[node] if n.score>0.05 and n.node in self.nodes]
            nodeset.sort()
            self.nodesets.add(tuple(nodeset))
        self.read_seq=self.blrsg.get_read_sequence(read_id)

    def plot_all_nodesets(min_cv):
        mp=LRMatchesPlotter(self.lorm,(15,5))
        for ns in self.nodesets:
            nm=[m for m in self.mappings if abs(m.node) in ns]
            if sum([m.qEnd-m.qStart for m in nm])>min_cv*max_map_end:
                mp.plot_matches(nm,xlim=max_map_end)

    def index_readwin_kmers(self, win_size, win_step):
        self.kidx={}
        self.win_size=win_size
        self.win_step=win_step
        win=0
        for wstart in range(0,len(self.read_seq)-win_size,win_step):
            wseq=self.read_seq[wstart:wstart+win_size]
            for kstart in range(len(wseq)-self.k):
                kn=abs(kmer_to_number(wseq[kstart:kstart+self.k]))
                if kn not in self.kidx: self.kidx[kn]=set()
                self.kidx[kn].add(win)
            win+=1

    def find_winners_per_window(self):
        win_votes=[{} for x in range((len(self.read_seq)-self.win_size)//self.win_step+1)] #every win has a dict of nodes with their votes
        for n in self.nodes:
            for kn in sequence_to_kmer_numbers(self.sg.nodes[n].sequence,self.k):
                if kn in self.kidx:
                    for w in self.kidx[kn]:
                        if n not in win_votes[w]: win_votes[w][n]=0
                        win_votes[w][n]+=1
        self.node_wins={}
        for n in self.nodes: self.node_wins[n]=0
        for w in range(len(win_votes)):
            #print("\nTop 5 nodes for win #%d (%d:%d): "%(w,w*self.win_step,w*self.win_step+self.win_size-1))
            votes=[(win_votes[w][n],n) for n in win_votes[w].keys()]
            votes.sort(reverse=True)
            for v in votes:
                if v[0]<votes[0][0] or v[0]<.2*self.win_size:break
                self.node_wins[v[1]]+=1
            #for v in votes[:5]:
            #    print ("%d: %d votes" % (v[1],v[0]))

    def score_nodesets(self):
        nss=[]
        for ns in self.nodesets:
            tvotes=sum([self.node_wins[n] for n in ns])
            nss.append([tvotes,ns])
        nss.sort(reverse=True)
        return nss

class LDGLineConsensus(object):
    """docstring for LDGLineConsensus."""
    def __init__(self, ldg, mldg, ws):
        self.ldg = ldg
        self.mldg = mldg
        self.ws = ws

    def get_lines(self,min_components=1,min_size=10000):
        self.lines=self.ldg.get_all_lines(min_components,min_size)

    def line_consensus(self, line_number):
        line=self.lines[line_number]
        completed_line=self.add_intermediate_nodes(line)
        return self.consensus_from_line(completed_line)

    def add_intermediate_nodes(self, line):
        return line

    def consensus_from_line(self, line):
        #print ("generating consensus for line ",line)
        n_stretches=0
        direct_connections=0
        short_indirect=0
        s=""
        #first node, no overlap possible to any previous
        n=line[0]
        if n>0:
            s+=self.ws.sg.nodes[n].sequence
        else:
            rnode=bsg.Node(self.ws.sg.nodes[-n].sequence)
            rnode.make_rc()
            s+=rnode.sequence
        #all other nodes, check for overlap/route to previous
        for ni in range(1,len(line)):
            n=line[ni]
            pn=line[ni-1]
            connected=False
            mlc=[x.dist for x in self.mldg.get_fw_links(pn) if x.dest==n]
            #first try a direct connection
            if median(mlc)<200:
                for fl in self.ws.sg.get_fw_links(pn):
                    if fl.dest==n:
                        #print("nodes %d and %d are connected on sg, distance %d !" % (pn,n,fl.dist))
                        if fl.dist<0:
                            s=s[:fl.dist]
                            connected=True
                            direct_connections+=1
            #Now an indirect short-range connection
            if not connected:
                conn_paths=self.ws.sg.find_all_paths_between(pn,n,max(2*max(mlc)+2*199,500),30)
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
            if not connected:
                n_stretches+=1
                #s+=''.join(['N' for x in range(200)])
                s+=''.join(['N' for x in range(max(200,int(median(mlc))))])
            if n>0:
                s+=self.ws.sg.nodes[n].sequence
            else:
                rnode=bsg.Node(self.ws.sg.nodes[-n].sequence)
                rnode.make_rc()
                s+=rnode.sequence
        print("consensus includes %d direct connections, %d short patches, and %d N stretches" %
        (direct_connections,short_indirect,n_stretches))
        return s
