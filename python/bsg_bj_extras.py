#!/usr/bin/env python3
import pybsg as bsg
import pylab

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

    def plot(self,read_id,signed_nodes=False,highlight_nodes=[]):
        self.plot_matches(self.lorm.filtered_read_mappings[read_id],signed_nodes,highlight_nodes)

    def plot_matches(self,matches,signed_nodes=False, highlight_nodes=[],light_highlight_nodes=[]):
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
                    pylab.plot(xpos,ypos,".-")

            last_y+=nlength
            pylab.plot([0,read_last_pos],[last_y,last_y],"k:")
        pylab.title('Filtered matches for read %d' % matches[0].read_id)
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







from graphviz import Graph
def plot_nodes_around(_nodes,ldg,radius,min_links):
    link_edges=[]
    initialnodes=[abs(x) for x in _nodes]
    s = Graph('structs', node_attr={'shape': 'plaintext'},engine='neato')
    s.attr(size='15,8',ratio='fill')
    nodes=initialnodes.copy()
    for x in range(radius):
        new_nodes=[]
        for n in nodes:
            for fnd in ldg.fw_neighbours_by_distance(n,min_links):
                if abs(fnd[1]) not in nodes: new_nodes.append(abs(fnd[1]))
            for bnd in ldg.fw_neighbours_by_distance(-n,min_links):
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
        for fnd in ldg.fw_neighbours_by_distance(n,min_links):
            if abs(fnd[1]) in nodes and abs(fnd[1])<=abs(n):
                v1='%d:-'%n
                if fnd[1]>0: v2='%d:+'%fnd[1]
                else: v2='%d:-'% -fnd[1]
                if (v1,v2) not in link_edges and (v2,v1) not in link_edges: link_edges.append((v1,v2))
        for bnd in ldg.fw_neighbours_by_distance(-n,min_links):
            if abs(bnd[1]) in nodes and abs(bnd[1])<=abs(n):
                v1='%d:+'%n
                if bnd[1]>0: v2='%d:+'%bnd[1]
                else: v2='%d:-'% -bnd[1]
                if (v1,v2) not in link_edges and (v2,v1) not in link_edges: link_edges.append((v1,v2))
    s.edges(link_edges)
    return s

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

#now check for every anchor how many reads connect to other lines!
