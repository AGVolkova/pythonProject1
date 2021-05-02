def PairedGraph(Patterns): # То же, что и DeBruijnGraph, но на вход берет не строку, а кучу кусочков
    slovar = {}
    s = []
    for i in set(Patterns):
        a, b = i.split('|')[0], i.split('|')[1]
        for j in Patterns:
            c, d = j.split('|')[0], j.split('|')[1]
            if a[1:]==c[:-1] and b[1:]==d[:-1] and i not in slovar.keys():
                slovar[i]=j
            elif a[1:]==c[:-1] and b[1:]==d[:-1] and i in slovar.keys():
                slovar[i]=slovar[i]+','+j
    Graph=[]
    for i in slovar:
        Graph.append(i+' -> '+slovar[i])
    return Graph

def String_from_paired_reads(Graph, k, d): # Сюда подставлять граф, полученный функцией PairedGraph
    eu=Eulerian_Path(Graph)
    E=eu.split('->')
    Prefixes=[]
    Suffixes=[]
    for i in E:
        a, b = i.split('|')[0], i.split('|')[1]
        Prefixes.append(a)
        Suffixes.append(b)
    for i in range (k+d+1, len(PathToGenome(Prefixes))):
        if PathToGenome(Prefixes)[i]!=PathToGenome(Suffixes)[i-(k+d)]:
            print('there is no string')
    print(PathToGenome(Prefixes)+PathToGenome(Suffixes)[-(k+d):])

print('Train')
