import random
import numpy as np
import pandas as pd
f=open('C:/Users/av/Desktop/Bioinformatics/2/dataset_205_5 (1).txt')

#print(f)
Data=f.read().splitlines()
#print(Dataset)

def StringComposition(S, k): #разбивает строку S на кусочки длины k, накладывающиеся друг на друга одним элементом
    l=[]
    for i in range(0, len(S)-k+1):
        print(S[i:i+k])
        l.append(S[i:i+k])


#StringComposition('ATTGCCTA', 3)

def PathToGenome(path): # Противоположно СтрингКомпозишн составляет строку из последовательных кусочков, пересекающихся на 1 элемент
    genome=path[0]
    for i in range(1, len(path)):
        genome=genome+path[i][-1]
    return genome

#print(PathToGenome(['ATT', 'TTG', 'TGC', 'GCC', 'CCT', 'CTA']))

def Overlap(Patterns): #Среди кучи кусочков находит пересечения и представляет в виде ААT->ATG, ATC
    l=list(set(Patterns))
    answer=[]
    for i in l:
        spisok=[]
        for j in Patterns:
            if i[1:]==j[:-1]:
                spisok.append(j)
        if spisok:
            answer.append(i+'->'+','.join(spisok))
            print(i+'->'+','.join(spisok))
    return answer

def All_Overlap(Patterns): #Среди кучи кусочков находит пересечения и представляет в виде ААT->ATG, ATC
    answer=[]
    for i in Patterns:
        spisok=[]
        for j in Patterns:
            if i[1:]==j[:-1]:
                spisok.append(j)
        if spisok:
            answer.append(i+'->'+','.join(spisok))
            print(i+'->'+','.join(spisok))
    return answer
print(All_Overlap(['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']))


def DeBruijnGraph(Text, k): # Строку представлет в виде пересекающихся кусочков длины К-1
    l = []
    for i in range(0, len(Text) - k + 1):
        l.append(Text[i:i + k])
    d = {}
    d[l[0][:-1]] = l[0][1:]
    a = Text[1:1 + k - 1]
    s = []
    for i in set(l):
        if i[1:] in s:
            continue
        for j in l[1:]:
            if i[1:] == j[:-1]:
                if i[1:] not in d.keys():
                    d[i[1:]] = j[1:]
                else:
                    d[i[1:]] = d[i[1:]] + ',' + j[1:]
            s.append(i[1:])

    for i in d:
        print(i + ' -> ' + d[i])

#DeBruijnGraph('ATTGGCTAAAGGCCCT', 4)


def DeBruijn(Patterns): # То же, что и DeBruijnGraph, но на вход берет не строку, а кучу кусочков
    d = {}
    s = []
    for i in set(Patterns):
        if i[:-1] in s:
            continue
        for j in Patterns:
            if i[:-1] + j[-1] == j and i[1:-1] == j[1:-1] and i[:-1] not in d.keys():
                d[i[:-1]] = j[1:]
            elif i[:-1] + j[-1] == j and i[1:-1] == j[1:-1] and i[:-1] in d.keys():
                d[i[:-1]] = d[i[:-1]] + ',' + j[1:]

            s.append(i[:-1])
    deBruijn = []
    for i in d:
        deBruijn.append(i + ' -> ' + d[i])
    return deBruijn

#print(DeBruijn(['ATT', 'GGC', 'TGG', 'GCC', 'AAA', 'CCC', 'CTA', 'TTG', 'GCT', 'TAA', 'AAG', 'AGG']))

def Eulerian_Path(Dataset):
    dictionary1={}
    dictionary2={}
    keys=[]
    values=[]
    for a in Dataset:
        u, v = a.split(' -> ')
        v=v.split(',')
        keys.append(u)
        for i in v:
            values.append(i)
        dictionary1[u]=v
    vse_uzly=[el for el in keys]
    for i in values:
        if i not in vse_uzly:
            vse_uzly.append(i)

    for i in vse_uzly:
        if i not in dictionary1.keys():
            dictionary2[i]=(0, values.count(i))
        else:
            dictionary2[i]=(len(dictionary1[i]), values.count(i))
    begin=0
    end=0
    for i in dictionary2.keys():
        if dictionary2[i][0]>dictionary2[i][1]:
            begin=i
        elif dictionary2[i][0]<dictionary2[i][1]:
            end=i

    all_edges=[]
    nod=vse_uzly
    for i in dictionary1:
        for j in dictionary1[i]:
            all_edges.append(i+'->'+j)
    all_edges.append(end+'->'+begin)

    dic=dictionary1
    if end not in dic.keys():
        dic[end]=[begin]
    else:
        dic[end]=dic[end]+[begin]
    n=nod
    Cycle=[]
    used_edges=[]

    Cycle.append(begin)
    n1=begin
    n2=random.choice(dic[n1])

    while n1+'->'+n2 not in used_edges:
        Cycle.append(n2)
        used_edges.append(n1+'->'+n2)
        dic[n1]=[el for el in dic[n1] if el !=n2]
        if dic[n1]==[]:
            n.remove(n1)
        if n2 not in dic.keys():
            break
        elif dic[n2]==[]:
            break
        else:
            n1, n2 = n2, random.choice(dic[n2])

        while n1+'->'+n2 not in used_edges:
            Cycle.append(n2)
            used_edges.append(n1+'->'+n2)
            dic[n1]=[el for el in dic[n1] if el !=n2]
            if dic[n1]==[]:
                n.remove(n1)
            if n2 not in dic.keys():
                break
            elif dic[n2]==[]:
                break
            else:
                n1, n2 = n2, random.choice(dic[n2])

    while len(used_edges)!=len(all_edges):
        unused_nodes=[]
        for i in Cycle:
            if i in n:
                unused_nodes.append(i)
        n1=random.choice(unused_nodes)
        i=Cycle.index(n1)
        difference=i
        new_edges=[]
        new_Cycle=[]
        new_Cycle.append(n1)
        while len(new_edges)!=len(used_edges):
            if len(new_Cycle)!=len(Cycle)-difference:
                new_edges.append(n1+'->'+Cycle[i+1])
                new_Cycle.append(Cycle[i+1])
                i=i+1
                n1=Cycle[i]
            else:
                new_edges.append(Cycle[0]+'->'+Cycle[1])
                new_Cycle.append(Cycle[1])
                n1, i=Cycle[1], 1
        Cycle=new_Cycle
        n2=random.choice(dic[n1])

        while n1+'->'+n2 not in used_edges:
            Cycle.append(n2)
            used_edges.append(n1+'->'+n2)
            dic[n1]=[el for el in dic[n1] if el !=n2]
            if dic[n1]==[]:
                n.remove(n1)
            if dic[n2]==[]:
                break
            else:
                n1, n2 = n2, random.choice(dic[n2])
    stroka='->'.join(Cycle[1:])
    s=stroka.split(end+'->'+begin)
    stroka_rearranged=begin+s[1]+'->'+s[0]+end
    return stroka_rearranged

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



def Find_contigs(Data):
    Dataset=All_Overlap(Data)
    dictionary1={}
    dictionary2={}
    keys=[]
    values=[]
    for a in Dataset:
         u, v = a.split('->')
         v=list(set(v.split(',')))
         print(v)
         keys.append(u)
         sp = []
         for i in v:
             if i not in sp:
                 sp.append(i)
             print(sp)
         values.append(','.join(sp))
         if u in dictionary1.keys():
             dictionary1[u] = dictionary1[u] + v
         else:
             dictionary1[u] = v

    values=(','.join(values)).split(',')
    print(dictionary1)
    vse_uzly=Data
    for i in vse_uzly:
        if i not in dictionary1.keys():
            dictionary2[i]=(0, values.count(i))
        else:
            dictionary2[i]=(len(dictionary1[i]), values.count(i))
    print(dictionary2)
    values=(','.join(values)).split(',')
    Single_in_single_out=[i for i in Data if dictionary2[i][0]==1 and dictionary2[i][1] == 1 and i[:-1]!=i[1:]]
    Non_branching=[]
    Last=[i for i in Data if dictionary2[i][0]==0]
    Used_last=[]
    for i in Data:
        if i not in Single_in_single_out and dictionary2[i][0]>0:
            for j in dictionary1[i]:
                Nonbranchingpath = i
                while j in Single_in_single_out:
                    Nonbranchingpath=Nonbranchingpath+j[-1]
                    j=dictionary1[j].copy().pop()
                if j in Last and len(dictionary1[i])==1:
                    Non_branching.append(Nonbranchingpath+j[-1])
                    Used_last.append(j)
                else:
                    Non_branching.append(Nonbranchingpath)
                break
    for i in Last:
        if i not in Used_last:
            Non_branching.append(i)
    for i in Single_in_single_out:
        Cycle=[]
        Cycle.append(i)

        while dictionary1[i].copy().pop() in Single_in_single_out:
            i=dictionary1[i].copy().pop()
            Cycle.append(i)

            if dictionary1[i].copy().pop() in Cycle:
                Non_branching.append(PathToGenome(Cycle))
                break
        break

    for i in Non_branching:
        print(i)


