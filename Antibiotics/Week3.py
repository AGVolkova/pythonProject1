RNA_codon_table=open('C:/Users/AVolkova/Desktop/Биоинформатика/RNA_codon_table.txt')

Data=RNA_codon_table.read().splitlines()
Codons={}
for i in Data:
    Codons[i[0:3]]=i[-1]
print(Codons)

def RNA_to_peptide(RNA):
    Peptide=Codons[RNA[0:3]]
    for i in range(3, len(RNA), 3):
        if Codons[RNA[i:i+3]]==' ':
            break
        else:
            Peptide=Peptide+Codons[RNA[i:i+3]]
        i=i+3
    return Peptide

DNA_reversion={'A':'T', 'T':'A', 'G':'C', 'C':'G'}

def DNA_to_RNA(DNA):
    reversed_DNA=DNA_reversion[DNA[0]]
    for i in DNA[1:]:
        reversed_DNA=reversed_DNA+DNA_reversion[i]
    reversed_DNA=reversed_DNA[::-1]
    RNA=[]
    if DNA[0]=='T':
        RNA1='U'
    else:
        RNA1=DNA[0]
    if reversed_DNA[0]=='T':
        RNA2='U'
    else:
        RNA2=reversed_DNA[0]
    for i in DNA[1:]:
        if i=='T':
            RNA1=RNA1+'U'
        else: RNA1=RNA1+i
    RNA.append(RNA1)
    for i in reversed_DNA[1:]:
        if i == 'T':
            RNA2 = RNA2 + 'U'
        else:
            RNA2 = RNA2 + i
    RNA.append(RNA2)
    return RNA

def Peptide_Encoding_Problem(dna, Peptide):
    print(Peptide)
    Length=len(Peptide)*3
    Genes=[]
    for i in range(0, len(dna)-Length+1):
        if RNA_to_peptide(DNA_to_RNA(dna[i: i+Length])[0])==Peptide or RNA_to_peptide(DNA_to_RNA(dna[i: i+Length])[1])==Peptide:
            Genes.append(dna[i: i+Length])
    return Genes

Bacillus_brevis_Genome=open('C:/Users/AVolkova/Desktop/Биоинформатика/Bacillus_brevis.txt')
Bacillus=Bacillus_brevis_Genome.read().replace('\n', '')
#print(len(Peptide_Encoding_Problem(Bacillus, 'VKLFPWFNQY')))

AminoAcids_Mass={'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
All_AminoAcids=['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
def LinearSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass=[0]
    for i in range(1, len(Peptide)+1):
        for j in Alphabet:
            if Peptide[i-1]==j:
                PrefixMass.append(PrefixMass[i-1]+AminoAcidMass[j])
    print(PrefixMass)
    LinearSpectrum=[0]
    for i in range(0, len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    Sorted_list=sorted(LinearSpectrum)
    return Sorted_list
# a=[str(i) for i in LinearSpectrum('LEQN', All_AminoAcids, AminoAcids_Mass)]
# print(" ".join(a))

def CyclicSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass=[0]
    for i in range(1, len(Peptide)+1):
        for j in Alphabet:
            if Peptide[i-1]==j:
                PrefixMass.append(PrefixMass[i-1]+AminoAcidMass[j])
    peptidemass=max(PrefixMass)
    print(PrefixMass)
    CyclicSpectrum=[0]
    for i in range(0, len(Peptide)):
        for j in range(i+1, len(Peptide)+1):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i >0 and j< (len(Peptide)):
                CyclicSpectrum.append(peptidemass-(PrefixMass[j]-PrefixMass[i]))
                print('Yes')
    Sorted_list=sorted(CyclicSpectrum)
    return Sorted_list


a=[str(i) for i in CyclicSpectrum('QCMQEQCIDAHI', All_AminoAcids, AminoAcids_Mass)]
print(" ".join(a))