import requests
#from bs4 import BeautifulSoup

def NucleotideCount(s):
    A=0
    C=0
    G=0
    T=0
    for i in s:
        if i == 'A':
            A+=1
        if i == 'C':
            C += 1
        if i == 'G':
            G += 1
        if i == 'T':
            T += 1
    file_name = "result.txt"
    with open(file_name, "w") as file:
        # Write data to the file
        file.write(str(A) + " " + str(C) + " "+ str(G) + " " + str(T))

def RNAtranscribe(file_name):
    with open(file_name, "r") as file:
        RNA = file.read()
    DNA = RNA.replace('U', 'T')
    with open("dna.txt", "w") as file:
        file.write(DNA)

def DNAtranscribe(file_name):
    with open(file_name, "r") as file:
        DNA = file.read()
    RNA=DNA.replace('T','U')
    with open("rna.txt", "w") as file:
        file.write(RNA)
def reverseComplement(file_name):
    with open(file_name, "r") as file:
        DNA = file.read()
    reverse=DNA[::-1]
    revcompl=''
    for i in reverse:
        if i == 'A':
            revcompl+='T'
        if i == 'T':
            revcompl+='A'
        if i == 'G':
            revcompl+='C'
        if i == 'C':
            revcompl+='G'

    with open("reversecomlement.txt", "w") as file:
        file.write(revcompl)

def RabbitFibonacci(n,k):
    #n=5 #months
    #k=3 #pairofrabbits/birth
    youngcount=1
    adultcount=0
    for i in range(n-1):
        d=adultcount*k
        adultcount+=youngcount
        youngcount=d
        print(youngcount+adultcount)
def RabbitFibonacciMortal(n,l):
    #n=5 #months
    #k=3 #pairofrabbits/birth
    #l: life expectancy in months
    list=[]
    for i in range(l):
        list.append(0)
    list[0]=1
    for i in range(n-1):
        adultcount=sum(list)-list[0]
        list[l-1]=0
        for i in range(l-1,0,-1):
            list[i]=list[i-1]
        list[0]=0
        list[0]+=adultcount
        #print(list)
        print(sum(list))
def GCpercent(file_name):
    codes=[]
    FASTAs=[]
    GCcontent=[]
    k=-1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                k+=1
                codes.append(line.strip())
                FASTAs.append('')
            if line.startswith('>')==False:
                FASTAs[k]+= line.strip()

    for i in FASTAs:
        GC = 0
        for d in i:
            if d == 'C':
                GC += 1
            if d == 'G':
                GC += 1
        GCcontent.append((GC/len(i))*100)

    maxGC=max(GCcontent)
    for i in range(len(GCcontent)):
        if GCcontent[i] == maxGC:
            GCindex=i
    print(codes)
    print(FASTAs)
    print( codes[GCindex], maxGC)
    with open("GCcontent.txt", "w") as file:
        file.write((str(codes[GCindex]).strip('>')) +'\n'+str(round(maxGC,6)))
def findMotif(file_name):
    strings=[]
    with open(file_name, "r") as file:
        for line in file:
            strings.append(line.strip())
    s=strings[0]
    t=strings[1]
    slices=[]

    result = ''
    for i in range(len(s)-len(t)):
        slices.append(s[i:i+len(t)])
        if slices[i]==t:
            result+=str(i+1) + ' '

    with open("findMotif.txt", "w") as file:
        file.write(result.strip())

def findPalindromes(file_name):
    FASTAs=[]
    k=-1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                k+=1
                code=line.strip()
                FASTAs.append('')
            if line.startswith('>')==False:
                FASTAs[k]+= line.strip()
    DNA=FASTAs[0]
    compl = ''
    for i in DNA:
        if i == 'A':
            compl += 'T'
        if i == 'T':
            compl += 'A'
        if i == 'G':
            compl += 'C'
        if i == 'C':
            compl += 'G'
    print(DNA + '\n' + compl)
    poz=[]
    length=[]
    print(len(DNA))
    for i in range(len(DNA)-1):
        if DNA[i] == compl[i+1]:
            d=0
            while (i-d)>=0 and (i+d+1)<len(DNA) and DNA[i-d]==compl[i+1+d]:
                d+=1
                print(i, d)
                if 12 >= 2*d >=4:
                    poz.append(i-d+2)
                    length.append(2*d)
    print(poz, length)
    with open("findPalindromes.txt", "w") as file:
        for i in range(len(poz)):
            file.write(str(poz[i]) + '\t' + str(length[i]) + '\t' + '\n')

def countMutations(file_name):
    strings=[]
    with open(file_name, "r") as file:
        for line in file:
            strings.append(line.strip())
    s=strings[0]
    t=strings[1]
    count=0
    for i in range(len(s)):
        if s[i] != t[i]:
            count+=1
    print(count)

def translateRNA(file_name):
    codons = []
    triplets = []
    with open('codon table.txt', "r") as file:
        codon_table=file.read()
    ct=codon_table.replace("\n", "      ")
    with open('ct.txt', "w") as file:
        file.write(ct)
    with open('ct.txt', "r") as file:

        while True:
            char=file.read(3)
            if not char: break
            triplets.append(char.strip())
            char=file.read(8)
            if not char: break
            codons.append(char.strip())
    with open(file_name, "r") as file:
        RNA=file.read()
    protein =''
    slice=''
    k=0
    for i in RNA:
        slice+=i
        if len(slice) == 3:
            index = triplets.index(slice)
            if codons[index]=='Stop':
                break
            protein+=codons[index]
            slice=''
    print(codons)
    print(triplets)
    print(protein)

def WebScrape(url):
    response = requests.get(url).text
    return response

def findN_glycMotif(file_name):
    IDs=[]
    with open(file_name,'r') as file:
        for line in file:
            IDs.append(line.strip())
    print(IDs)
    new_IDs=[]
    for i in range(len(IDs)):
        new_IDs.append('')
        for d in range(6):
            new_IDs[i]+=IDs[i][d]

    urls=[]
    for i in new_IDs:
        urls.append('http://www.uniprot.org/uniprot/'+ i +'.fasta')
    print(urls)
    fastas=[]
    for url in urls:
        fastas.append(WebScrape(url))
    print(fastas)
    with open('temp.txt','w') as file:
        for i in fastas:
            file.write(i)
            file.write('\n')

    FASTAs=[]
    codes=[]
    k=-1
    with open('temp.txt', "r") as file:
        for line in file:
            if line.startswith('>'):
                k+=1
                code=line.strip()
                FASTAs.append('')
                codes.append(code)
            if line.startswith('>')==False:
                FASTAs[k]+= line.strip()
    slices=[]
    slice_length=4
    for f in range(len(FASTAs)):
        slices.append([])
        for i in range(len(FASTAs[f])-(slice_length-1)):
            slices[f].append(FASTAs[f][i:i+slice_length])
    print(slices)
    locations=[]
    logic=[]
    for i in range(len(slices)):
        logic.append([])
        locations.append([])
        for s in range(len(slices[i])):
            logic[i].append([])
            if slices[i][s][0]== 'N':
                logic[i][s].append(True)
            else: logic[i][s].append(False)

            if slices[i][s][1]!= 'P':
                logic[i][s].append(True)
            else: logic[i][s].append(False)

            if slices[i][s][2]== 'S' or slices[i][s][2]== 'T':
                logic[i][s].append(True)
            else: logic[i][s].append(False)

            if slices[i][s][3] != 'P':
                logic[i][s].append(True)
            else: logic[i][s].append(False)

    for i in range(len(slices)):
        for s in range(len(slices[i])):
            if False not in logic[i][s]:
                locations[i].append(s+1)
    print(locations)
    with open('findprotmotif.txt','w') as file:
        for i in range(len(locations)):
            if locations[i]:
                file.write(IDs[i]+'\n')
                for l in locations[i]:
                    file.write(str(l)+'\t')
                file.write('\n')

def reverseTranslate(file_name):
    with open(file_name, 'r') as file:
        prot = file.read().replace('\n', '')
    codons = []
    triplets = []
    with open('codon table.txt', "r") as file:
        codon_table = file.read()
    ct = codon_table.replace("\n", "      ")
    with open('ct.txt', "w") as file:
        file.write(ct)
    with open('ct.txt', "r") as file:
        while True:
            char = file.read(3)
            if not char: break
            triplets.append(char.strip())
            char = file.read(8)
            if not char: break
            codons.append(char.strip())
    print(triplets)
    print(codons)
    unique_codons=[]
    for i in codons:
        if i not in unique_codons:
            unique_codons.append(i)

    triplets_by_codons=[]
    for d in range(len(unique_codons)):
        triplets_by_codons.append([])
    for i in range(len(triplets)):
        for d in range(len(unique_codons)):
            if codons[i]==unique_codons[d]:
                triplets_by_codons[d].append(triplets[i])
    print(unique_codons)
    print(triplets_by_codons)
    print(len(triplets_by_codons[unique_codons.index('M')]))
    print(len(triplets_by_codons[unique_codons.index('A')]))
    possible_RNA_count=3
    for l in prot:
        possible_RNA_count=possible_RNA_count*(len(triplets_by_codons[unique_codons.index(l)]))
        possible_RNA_count=possible_RNA_count%1000000
    print('possible RNA count for this protein = '+str(possible_RNA_count))
def Dictionaries():
    toddlers = dict(age=2,name='Joe')
    toddlers['F'] = 'UUG'
    toddlers['F'] = ['UUG','UUA','UUC']
    toddlers['D']=[]
    D=toddlers['D']
    toddlers['D'].append('ASD')

    print(toddlers)

def consensusSequence(file_name):
    titles=[]
    sequences=[]
    k = -1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                k += 1
                code = line.strip()
                sequences.append('')
                titles.append(code)
            if line.startswith('>') == False:
                sequences[k] += line.strip()
    print(titles)
    print(sequences)
    profile=[]
    for i in range(4):
        profile.append([])
        for d in range(len(sequences[i])):
            profile[i].append(0)
    for i in range(len(sequences)):
        for d in range(len(sequences[i])):
            if sequences[i][d] == 'A':
                profile[0][d]+=1
            if sequences[i][d] == 'C':
                profile[1][d]+=1
            if sequences[i][d] == 'G':
                profile[2][d]+=1
            if sequences[i][d] == 'T':
                profile[3][d]+=1
    print(profile)
    profile_flip=[]
    for i in range(len(sequences[0])):
        profile_flip.append([0,0,0,0])
        for d in range(4):
            profile_flip[i][d]+=profile[d][i]
    print(profile_flip)
    consensus=''
    ACGT='ACGT'
    for i in range(len(profile_flip)):
        x=0
        for d in range(4):
            if profile_flip[i][d] > x:
                x=profile_flip[i][d]
        consensus+=ACGT[profile_flip[i].index(x)]
    print(consensus)
    with open('consensus.txt','w') as file:
        file.write(consensus+'\n')
        for i in range(4):
            file.write(ACGT[i]+':'+'\t')
            for d in range(len(profile[0])):
                file.write(str(profile[i][d])+'\t')
            file.write('\n')

def proteinMass(file_name):
    with open(file_name,'r')as file:
        protein = file.read().strip()
    masstable={}
    with open('monomass.txt', 'r') as file:
        content=file.read().split()
    print(content)
    for i in range(len(content)//2):
        masstable[content[2*i]]=content[(2*i)+1]
    print(masstable)
    mass=0
    for p in protein:
        mass+=float(masstable[p])
        mass=round(mass,5)
    print(mass)

def overlapGraph(file_name,k):
    titles=[]
    sequences=[]
    x = -1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                x += 1
                code = line.strip()
                sequences.append('')
                titles.append(code)
            if line.startswith('>') == False:
                sequences[x] += line.strip()
    print(titles)
    print(sequences)
    result=[]
    for i in range(len(sequences)):
        s=sequences[i]
        for t in range(len(sequences)):
            print(s[-k:])
            print(sequences[t][:(k)])
            if s[-k:] == sequences[t][:(k)]:
                if t!=i:
                    result.append(titles[i].strip('>')+' '+titles[t].strip('>'))
                    print(titles[i])
                    print(titles[t])
    print(result)
    with open('ovGraph.txt', 'w') as file:
        for i in result:
            file.write(str(i)+'\n')

def mendelInheritance(k,m,n): # a population containing k+m+n organisms:
    # #k individuals are homozygous dominant for a factor,
    # m are heterozygous, and
    # n are homozygous recessive
    p= k + m + n
    kk_chance= (k / p)*((k-1)/(p-1))
    mm_chance = (m / p)*((m-1)/(p-1))
    nn_chance= (n / p)*((n-1)/(p-1))
    kn_chance=((k / p)*(n/(p-1)))+((n / p)*(k/(p-1)))
    km_chance=((k / p)*(m/(p-1)))+((m / p)*(k/(p-1)))
    mn_chance=((m / p)*(n/(p-1)))+((n / p)*(m/(p-1)))
    d=kn_chance + (mn_chance*0.5) + (mm_chance*0.75) + km_chance + kk_chance
    r=nn_chance + (mn_chance* 0.5) + (mm_chance*0.25)
    t=d+r
    print(d)
    print(r)
    print(t)

def expectedOffspring(file_name):
    with open(file_name,'r') as file:
        data=file.read().split()
    print(data)
    result=0
    for i in range(6):
        print(result)
        if i <3:
            result+=2*int(data[i])

        if i == 3:
            result+= 1.5*int(data[i])
        if i == 4:
            result+=int(data[i])
    print(result)
def longestCommonSubstring(file_name): #it is slow af but works
    # fasta reader
    titles = []
    sequences = []
    k = -1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                k += 1
                code = line.strip()
                sequences.append('')
                titles.append(code)
            if line.startswith('>') == False:
                sequences[k] += line.strip()
    substrings=[]
    for i in range(len(sequences[0])):
        k=0
        while k+i+1<= len(sequences[0]):
            substrings.append(sequences[0][k:k+i+1])
           # print(substrings)
            k+=1
    common_substrings=[]
    for sub in substrings:
        d=0
        for i in range(len(sequences)):
          #  print(common_substrings)
            if sub in sequences[i]:
                d+=1
        if d==len(sequences):
            common_substrings.append(sub)


    print(common_substrings[-1])

def fastaReader(file_name,output_type=None):
    # fasta reader
    titles = []
    sequences = []

    k = -1
    with open(file_name, "r") as file:
        for line in file:
            if line.startswith('>'):
                k += 1
                code = line.strip()
                sequences.append('')
                titles.append(code)
            if line.startswith('>') == False:
                sequences[k] += line.strip()

    if output_type == 'titles':
        return titles
    else:
        return sequences

def fastaXD(file_name):
    sequnces=fastaReader(file_name,)
    titles=fastaReader(file_name,'titles')
    print(titles[0])
    print(sequnces[0])
fastaXD('rosalind_lcsm.txt')