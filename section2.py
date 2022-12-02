#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyopenms import *
seq = AASequence.fromString("DFPIANGER") 
prefix = seq.getPrefix(4)                  
suffix = seq.getSuffix(5)                   
concat = seq + seq                         
print("Sequence:", seq)
print("Prefix:", prefix)
print("Suffix:", suffix)
print("Concatenated:", concat)


# In[2]:


mfull = seq.getMonoWeight()                           
mprecursor = seq.getMonoWeight(Residue.ResidueType.Full, 2)   
mz = seq.getMonoWeight(Residue.ResidueType.Full, 2) 
mz = seq.getMZ(2) 
print("Monoisotopic mass of peptide [M] is", mfull)
print("Monoisotopic mass of peptide precursor [M+2H]2+ is", mprecursor)
print("Monoisotopic m/z of [M+2H]2+ is", mz)


# In[3]:


seq = AASequence.fromString("DFPIANGER")
print("The peptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())


# In[4]:


seq = AASequence.fromString("DFPIANGER")
seq_formula = seq.getFormula()
print("Peptide", seq, "has molecular formula", seq_formula)


# In[5]:


suffix = seq.getSuffix(3)                              
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2)       
suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0            
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0            
print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
print("y3 molecular formula:", y3_formula)


# In[6]:


seq = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
print(seq.toUnmodifiedString())
print(seq.toString())
print(seq.toUniModString())
print(seq.toBracketString())
print(seq.toBracketString(False))
print("="*35)
print(AASequence.fromString("DFPIAM(UniMod:35)GER"))
print(AASequence.fromString("DFPIAM[+16]GER"))
print(AASequence.fromString("DFPIAM[+15.99]GER"))
print(AASequence.fromString("DFPIAM[147]GER"))
print(AASequence.fromString("DFPIAM[147.035405]GER"))


# In[7]:


bsa = FASTAEntry() 
bsa.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
bsa.description = "BSA Bovine Albumin (partial sequence)"
bsa.identifier = "BSA"
alb = FASTAEntry()
alb.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE"
alb.description = "ALB Human Albumin (partial sequence)"
alb.identifier = "ALB"
entries = [bsa, alb]
f = FASTAFile()
f.store("example.fasta", entries)


# In[8]:


entries = []
f = FASTAFile()
f.load("example.fasta", entries)
print( len(entries) )
for e in entries:
    print (e.identifier, e.sequence)


# In[ ]:




