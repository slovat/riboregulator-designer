#!/usr/bin/env python
# coding: utf-8

# # Designing Riboregulators

#  

# This programme automates the design and analysis of riboregulators. The current design is a riboregulator with a target binding site in the loop. The target RNA binds to the loop and unravels the stem of the hairpin structure, analagous to molecular beacons. Because the target RNA has to compete with the stem for binding, such riboregulators should be highly specific as mismatches are disfavoured as they are unable to outcompete the base pairing in the stem. Note mismatches are disfavoured to a greater extent than for a mismatch binding a linear toehold region as seen with toehold switches as there is no competition for binding in the linear toehold region in toehold switches.  
# 
# This programme uses python wrappers for NUPACK, to allow for easy interfacing
# between Python and calls to the NUPACK core executables. These wrappers were written for python 2, but have been modified for python 3. The wrappers work for all NUPACK executables except: energy, complexes, concentrations, design and distributions

# # Steps 
# 1. Define a specific target sequence and other parameters.
# 2. Divide the sequence into windows of the same length to form the possible triggers and if there are given homolog sequences, check that the trigger overlaps a mismatch site on the homolog to ensure the riboregulator can distinguish between the target and homologs. 
# 3. For each possible trigger, define a riboregulator sequence.
# 4. Evaluate sequence properties using NUPACK and rank them.

# ## 1. Define inputs
# 
# Enter the appropriate values into the cells. 

# #### Enter target sequence here
# 
# Currently, the tool designs riboregulators for RNA targets. Please enter an RNA sequence here using FASTA format (single letter code using A, U, C and G, T is also accepted, but will be converted to U), upper or lower case is fine. The RNA sequence must be surrounded with quotation marks.
# 
# Example entry: "gtactgccaactggatccttcgcgagagcgagtgtgtgtgagcgggatgctgatcgatctaaacgtttagctagctaga" 

# In[1]:


target_seq = "gtactgccaactggatccttcgcgagagcgagtgtgtgtgagcgggatgctgatcgatctaaacgtttagctagctaga" 


# #### Enter homolog sequences here
# 
# The riboregulators are created from the target sequence, but sequences with a small number of mismatches (homologs) will also be able to bind to and activate the riboregutor, albeit at a reduced rate. This option enables you to include homolog sequences and check how well they bind to the riboregulators designed, providing a measure of the specificity of the the riboregulators. Homologs can substantially alter the concentration measured, meaning the true concentration of the target cannot be accurately measured, it is therefore important to check  whether homologs may prevent accurate quantification of the target RNA. Homologs can be identified using BLAST - we recommend searching for similar sequences in all RNA molecules if the target RNA is being quanitified from a total RNA sample (as any RNA species could bind), although for RNAs present in extremely low concentrations, then it may not be worth including them as they won't impact the measured value much. If a specific kit is used to purify RNAs of a specific type, we recommend looking for homologs only of that type of RNA molecule. 
# 
# Sequences must be entered as a list of strings using FASTA format for the sequences. If no homologs are to be entered, give an empty list (i.e. []).  
# 
# Example entry with 2 homolog sequences: 
# homologs = [
#     "gtactgccaactggatccttcgcgagagcg",
#     "atcgatctaaacgtttagctagctaga"
#     ]

# In[2]:


homologs = [
    "gtactgccaactggatccttcgcgagagcg",
    "atcgatctaaacgtttagctagctaga"
    ]


# #### Enter length of trigger here
# 
# The target sequence is divided into trigger sequences by a sliding window. A riboregulator is then generated for each trigger seuqence and the riboregulator is then analysed. The size of the trigger sequence determines how long the region of complementarity between the riboregulator and the target sequence is. The larger the trigger size, the longer the region of complementarity. Note too long a trigger binding site in the loop may result in secondary structures forming in the loop, which can reduce binding kinetics, leading to poor riboregulator activation. Additionally, the binding may be very stable to the point where homologs may be able to bind effectively as the larger number of base pairing can account for the instability introduced by mismatches. However, too short a trigger binding site can result in poor binding kinetics and an inability to unfold the stem. We recommend starting with a length of 15-20 nucleotides.

# In[3]:


trigger_length = 18


# #### Enter simulation temperature here
# 
# Temperature is in degrees Celcius. This determines what temperature should be used for the thermodynamic analysis in NUPACK's calculations.

# In[4]:


simulation_temp = 37


# #### Check inputs

# In[5]:


if type(target_seq) != str:
    raise Exception('Error, target_seq has to be a string')
    

if len(target_seq) == 0: 
    raise Exception("Error, please enter a sequence for target")
else:
    for base in target_seq:
        if base.upper() not in ['A', 'C', 'G', 'U', 'T']:
            raise Exception("""
            Error, please make sure target sequence consists of only A, C, 
            G, U or T. Check that there are not spaces in the sequence.
            """)


# In[6]:


if type(homologs) != list:
    raise Exception("""
    Error, homologs has to be a list, if there are not entries, leave the list empty i.e. []
    """)
    
for i in range(len(homologs)):
    for base in homologs[i]:
        if base.upper() not in ['A', 'C', 'G', 'U', 'T']:
            print(homologs[i])
            raise Exception("""
            Error, please make sure target sequence consists of only
            A, C, G, U or T. Check that there are not spaces in the sequence.
            """)


# In[7]:


if type(trigger_length) != int:
    raise Exception('trigger_length must be an integer')


# In[8]:


if type(simulation_temp) == int or type(simulation_temp) == float:
    pass
else:
    raise Exception('simulation_temp must be an integer or a floating point number')


# ### Import dependencies
# Do not change this section

# In[9]:


import numpy as np
import pandas as pd
from tqdm import tqdm
from nupack import nupack_wrapper as nupack


#  

# ### Functions

# In[10]:


def reversed_complement(sequence):
    """
    Creates the reverse complement sequence of an RNA sequence
    
    Args:
        sequence(str): input RNA strand
    
    Returns: 
        rev_complement(str): reversed complementary sequence of the input RNA 
    """
    nt_pairing = {
        'A': 'U',
        'G': 'C',
        'U': 'A',
        'C': 'G',
        }
    
    sequence_upper = sequence.upper()
    complement = ''
    for nt in sequence_upper:
        complement += nt_pairing[nt]
    rev_complement = complement[::-1]
    
    return rev_complement


def split_sequence(sequence, window_len):
    '''
    Sliding window that splits the sequence into smaller sequences
    
    Args: 
        sequence(str): target RNA strand
        window_len(int): length of sliding window
        
    Returns:
        sequences (list of strs): list of sequence windows from the original sequence 
    '''
    sequences = []
    final_window_start_pos = (len(sequence) - window_len) + 1  
    # +1 as final index value is value specified -1
    for i in range(0, final_window_start_pos):
        sequences.append(sequence[i:window_len + i])

    return sequences


def no_stop(sequence):
    '''Check for stop codons
    
    Assumes that the first nucleotide is also the first position of a codon.
    When no_stop is called in make_riboregulator, it starts checking for 
    stop codons from the start codon.
    
    sequence(str): input RNA strand
    
    returns(Bool): whether the input strand contains a stop codon
    '''
    stop = ['UAA', 'UAG', 'UGA']
    for i in range(0, len(sequence), 3):
        if sequence[i:i + 3] in stop:
            return False

    return True


def make_riboregulators(loop):
    '''Creates list of riboregulators lacking stop codons after the start codon
    
    Make sure to write sequences in RNA form for correct formatting.
    This function can be easily changed for a different design
    
    loop(str): target RNA binds here, it is the reverse complement of the target RNA
    
    returns(lst of strs): list of riboregulators for each trigger
    '''
    riboregs = []
    if no_stop(start_codon + base_3 + linker) == True:
        riboregs.append(stem_5_prime + loop + loop_extension + stem_3_prime + linker)

    return riboregs


# 
# 
#  

# ### Define the sequences of the riboregulator subcomponents
# 
# By explicitly stating all of the subcomponents of the stem it allows for sequences to be easily changed, such as the RBS, without requiring the other parts of the code to be changed (such as where the reading frame starts for the no_stop function and ribo_struct used below).

# #### 3' stem sequence
# Changes are most likely to occur in this part of the stem as this is where the start codon and RBS are.

# In[11]:


RBS = 'AGGAGA' #this is at the top of the stem 
mid_3 = 'G'
internal_loop_3 = 'AAAA'
#internal loop reduces stem stability to increase binding rate 
#and enable unfolding for shorter triggers
start_codon = 'AUG'
base_3 = 'CAG'

stem_3_prime = RBS + mid_3 + internal_loop_3 + start_codon + base_3


# #### 5' stem sequence

# In[12]:


base_5 = reversed_complement(base_3)
start_codon_complement = reversed_complement(start_codon)
internal_loop_5 = 'CCGC'  

mid_5 ='C'
RBS_complement = reversed_complement(RBS)

stem_5_prime = (
    base_5 + start_codon_complement + internal_loop_5 
    + mid_5 + RBS_complement
    )


# #### loop extension
# Currently, this is empty, but it is a placeholder for a sequence between the trigger binding site and RBS, if this ever wants to be added. See notes at the end.

# In[13]:


loop_extension = ""


# #### linker sequence
# The current linker is based on a sequence from doi.org/10.1038/s41589-019-0388-1 supplementary info. This sequence can be optimised to reduce the stability of kinetic intermediates and stabilise the correctly folded structure. The linker sequence will also depend on the reporter used as this has a different sequence.

# In[14]:


linker = 'ACCUGGCGGCAGCGCAAGAAG' 


# ## Step 2: Determine suitable trigger sequences
# Split the target sequence into smaller trigger sequences with a sliding window of length designated by trigger_length.
# 
# Note that by convention, all sequences are written in the 5' to 3' direction (right to left).

# ### Determine all trigger sequences

# In[15]:


if trigger_length > len(target_seq):
    print('Error: trigger length has to be smaller than the target sequence length')


# In[16]:


formatted_seq = target_seq.upper().replace('T', 'U').replace(' ', '')
potential_triggers = split_sequence(formatted_seq, trigger_length) 


# ### Selecting unique trigger sequences
# If the trigger has the same sequence as part of the homolog sequence, then the riboregulator created won't be able to discriminate between the target and homolog RNAs because they will both bind with the same efficacy to the trigger binding site (ignoring binding of the RNAs to regions outside of the trigger binding site), thus K_target/K_homolog is not >> 1. Hence, non-unique trigger sequence (trigger seqeunces that are identical to part of the homolog sequence) are discarded. Note, this is only the first step of removing homologs, once the  K_target/K_homolog has been calculated, it is possible to remove more riboregulators that poorly discriminate against homologs.

# In[17]:


triggers_unique_seq = []
if homologs != []:
    for trigger in potential_triggers:
        if trigger not in str(homologs).upper().replace('T', 'U').replace(' ', ''): 
            triggers_unique_seq.append(trigger)


#  

# ## Step 3: Design a riboregulator for each trigger sequence
# Riboregulators are designed for valid trigger sequence. The loop of the riboregulator is the reverse complement of the trigger RNA.

# In[18]:


if homologs != []:
    riboregulators = [
        make_riboregulators(reversed_complement(trigger))[0] 
        for trigger in triggers_unique_seq
        ]
    df_riboregulators = pd.DataFrame(data=triggers_unique_seq,columns=['triggers'])
elif homologs == []:
    riboregulators = [
        make_riboregulators(reversed_complement(trigger))[0] 
        for trigger in potential_triggers
        ]
    df_riboregulators = pd.DataFrame(data=potential_triggers,columns=['triggers'])
df_riboregulators['riboregulator_sequences'] = riboregulators


#  

#  

# ## Step 4: Evaluate sequence properties using NUPACK 
# NUPACK functions are used for thermodynamic analysis of the riboregulators. 
#  
# 
#  

# ### Specificy the ideal secondary structure of  riboregulator components
# Dot notation is used for RNA secondary structure (parantheses represent base pairing, whilsts dots represent unpaired bases).

# #### Ideal structure (structureless) of the trigger binding site in the loop region:

# In[19]:


binding_struct = trigger_length*'.'
#only if no additional nts are added to the trigger binding site in loop, otherwise adjust


# #### Ideal structure of the riboregulator in the OFF state:

# In[20]:


stem_5_prime_structure = (
    len(base_5)*'(' + len(start_codon_complement)*'(' 
    + len(internal_loop_5)*'.'
    + len(mid_5)*'(' + len(RBS_complement)*'('
    )

stem_3_prime_structure = (
    len(RBS)*')' + len(mid_3)*')'+ len(internal_loop_3)*'.' 
    + len(start_codon)*')' + len(base_3)*')'
    )

loop_structure = binding_struct + len(loop_extension)*'.'  

linker_structure = len(linker)*'.'

ribo_struct = (
    stem_5_prime_structure + loop_structure 
    + stem_3_prime_structure + linker_structure
    )


# #### Ideal structure (structureless) after the trigger binding site once the trigger RNA is bound:

# In[21]:


struct_after_complex = len(loop_extension + stem_5_prime + linker)*'.'


# ### Specificy values for indexing

# In[22]:


#first position after trigger-binding site complex
complex_end = len(stem_5_prime) + trigger_length

#first position of the trigger-binding site complex
complex_start = len(stem_5_prime)


# ### Calculation of the thermodynamic parameters for each design

# In[23]:


num_riboregs = len(df_riboregulators)


# Ensemble Defect: Represents the average number of incorrectly paired nucleotides at equilibrium, evaluated over the ensemble of  a part of the/the whole complex.
# 
# dfull_sensor: Ensemble defect for the full toehold switch sequence and structure.

# In[24]:


dfull_sensor = [
    nupack.defect([df_riboregulators.iloc[i, 1]], ribo_struct, T=simulation_temp) 
    for i in tqdm(range(num_riboregs),desc='dfull_sensor')
    ]


#  

# dactive_sensor: Ensemble defect from the end of the trigger binding site-target RNA complex. A completely single-stranded secondary structure is used for assessing design quality for dactive_sensor.

# In[25]:


dactive_sensor = [
    nupack.defect([df_riboregulators.iloc[i,1][complex_end::]],struct_after_complex, T=simulation_temp)
    for i in tqdm(range(num_riboregs), desc='dactive_sensor')
    ]


#  

# dbinding_site: Ensemble defect for the trigger binding site, specifying a completely single-stranded structure as the ideal structure for the binding site region as an unstructured loop binds the target RNA more efficiently.

# In[26]:


dbinding_site = [
    nupack.defect([df_riboregulators.iloc[i,1][complex_start:complex_end]],binding_struct, T=simulation_temp) 
    for i in tqdm(range(num_riboregs), desc='dbinding_site')
    ]


#  

# dG_RBS_linker: Minimum free energy of the structure between the ribosome binding site and the end of the linker region. In Green et al. 2014, this parameter is the single best predictor of toehold switch function.

# In[27]:


dG_RBS_linker = [
    nupack.mfe([df_riboregulators.iloc[i,1][complex_end+len(loop_extension)::]], T=simulation_temp) 
    for i in tqdm(range(num_riboregs), desc='dG_RBS_linker')
    ]


#  

# dG_target_binding: minumum free energy structure of target bound riboregulator 

# In[28]:


dG_target_binding = [
    nupack.mfe([target_seq,(df_riboregulators.iloc[i][1])], T=simulation_temp) 
    for i in tqdm(range(num_riboregs),desc ='dG_target_binding')
    ]


#  

# dG_homolog_binding: minumum free energy structure of homolog bound riboregulator

# In[29]:


if homologs != []:
    dict_dG_homolog_bindings = {}
    for seq in homologs:
        dict_dG_homolog_bindings[seq] = [
            nupack.mfe([seq,(df_riboregulators.iloc[i][1])], T=simulation_temp) 
            for i in tqdm(range(num_riboregs),desc ='dG_homolog_binding')
            ]
else: 
    dict_dG_homolog_bindings = {}


# In[30]:


df_dG_homolog_bindings = pd.DataFrame([dict_dG_homolog_bindings[key] for key in dict_dG_homolog_bindings]).transpose()


# In[31]:


df_dG_homolog_bindings['min_dG'] = df_dG_homolog_bindings.min(axis=1)
#selects the strongest binding, if the binding is too strong then don't use that riboregulator


#  

# ddG_binding: ΔΔG for target binding. Have to calculate the mfe structures for unbound riboregulator and target RNA and bound riboregulator and subtract the unbound riboregulator and target mfes from the bound mfe.

# In[32]:


dG_closed_riboregulator = [
    nupack.mfe([df_riboregulators.iloc[i,1]], T=simulation_temp) 
    for i in tqdm(range(num_riboregs), desc = 'dG_closed_riboregulator')
    ]

dG_free_target = nupack.mfe([target_seq], T=simulation_temp)

ddG_binding = [
    (dG_target_binding[i]-(dG_closed_riboregulator[i]+dG_free_target)) 
    for i in range(num_riboregs)
    ]


# target_homolog_ratio: ratio of equilibrium constants for target and homolog binding to the riboregulator
# 
# The ratio of equilibrium constants gives an insight into how well a riboregulator discriminates between target and homolog RNAs.

# In[33]:


cal_to_joule = 4.184
temp_kelvin = simulation_temp + 273.15
kcal_to_cal = 1000

if homologs != []:
    K_target = [
        np.exp((dG_target_binding[i]*kcal_to_cal*cal_to_joule)/(-8.31*temp_kelvin)) 
        for i in range(num_riboregs)
    ]
    K_homolog = [
        np.exp((df_dG_homolog_bindings['min_dG'].iloc[i]*kcal_to_cal*cal_to_joule)/(-8.31*temp_kelvin)) 
        for i in range(num_riboregs)
    ]
    target_homolog_ratio = [K_target[i]/K_homolog[i] for i in range(num_riboregs)]
else:
    target_homolog_ratio = []


#  

# Scoring from Ma, D, et al. 2018):
# - Three-parameter fit (R2 = 0.57):
# - Fold change = –71.7 dfull_sensor  – 49.1 dactive_sensor – 22.6 dbinding_site + 54.3
# 
# 
# Note that the R squared values are calculated using 3 or 4 parameters to divide up a data set of 6 toehold switches. Care should be taken with overreliance of these scores as the data has been overfitted (too many parameters are used to distinguish between just 6 toheold switches). A larger dataset is needed to create a better scoring system. A number of different calculations are performed to give insights into the riboregulators, many of which are not included in the current scoring calculation, but may offer insights into which riboregulators are desirable (e.g. K_target/K_homolog is a measure of specificty and shows whether a homolog RNA species may bind well to the riboregulator, the smaller the value, the poorer the specificity is likely to be - note binding kinetics are not included in this).  

# In[34]:


score = (
    54.3 - 71.7*np.array(dfull_sensor) - 49.1*np.array(dactive_sensor)
    - 22.6*np.array(dbinding_site)
    )

target_data = [
    dfull_sensor,
    dactive_sensor,
    dbinding_site,
    dG_RBS_linker,
    ddG_binding,
    dG_target_binding,
    score
    ]

target_columns = [
    'dfull_sensor',
    'dactive_sensor',
    'dbinding_site',
    'dG_RBS_linker (kcal/mol)',
    'ddG_binding (kcal/mol)',
    'dG_target_binding (kcal/mol)',
    'score'
    ]
   
if homologs != []:
    homolog_data = target_data
    homolog_data.insert(-1, list(df_dG_homolog_bindings['min_dG']))
    homolog_data.insert(-1, target_homolog_ratio)
    homolog_columns = target_columns
    homolog_columns.insert(-1, 'dG_homolog_binding (kcal/mol)')
    homolog_columns.insert(-1, 'K_target/K_homolog')
    score_frame = pd.DataFrame(np.stack(homolog_data).T, columns=homolog_columns)
else:
    score_frame = pd.DataFrame(np.stack(target_data).T, columns=target_columns)

sorted_riboregulators = df_riboregulators.join(score_frame).sort_values('score', ascending=False)
pd.set_option('display.max_colwidth', None)
sorted_riboregulators.head()


# ## Save the output to a CSV

# In[37]:


sorted_riboregulators.to_csv("outputs/ranked_designs.csv")


# ## Notes

# - Shorter loops (and hence smaller window) would give a higher specificity as the relative contribution of each nt to the overall binding is greater so it matters more if there is a mismatch. However, a shorter miRNA binding platform would require a shorter stem, which would increase leakage. Could increase trigger binding site by putting some of the trigger binding site in the stem and this won't affect loop size, the RNA would then have to undergo a short strand displacement reaction as well.
# 
# 
# - Could alter loop size e.g add buffer between binding site and RBS. This may be desirable as the double stranded RNA region could interfere with access to the RBS if the dsRNA is too close. 
# 
# 
# - Could change stem lengths too to alter competition between the stem and target binding - alter specificity.
