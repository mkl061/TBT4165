#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import cobra
from cobra import Model, Reaction, Metabolite
import json
import matplotlib as plt
from matplotlib import pyplot





#%% TASK 1

#%% i)
"""
A is imported into the cell with rate R1 and D is exported from the cell with rate R8

"""
#%% ii)
rates=["R"+str(i) for i in range(1,9)]
species=["A", "B", "C", "D", "E", "ATP", "ADP", "Pi"]

df=len(rates)-len(species)

print("The degree of freedom is {} and therefore, S has a zero dimensional solution space".format(df))


#%% iii)
tm = Model('toymodel') #creating the model

S= np.array([[1, -1, 0, 0, -1, 0, 0, 0], #creating the stochiometric matrix
             [0, 1, 1, -1, 0, 0, 0, 0],
             [0, 1, -1, 0, 0, 1, 0, 0],
             [0, 0, 0, 1, -1, 0, 0, -1],
             [0, 0, 0, 0, 1, -1, 0, 0],
             [0, 1, -1, 0, 1, 0, -1, 0], 
             [0, -1, 1, 0, -1, 0, 1, 0],
             [0, -1, 1, 0, -1, 0, 1, 0]])

species = [Metabolite(i) for i in species] #creating a list with all the metabolites

Rs = [Reaction("R"+str(i+1)) for i in range(np.shape(S)[1])] #creating a list with all the reaction rates

# Looping through all reactions Ri's, add their metabolits including the stoch. coeffe. and add them to the model
for i, R in enumerate(Rs):
    # Picking 
    pos=np.argwhere(S[:, i]) #np.argwhere gives all indices in column i (so for Ri) for which the argument is not zero
    pos=pos[:,0] #necessary since we want to have the indices as a list (where they come from a column)
    spec=[species[position] for position in pos ] # picks the metabolites that take part in reaction Ri
    stocs=np.array([S[k,i] for k in pos]) # gets all of the stochiometric coefficients Of Ri
    R.add_metabolites(dict(zip(spec, stocs))) # adding the metabolites to with their stoch coefficients to the reactions
    tm.add_reactions([R])
#%%
tm.reactions.get_by_id("R1").upper_bound=10 # sets the boundary for R1
tm.objective ="R8"# sets reacdtion 8 as the reaction to optimize (max. flux) -> objective function

print("Flux distribution for the model with R1max=10 mmol/(gDWh) and maximizing flux through R8:")
print(tm.optimize().fluxes)
print("\nSo the maximal flux through R8 is {} mmol/(gDWh) and ATP net production is {} mmol/(gDWh).".format(tm.optimize().fluxes["R8"],tm.optimize().fluxes["R7"]))

#%% iv)
"""
When looking at the reactions, one can see that the ATP produced in R2 is directly consumed in R3. From the result above, this implies that NO net ATP production occurs since
the only other ATP producing reaction besides R2 is R5, which has a flux of 0. So even if one optimizes for ATP, I think
one shouldn't observe a flux for R7.
"""

tm.reactions.get_by_id("R1").upper_bound=10 # sets the boundary for R1
tm.objective ="R7"# sets reacdtion 8 as the reaction to optimize (max. flux) -> objective function
print(tm.optimize().fluxes)
print("\nApparently, the result doesn't change when optimizing for R7 and the ATP production, so the flux through R7 still remains {} mmol/(gDWh)".format(tm.optimize().fluxes["R7"]))


#%% v)
"""
I'm not sure about this one. My thinking: By doing this, we "faktisk" add a pseudo-species that doesn't exist in our reaction scheme.
Because we add only a 1 in column 6, this means that this pseudo-species gets produced only in reaction 6 and more important,
is NOWHERE consumed -> in order to follow the mass balance and steady state, the flux through R6 therefore MUST be zero,
otherwise this pseudo-species would accumulate and d[pseudo]/dt != 0.

"""



#%% TASK 2

# i)
# Reading the model and printing how many genes, reactions and metabolic species we have
ecore = cobra.io.read_sbml_model("ecoli_core_model/ecoli_core_model.xml")
print("In our E. coli core model we have {} metabolic species, {} reactions and {} genes!".format(len(ecore.metabolites),len(ecore.reactions),len(ecore.genes)))

# Metabolic subsystems
#print(ecore.groups())

print("Subsystems:\n")
subsystems=[]
for i in ecore.reactions:
    if i.subsystem not in subsystems:
        subsystems.append(i.subsystem)
        print(i.subsystem)
print("\nUnfortunately I couldn't manage with the groups function.")


#%% ii)
ecore.reactions.get_by_id("EX_glc__D_e").lower_bound = -18.5 # set the glucose uptake lower bound to -18.5
ecore.reactions.get_by_id("EX_o2_e").lower_bound = -1000 # set the oxygen uptake lower bound to -1000
print(ecore.summary())
print("So the max. spec. growth rate is {:.2f} (since this is the objective, biomass), the flux of glucose etc. can be seen in the table.".format(ecore.optimize().objective_value))


#%% iii)
ecore.reactions.get_by_id("EX_o2_e").lower_bound = 0 # anaerobic conditions
print(ecore.summary())
"""
By comparing the summaries, one can see that now, under anaerobic conditions, E. coli starts secreting acetat, ethanol, formiat and H+ instead of CO2 and H2O and H+ like under aerobic conditions.
This makes sense since E. coli is facultative anaerob, so under anaerobic conditions it will start fermentation, leading to the excretions described above.

"""

#%% iv
ecore.reactions.get_by_id("EX_o2_e").lower_bound = -1000 # set the oxygen uptake lower bound to -1000
aerob_dic = ecore.optimize().fluxes.to_dict() # writes all fluxes to the dictionary, aerobic conditions
ecore.reactions.get_by_id("EX_o2_e").lower_bound = 0 # set the oxygen uptake lower bound to -1000
anaerob_dic = ecore.optimize().fluxes.to_dict()

aerob_object= json.dumps(aerob_dic, indent=4)
anaerob_object = json.dumps(anaerob_dic, indent=4)

# Writing to json
for obje, label in zip([aerob_object, anaerob_object],["aerob", "anaerob"]):
    with open(label+".json", "w") as outfile:
        outfile.write(obje)

"""
Unfortunately, I have a problem with the Escher visualization: I can upload the .json file, but then the screen remains white, 
so no visualization happens. 
"""

#%%

substrate_list=["EX_ac_e", "EX_acald_e", "EX_akg_e", "EX_etoh_e", "EX_fru_e", "EX_fum_e", "EX_glc__D_e", "EX_gln__L_e", "EX_glu__L_e", "EX_lac__D_e", "EX_mal__L_e", "EX_pyr_e", "EX_succ_e"]

sub_dic={}
substratesnotinmodel=[]

#%%
for substrate in substrate_list:
    if substrate in aerob_dic.keys(): #apparently there are substrates which are not in the model. So we check if they're in the flux dictionary from above
        
        # setting the substrate uptake rate to 20 for substrate i
        ecore.reactions.get_by_id(substrate).lower_bound = -20 # substrate rate = 20
        
        sub_dic[substrate]={}
        
        # adding the optimized growth rate for each condition
        for oxygen, cond in zip([-1000,0],["aerob","anaerob"]):
            ecore.reactions.get_by_id("EX_o2_e").lower_bound = oxygen
            sub_dic[substrate]|={cond: ecore.optimize().objective_value} # update the dictionary
            
        ecore.reactions.get_by_id(substrate).lower_bound = 0 # setting the substrate uptake back again for substrate i

    # Keeps track for which substrates 
    else:
        substratesnotinmodel.append(substrate)
        

# Result
aero_l=[]
anaero_l=[]
sub_l=list(sub_dic.keys())

for substrate in substrate_list:
    if substrate in sub_dic:
        for cond in sub_dic[substrate].keys():
            if cond=="aerob":
                aero_l.append(sub_dic[substrate][cond])
            if cond == "anaerob":
                anaero_l.append(sub_dic[substrate][cond])
                
            print("Optimal growth rate for {} under {} conditions is: \t {:.2f}".format(substrate, cond, sub_dic[substrate][cond]))

    else:
        print("Unfortunately, {} was not found in the E. Coli model as a substrate :-(".format(substrate))
# Well this looks horrible -> plot it!


# Plotting optimal growth rate vs. substrate, categorical according to conditions
fig, ax = pyplot.subplots()
ax.scatter(sub_l, aero_l, label="aerob")
ax.scatter(sub_l, anaero_l, label="anaerob")
ax.legend()
pyplot.title("Optimal growth rate for different substrates")
pyplot.ylabel("Optimal growth rate")
pyplot.xlabel("Substrate")
pyplot.show()

