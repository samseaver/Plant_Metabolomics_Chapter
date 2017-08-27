# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:17:41 2017
@author: Nadine Toepfer
"""
### 0. Run this first: import modules and define functions  
##
#

from __future__ import print_function
import cobra
import pandas as pd
import numpy as np
from cobra import Reaction
import cobra.core.model
import matplotlib.pyplot as plt
import collections
    
def format_long_string(string, max_length):
    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string    
    
def _process_flux_dataframe(flux_dataframe, threshold, floatfmt):
    flux_dataframe = flux_dataframe[flux_dataframe.flux.abs() > threshold].copy()
    flux_dataframe['is_input'] = flux_dataframe.flux >= 0
    flux_dataframe.flux = \
        flux_dataframe.flux.abs().astype('float')
    return flux_dataframe

### 1. read generic tomato model
##
#
print("Reading the tomato model...")
model=cobra.io.read_sbml_model("Tomato_PlantSEED_Model.sbml")
# set solver, this can be changed to e.g. GUROBI or glpk
model.solver = "cplex"


### 2. print basic model properties
##
#
print("number of reactions: %i" % (len(model.reactions)))
print("number of metabolites: %i" % (len(model.metabolites)))
print("number of genes: %i \n" % (len(model.genes)))

### 3. print exchange reactions
##
#
print("Exchange Reactions: \n")
print("Reaction ID","\t", "Reaction Name","\t","Lower Bound","\t","Upper Bound")
print("---------------------------------------------------------------")
for x in model.exchanges:
    print(x.id,"\t", x.name,"\t",x.lower_bound,"\t",x.upper_bound)
print("\n")

### 4. print biomass reaction
##
#   
print("Biomass composition: \n") 
bio=model.reactions.get_by_id("bio1")   
coeffs=bio.get_coefficients(bio.reactants)
print("Metabolite name","\t","coefficient")
print("-------------------------------------")
for i in range(0,len(coeffs)):
    print(bio.reactants[i].name,"\t",coeffs[i])
print("\n")
    
### 5. run FBA
##
#  
print("Run test FBA: \n")    
model.optimize(new_objective=bio)
sol = model.optimize()
model.summary()
print("\n")

### 6.change uptakes and run FBA
##
#
print("Change uptake rates: \n")  
print(model.exchanges[6].name) #  urea
model.exchanges[6].lower_bound=0
print("New lower exchange boundary:",model.exchanges[6].lower_bound)

print(model.exchanges[11].name) # sucrose
model.exchanges[11].upper_bound=0 # not out
print("New lower exchange boundary:",model.exchanges[11].upper_bound)

print(model.exchanges[14].name) # light
model.exchanges[14].lower_bound=0
print("New lower exchange boundary:",model.exchanges[14].lower_bound)

print(model.exchanges[16].name) # biomass
model.exchanges[16].lower_bound=0
print("New lower exchange boundary:",model.exchanges[16].lower_bound)
print("\n")

# these values will serve as default values for comparision
print("Run FBA with new boundaries: \n")
def_sol = model.optimize()
def_f=def_sol.objective_value
model.summary()
print("\n")

# run parsimonious FBA (pFBA)
print("Run pFBA: \n")
def_pfba_sol=cobra.flux_analysis.pfba(model)
def_pfba_f = cobra.flux_analysis.pfba(model).f
print("pFBA solution: %.2f" % def_pfba_f)
print("\n")

# extract active subnetworks and save reaction and metabolite lists
print("\n")
S=cobra.util.array.create_stoichiometric_matrix(model, array_type="dense", dtype=None) # get S

default_active_rxns=[i for i, x in enumerate(def_pfba_sol.fluxes!=0) if x] # find index of active reactions
default_active_rxn_names=[model.reactions[i].name for i in default_active_rxns] 
S_reduced_temp=S[:,default_active_rxns]

default_active_mets=[i for i, x in enumerate(~(S_reduced_temp==0).all(1)) if x]  # remove inactive reactions and metabolites from S
default_active_met_names=[model.metabolites[i].name for i in default_active_mets] 

print("Active network size: %i metabolites %i reactions. \n" % (len(default_active_mets), len(default_active_rxns)))

# get default exchange fluxes
def_ex_fluxes=def_sol.fluxes
ex_fluxes=def_ex_fluxes[['EX_cpd00001_e0',
              'EX_cpd00007_e0',
              'EX_cpd00009_e0', 
              'EX_cpd00011_e0',
              'EX_cpd00013_e0',
              'EX_cpd00048_e0',
              'EX_cpd00067_e0', 
              'EX_cpd00073_e0',
              'EX_cpd00076_e0',
              'EX_cpd00099_e0',
              'EX_cpd00204_e0',
              'EX_cpd00205_e0',
              'EX_cpd00209_e0',
              'EX_cpd11632_e0',
              'EX_cpd02701_c0', 
              'EX_cpd11416_c0']]
          
### 7. read the matched metabolites and data
##
#   
data = []
with open("Matched_Metabolomics_Data.txt") as file:
    for line in file:
        list = [element.strip() for element in line.split("\t")]
        data.append(list)        

### 8. perform data integration
##
#

context_model_properties=collections.OrderedDict()

threshold=1E-8 # threashold for calling a flux zero
floatfmt='.3g' # format method for floats, passed to tabulate. Default is '.3g'

# add empty col, this is just for the plot later             
zero_col=ex_fluxes.copy()
zero_col[0:len(ex_fluxes)]=0    
ex_fluxes=pd.concat([ex_fluxes,zero_col],axis=1)
  
for j in range(3,8): # flesh only

    context_model = model.copy() # make a copy of the original model
    print("integrating metabolite data for conditon: %i \n" % (j-2))    
    measured_metabolite_list=[]    
    
    for i in range(0,len(data)):    
        if  int(data[i][j]):
            measured_metabolite_list.append(data[i][2])
         
    for metabolite in measured_metabolite_list:        
        context_model_tmp=context_model.copy()           
        # get list of a metabolites in different compartments that match the substr       
        metabolites_matched = [met_in_compartment for met_in_compartment in context_model_tmp.metabolites
                                if metabolite in str(met_in_compartment)]      
        
        # if metabolites are in context_model_tmp
        if (len(metabolites_matched)>0):            
            # add sink reaction for each compartment-specific metabolite and add reaction to list
            sink_rxn_list=[]
            
            for met in metabolites_matched:
                reaction = Reaction(str(met) +"_sink")
                reaction.add_metabolites({met: -1.0})
                sink_rxn_list.append(reaction)
                context_model_tmp.add_reaction(reaction)
                
            # generate overall constraint the for overall sink flux through metabolite       
            for i in range(0,len(sink_rxn_list)):
                constraint_expression=eval("context_model_tmp.reactions."+ sink_rxn_list[i].id +".flux_expression")
                if (i==0):
                    constraints_all= constraint_expression
                else:
                    constraints_all+= constraint_expression                   
                    
            sink_flux = context_model_tmp.problem.Constraint(constraints_all,lb=0.01,ub=len(sink_rxn_list)*1000)        
            context_model_tmp.add_cons_vars(sink_flux)
            solution = context_model_tmp.optimize()
           
            if solution.status=="optimal":
                context_model=context_model_tmp.copy()
  
    solution = context_model.optimize()     
   
    # save exchange
    exchange_fluxes = {}
    for rxn in context_model.exchanges:
        for met, stoich in dict.iteritems(rxn.metabolites):
            exchange_fluxes[met] = {'id': format_long_string(met.id, 15),'flux': stoich * rxn.flux}

    exchange_fluxes = pd.DataFrame(exchange_fluxes).T 
    exchange_fluxes = _process_flux_dataframe(exchange_fluxes, threshold, floatfmt)
    
    # FBA
    fba_solution = context_model.optimize()
    print("FBA solution: %.2f" % fba_solution.f)
    
    # pFBA
    pfba_solution = cobra.flux_analysis.pfba(context_model)
    print("pFBA solution: %.2f" % pfba_solution.f)   
   
    # extract active subnetworks and save reaction and metaboolite lists
    S=cobra.util.array.create_stoichiometric_matrix(context_model, array_type="dense", dtype=None)
    active_rxns=[i for i, x in enumerate(pfba_solution.fluxes!=0) if x] # check  # find index of active reactions
    active_rxn_names=[context_model.reactions[i].name for i in active_rxns] 
    S_reduced_temp=S[:,active_rxns]
    active_mets=[i for i, x in enumerate(~(S_reduced_temp==0).all(1)) if x] # check # remove inactive reactions and metabolites from stoichiometric matrix    
    active_met_names=[context_model.metabolites[i].name for i in active_mets] 
    print("Active network size: %i metabolites %i reactions. \n" % (len(active_mets), len(active_rxns)))

    context_model_properties[context_model]={'model': context_model,
                                            'objective value':solution.f,
                                            'pFBA value': pfba_solution.f,
                                            'solution fluxes': solution.fluxes,
                                            'pFBA fluxes': pfba_solution.fluxes,
                                            'Nr. active rxn':len(active_rxns),
                                            'Nr. active mets': len(active_mets),
                                            'active rxns': active_rxns,
                                            'active mets':active_mets,
                                            'active rxn names':active_rxn_names,
                                            'active met names':active_met_names,
                                            'exchange fluxes': exchange_fluxes}
                                            
### 9. analyse results
##
#
for context_model in context_model_properties:  
  
    cm_fluxes=context_model_properties[context_model]['solution fluxes'][0:len(def_pfba_sol.fluxes)] # newly added sink fluxes are at the end
    df=cm_fluxes[['EX_cpd00001_e0',
                  'EX_cpd00007_e0',
                  'EX_cpd00009_e0', 
                  'EX_cpd00011_e0',
                  'EX_cpd00013_e0',
                  'EX_cpd00048_e0',
                  'EX_cpd00067_e0', 
                  'EX_cpd00073_e0',
                  'EX_cpd00076_e0',
                  'EX_cpd00099_e0',
                  'EX_cpd00204_e0',
                  'EX_cpd00205_e0',
                  'EX_cpd00209_e0',
                  'EX_cpd11632_e0',
                  'EX_cpd02701_c0', 
                  'EX_cpd11416_c0']]  
                  
    ex_fluxes=pd.concat([ex_fluxes,df],axis=1)

# drop all exchange fluxes that are always zero
ex_fluxes = ex_fluxes[(ex_fluxes.T != 0).any()] 

# this is manuallly asigned from the output
names=['H2O', # cpd00001
       'O2', # cpd00007
       'Phosphate', #cpd00009
       'CO2_e0', #cpd00011
       'NH3', #cpd00013
       'Sulfate', # cpd00048
       'H_plus', # cpd00067
       'Sucrose', # cpd00076
       'Biomass c']  # cpd11416
   
exchange_names=[dstr for dstr in ex_fluxes.index] 

## plot heatmap for exchange fluxes
#plt.show() # avoids error message when run in console
fig, img = plt.subplots(1,1)
fg=img.imshow(ex_fluxes, cmap="RdBu", interpolation="nearest",vmin=-1000, vmax=1000)
img.set_xticklabels(['','Def','', 'IG','MG','Br','Or','Rd'])
img.set_yticks(np.arange(len(exchange_names)))
img.set_yticklabels(names)
plt.title('Exchange Fluxes')
plt.colorbar(fg)
fig.savefig("Exchange_fluxes.svg",dpi=600)

# get values for max biomass
biomass=[def_f]+[0]+[context_model_properties[cm]['objective value'] for cm in context_model_properties][0:5]
x1=[0.5,1.5,2.5,3.5,4.5,5.5,6.5]
x2=[0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9]

# plot values for max biomass
fig=plt.figure()
plt.bar(x1,biomass,1/1.5,color=['#0a335c','#003366','#0a5c0a','#8fb814','#e8c330','#e68019','#b81414'],linewidth=0) 
ax = plt.gca()  
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.xticks(x2, ['Def','','IG','MG','Br','Or','Rd'])
plt.tick_params(axis='x',which='both',bottom='off') 
plt.title('Maximum biomass')
plt.xlabel('Ripening Stage')
plt.ylabel('Objective value')
fig.savefig("Biomass.svg",dpi=600)


# get values for pFBA
pFBA=[def_pfba_f]+[0]+[context_model_properties[cm]['pFBA value'] for cm in context_model_properties][0:5]

# plot values for pFBA
fig=plt.figure()
plt.bar(x1,pFBA,1/1.5,color=['#0a335c','#003366','#0a5c0a','#8fb814','#e8c330','#e68019','#b81414'],linewidth=0) 
ax = plt.gca()  
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.xticks(x2, ['Def','','IG','MG','Br','Or','Rd'])
plt.tick_params(axis='x',which='both',bottom='off') 
plt.title('pFBA solutions')
plt.xlabel('Ripening Stage')
plt.ylabel('Sum of fluxes')
fig.savefig("pFBA_solutions.svg",dpi=600)

############################################################################

# what causes the change between time 1 and 2?

# get model
cm_1=context_model_properties[context_model_properties.keys()[0]]
cm_2=context_model_properties[context_model_properties.keys()[1]]

## 50 most changing reactions 
cm_1_fluxes=cm_1['solution fluxes'][0:len(def_pfba_sol.fluxes)] # newly added sink fluxes are at the end
cm_2_fluxes=cm_2['solution fluxes'][0:len(def_pfba_sol.fluxes)] # newly added sink fluxes are at the end

cm_1_all_rxns=[r for r in cm_1['model'].reactions]
cm_2_all_rxns=[r for r in cm_2['model'].reactions]

cm_1_all_rxns_names=[r.name for r in cm_1['model'].reactions]
cm_2_all_rxns_names=[r.name for r in cm_2['model'].reactions]

diff=abs(cm_1_fluxes-cm_2_fluxes)

# indexes of ordered list of differences between conditions 1 and 2
order=[i[0] for i in sorted(enumerate(abs(cm_1_fluxes-cm_2_fluxes)), reverse=True, key=lambda x:x[1])]

# print most changing reactions
print("print most changing reactions:")
print("\n")
print("Reaction ID","\t", "Reaction Name","\t","Flux 1","\t","Flux 2")
print("---------------------------------------------------------------")
for i in range(0,50):
    print(cm_1_all_rxns[order[i]],"\t", cm_1_all_rxns_names[order[i]],"\t",cm_1_fluxes[order][i],"\t",cm_2_fluxes[order][i], diff[order[i]])
print("\n")



