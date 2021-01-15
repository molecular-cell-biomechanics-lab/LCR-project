#!/usr/bin/env python
# coding: utf-8

# In[11]:


#! /usr/bin/python

from __future__ import division
from collections import defaultdict
import matplotlib.pyplot as plt
import sys
import numpy
import math
import re
import operator
from Bio import SeqIO
import numpy as np
import json
from matplotlib.pyplot import hist, xticks, boxplot
from matplotlib import image


title = ""
# All proteins are organized into a dict of dicts with disprot ID as keys
# Inner dicts contain pos_LCR_length, pos_LCR_charge_count, same for negMK,M.M
# Thus four items total for each ID
# Dictionary is because separate disprot regions can come from same ID
prot_dict = {}

fgnups = open("NupsDisRegionsWithName.txt")
for record in fgnups:
	if record[0] == '>':
		ID = record[1:10]
        #print(ID)
		seqlen = 0
	else:
		seqlen += len(record)
		prot_dict[ID] = {'pos_LCR_length':0, "neg_LCR_length":0, "seqlen": seqlen}

IDnum = -1
fgnups = open("NupsDisRegionsWithName.txt")
for record in fgnups:
	if record[0] == '>':
		ID = record[1:10]
		IDnum += 1
	else:
		max_pos_len = prot_dict[ID]['pos_LCR_length']
		max_neg_len = prot_dict[ID]['neg_LCR_length']
		curr_pos, curr_neg = False, False
		curr_start, curr_end = 0, 0
		curr_charge_count = 0
		sequence = record
		for i in range(len(sequence)):
			if sequence[i] == 'D' or sequence[i] == 'E':
				if curr_pos:
					curr_length = curr_end - curr_start + 1
					if curr_length > max_pos_len and curr_charge_count > 2 and curr_length > 1:
						max_pos_len = curr_length
						max_pos_start_position = curr_start
						max_pos_end_position = curr_end
						max_pos_charge = curr_charge_count
					curr_pos = False
				if not curr_neg:
					curr_neg = True
					curr_start = i
					curr_end = i
					curr_charge_count = 1
					curr_length = 0
				else:
					curr_end = i
					curr_charge_count+=1
			elif sequence[i] == 'K' or sequence[i] == 'R':
				if curr_neg:
					curr_length = curr_end - curr_start + 1
					if curr_length > max_neg_len and curr_charge_count > 2 and curr_length > 1:
						max_neg_len = curr_length
						max_neg_start_position = curr_start
						max_neg_end_position = curr_end
					curr_neg = False
				if not curr_pos:
					curr_pos = True
					curr_start = i
					curr_end = i
					curr_charge_count = 1
					curr_length = 0
				else:
					curr_end = i
					curr_charge_count+=1
		if curr_pos:
			curr_length = curr_end - curr_start + 1
			if curr_length > max_pos_len and curr_charge_count > 2 and curr_length > 1:
				max_pos_len = curr_length
				max_pos_start_position = curr_start
				max_pos_end_position = curr_end
			curr_pos = False
		elif curr_neg:
			curr_length = curr_end - curr_start + 1
			if curr_length > max_neg_len and curr_charge_count > 2 and curr_length > 1:
				max_neg_len = curr_length
				max_neg_start_position = curr_start
				max_neg_end_position = curr_end
				max_neg_charge = curr_charge_count
			curr_neg = False
		prot_dict[ID]['pos_LCR_length'] = max_pos_len
		prot_dict[ID]['neg_LCR_length'] = max_neg_len
        
neg_LCR_lengths = np.array([])
norm_neg_LCR_lengths = np.array([])

pos_LCR_lengths = np.array([])
norm_pos_LCR_lengths = np.array([])

for ID in prot_dict.keys():
	if prot_dict[ID]['neg_LCR_length'] > 0: 
		neg_LCR_lengths = np.append(neg_LCR_lengths, prot_dict[ID]['neg_LCR_length'])
		norm_neg_LCR_lengths = np.append(norm_neg_LCR_lengths, prot_dict[ID]['neg_LCR_length']/prot_dict[ID]['seqlen'])
	if prot_dict[ID]['pos_LCR_length'] > 0:
		pos_LCR_lengths = np.append(pos_LCR_lengths, prot_dict[ID]['pos_LCR_length'])
		norm_pos_LCR_lengths = np.append(norm_pos_LCR_lengths, prot_dict[ID]['pos_LCR_length']/prot_dict[ID]['seqlen'])


fig = plt.figure(figsize = (12,6))
axes = fig.add_subplot(1,2,1)

axes.scatter(np.array(pos_LCR_lengths), pos_ND, label = 'largest positive LCRs')
axes.set_xlim(xmin = 0, xmax= 700)
axes.set_ylim(ymin = 0, ymax = 1.0)
axes.scatter(np.array(neg_LCR_lengths), neg_ND, facecolor = 'orange', label = 'largest negative LCRs')
axes.set_xlim(xmin = 0, xmax = 700)
axes.set_ylim(ymin = 0, ymax = 1.0)
axes.set_xlabel('largest LCR length', fontsize = 14)
axes.set_ylabel('largest LCR charge content (Number Density)', fontsize = 12)
axes.set_title(title)
axes.legend()

# Scatter plot with normalized LCR lengths
axes = fig.add_subplot(1,2,2)
axes.scatter(np.array(norm_pos_LCR_lengths), pos_ND, facecolor = 'blue',label = 'Positive LCRs')
axes.scatter(np.array(norm_neg_LCR_lengths), neg_ND, facecolor = 'orange', label = 'Negative LCRs')
#axes.set_xlim(xmin = 0, xmax = 1)
#axes.set_ylim(ymin = 0, ymax = 1.0)
axes.set_xlabel('Normalized LCR length', fontsize = 14)
axes.set_ylabel('LCR charge content (Number Density)', fontsize = 12)

axes.legend()
axes.set_title(title)

plt.savefig("Length_Density_Scatterplot.png")
plt.show()


# In[ ]:




