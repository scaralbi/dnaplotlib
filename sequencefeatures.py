#!/usr/bin/env python
"""
	Visualize sequence features
"""

import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

# Required for drawing shapes
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches


import matplotlib as mpl
from matplotlib import cm
from colorspacious import cspace_converter
from collections import OrderedDict



cmaps = OrderedDict()

# Color maps
col_map = {}
col_map['red']     = (0.95, 0.30, 0.25)
col_map['green']   = (0.38, 0.82, 0.32)
col_map['blue']    = (0.38, 0.65, 0.87)
col_map['orange']  = (1.00, 0.75, 0.17)
col_map['purple']  = (0.55, 0.35, 0.64)
col_map['yellow']  = (0.98, 0.9, 0.55)
col_map['black']  = (0, 0, 0)
col_map['gray']  = (0.41, 0.41, 0.41)


Cmap = cm.get_cmap('Pastel2', 2)



##############################################################################
# HELPER FUNCTIONS
# Could be integrated into DNAplotlib in the future. You can create renderer
# functions for each type of annoation and these are then easily used by
# DNAplotlib. I created a very simple one for what we discussed - these can
# be customised further by the opts dictionary in the later section to capture
# metadata such as strength, etc.
##############################################################################

# Function for writing a sequence to axis
def draw_sequence (ax, start_pos, y_offset, seq):
	# For each letter position in center of range
	for c_idx, c in enumerate(seq):
		ax.annotate(c,(start_pos+c_idx+0.5, y_offset), ha='center', va='bottom', fontsize=6)

# Custom renderer for promoter -35 to -10 site (this is used by DNAplotlib to render the sites in a design)
def promoter_region (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	# Default parameters - these can be added to, but we usually use this style (probably should simplify in future)
	y_offset = 0.0
	color_35 = (0.5,0.5,0.5)
	color_10 = (0.5,0.5,0.5)
	color_connector = (0,0,0)
	linewidth_connector = 1.0
	len_35 = 4
	len_10 = 2
	y_extent = 2.0
	# Update default parameters if provided
	if opts != None:
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'y_offset' in list(opts.keys()):
			y_offset = opts['y_offset']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'color_35' in list(opts.keys()):
			color_35 = opts['color_35']
		if 'color_10' in list(opts.keys()):
			color_10 = opts['color_10']
		if 'color_connector' in list(opts.keys()):
			color_connector = opts['color_connector']
		if 'linewidth_connector' in list(opts.keys()):
			linewidth_connector = opts['linewidth_connector']
		if 'len_35' in list(opts.keys()):
			len_35 = opts['len_35']
		if 'len_10' in list(opts.keys()):
			len_10 = opts['len_10']
	# Check direction (we don't use at moment)
	fwd = True
	if start > end:
		fwd = False
	# Draw the -35 site (from start to start + length of -35 site)
	p35 = Polygon([(start, y_offset),
		           (start, y_offset+y_extent),
		           (start+len_35,y_offset+y_extent),
		           (start+len_35,y_offset)],
		            edgecolor=(0,0,0), facecolor=color_35, linewidth=linewidth, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(p35)
	# Draw the -10 site (from end-length of -10 site to end)
	p10 = Polygon([(end-len_10, y_offset),
		           (end-len_10, y_offset+y_extent),
		           (end,y_offset+y_extent),
		           (end,y_offset)],
		            edgecolor=(0,0,0), facecolor=color_10, linewidth=linewidth, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(p10)
	l1 = Line2D([start+len_35, end-len_10],
                [y_offset+(y_extent/2.0), y_offset+(y_extent/2.0)], linewidth=linewidth_connector,
                color=color_connector, zorder=10)
	ax.add_line(l1)

	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			dpl.write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			dpl.write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer


# Custom renderer for promoter -35 to -10 site (this is used by DNAplotlib to render the sites in a design)
def operator_region (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	# Default parameters - these can be added to, but we usually use this style (probably should simplify in future)
	y_offset = 0
	color_O1 = col_map['gray']
	color_O2 = col_map['gray']
	color_connector = (0,0,0)
	linewidth_connector = 1
	len_O1 = 17
	len_O2 = 17
	y_extent = 1
	# Update default parameters if provided
	if opts != None:
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'y_offset' in list(opts.keys()):
			y_offset = opts['y_offset']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'color_O1' in list(opts.keys()):
			color_O2 = opts['color_O1']
		if 'color_O2' in list(opts.keys()):
			color_O2 = opts['color_O2']
		if 'color_connector' in list(opts.keys()):
			color_connector = opts['color_connector']
		if 'linewidth_connector' in list(opts.keys()):
			linewidth_connector = opts['linewidth_connector']
		if 'len_O1' in list(opts.keys()):
			len_O1 = opts['len_O1']
		if 'len_O2' in list(opts.keys()):
			len_O2 = opts['len_O2']
	# Check direction (we don't use at moment)
	fwd = True
	if start > end:
		fwd = False
	# Draw the -35 site (from start to start + length of -35 site)
	pO1 = Polygon([(start, y_offset),
		           (start, y_offset+y_extent),
		           (start+len_O1,y_offset+y_extent),
		           (start+len_O1,y_offset)],
		            edgecolor=(1,1,1), facecolor=color_O1, linewidth=0.2, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(pO1)
	# Draw the -10 site (from end-length of -10 site to end)
	pO2 = Polygon([(end-len_O2, y_offset),
		           (end-len_O2, y_offset+y_extent),
		           (end,y_offset+y_extent),
		           (end,y_offset)],
		            edgecolor=(1,1,1), facecolor=color_O2, linewidth=0.2, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(pO2)
	lO1 = Line2D([start+len_O1, end-len_O2],
                [y_offset+(y_extent/2.0), y_offset+(y_extent/2.0)], linewidth=linewidth_connector,
                color=color_connector, zorder=10)
	ax.add_line(lO1)

	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			dpl.write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			dpl.write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer
	return start, end

# Custom renderer for promoter -35 to -10 site (this is used by DNAplotlib to render the sites in a design)
def palindrome_region (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	# Default parameters - these can be added to, but we usually use this style (probably should simplify in future)
	y_offset = 0
	color_1 = Cmap(1)
	color_2 = Cmap(1)
	color_connector = (0,0,0)
	linewidth_connector = 1
	len_1 = 10
	len_2 = 10
	y_extent = 1
	# Update default parameters if provided
	if opts != None:
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'y_offset' in list(opts.keys()):
			y_offset = opts['y_offset']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'color_1' in list(opts.keys()):
			color_2 = opts['color_1']
		if 'color_2' in list(opts.keys()):
			color_2 = opts['color_2']
		if 'color_connector' in list(opts.keys()):
			color_connector = opts['color_connector']
		if 'linewidth_connector' in list(opts.keys()):
			linewidth_connector = opts['linewidth_connector']
		if 'len_1' in list(opts.keys()):
			len_1 = opts['len_1']
		if 'len_2' in list(opts.keys()):
			len_2 = opts['len_2']
	# Check direction (we don't use at moment)
	fwd = True
	if start > end:
		fwd = False
	# Draw the -35 site (from start to start + length of -35 site)
	pP1 = Polygon([(start, y_offset),
		           (start, y_offset+y_extent),
		           (start+len_1,y_offset+y_extent),
		           (start+len_1,y_offset)],
		            edgecolor=(1,1,1), facecolor=color_1, linewidth=0.2, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(pP1)
	# Draw the -10 site (from end-length of -10 site to end)
	pP2 = Polygon([(end-len_2, y_offset),
		           (end-len_2, y_offset+y_extent),
		           (end,y_offset+y_extent),
		           (end,y_offset)],
		            edgecolor=(1,1,1), facecolor=color_2, linewidth=0.2, zorder=11,
		            path_effects=[Stroke(joinstyle="miter")])
	ax.add_patch(pP2)
	lP1 = Line2D([start+len_1, end-len_2],
                [y_offset+(y_extent/2.0), y_offset+(y_extent/2.0)], linewidth=linewidth_connector,
                color=color_connector, zorder=10)
	ax.add_line(lP1)

	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			dpl.write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			dpl.write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer
	return start, end


##############################################################################
# CREATES THE FIGURE
##############################################################################

# Create the figure
fig = plt.figure(figsize=(8.0,3))
gs = gridspec.GridSpec(1, 1)
ax_dna = plt.subplot(gs[0])

# This is the main sequence
PprqR = 'ggacggtatcttgaaagtgtccaaccaactggacggtttaatTATTTCcaggtcttctgtcaattTAATATttgcattatctgcGGAGgcgaaattatATG'
draw_sequence(ax_dna, 0, 30, PprqR)

# Another sequence to show you can overlay them if needed with an offset
J23119 = 'TTGACAgctagctcagtcctaggTATAATgctagc'
draw_sequence(ax_dna, 42, 16, J23119)
Ptic10= 'TTGACAattaatcatcgcggctcgTATAATgtgtgg'
draw_sequence(ax_dna, 42, 13, Ptic10)
BBa_0034 = 'aaagaGGAGaaa'
draw_sequence(ax_dna, 79, 16, BBa_0034)
BBa_R0040 = 'tccctatcagtgatagagaTTGACAtccctatcagtgatagaGATACTgagcac'
draw_sequence(ax_dna, 23, 7, BBa_R0040 )
BBa_R0011 = 'ttgtgagcggataacaaTTGACAttgtgagcggataacaaGATACTgagcac'
draw_sequence(ax_dna, 25, 2, BBa_R0011)
SynRNA = 'aCCUCcuuu'
draw_sequence(ax_dna, 83, 32, SynRNA)
EcoRNA = 'aCCUCcuaa'
draw_sequence(ax_dna, 83, 18, EcoRNA)

# Create the DNAplotlib renderer and map of part types to renderering functions
dr = dpl.DNARenderer()
reg_renderers = dr.std_reg_renderers()
part_renderers = dr.SBOL_part_renderers()
part_renderers['PromoterRegion'] =  promoter_region
part_renderers['OperatorRegion'] = operator_region
part_renderers['Palindromes'] = palindrome_region


# Create the sites to draw
SynPromoters = {'type':'PromoterRegion', 'name':'region1', 'start': 42, 'end': 71, 'fwd':True, 'opts':{'y_offset':32, 'len_35':6, 'len_10':6, 'color_35':(1,1,1), 'color_10':(0.9,0.2,0.7)}}
EcoPromoters = {'type':'PromoterRegion', 'name':'region1', 'start': 42, 'end': 71, 'fwd':True, 'opts':{'y_offset':19, 'len_35':6, 'len_10':6, 'color_35':(0,0.5,0.2), 'color_10':(0.9,0.2,0.7)}}
LacO = {'type':'OperatorRegion', 'name':'operatorlac', 'start': 23, 'end': 65, 'fwd':True, 'opts':{'y_offset':6, 'len_O1':17, 'len_O2':17}}
TetO = {'type':'OperatorRegion', 'name':'operatortet', 'start': 25, 'end': 65, 'fwd':True, 'opts':{'y_offset':1, 'len_O1':17, 'len_O2':17}}
PrqRO = {'type':'OperatorRegion', 'name':'operatorprqr', 'start': 0, 'end': 37, 'fwd':True, 'opts':{'y_offset':29, 'len_O1':7, 'len_O2':7}}
Pal1 = {'type':'Palindromes', 'name':'palindrome1', 'start': 13, 'end': 39, 'fwd':True, 'opts':{'y_offset':27, 'len_O1':10, 'len_O2':10}}
Pal2 = {'type':'Palindromes', 'name':'palindrome1', 'start': 10, 'end': 51, 'fwd':True, 'opts':{'y_offset':25, 'len_O1':7, 'len_O2':7}}

# We don't use a design, so just annotate an axis
dr.annotate(ax_dna, part_renderers, SynPromoters)
dr.annotate(ax_dna, part_renderers, EcoPromoters)
dr.annotate(ax_dna, part_renderers, LacO)
dr.annotate(ax_dna, part_renderers, TetO)
dr.annotate(ax_dna, part_renderers, PrqRO)
dr.annotate(ax_dna, part_renderers, Pal1)
dr.annotate(ax_dna, part_renderers, Pal2)



# Set up bounds of the axis
ax_dna.set_xlim([-0.5, len(PprqR)+0.5])
ax_dna.set_ylim([0,40])
ax_dna.set_xticks([])
ax_dna.set_yticks([])
ax_dna.axis('off')

# Update subplot spacing
plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

# Save the figure
fig.savefig('sequence_features.pdf', transparent=True)
fig.savefig('sequence_features.png', dpi=600)

# Clear the plotting cache
plt.close('all')
