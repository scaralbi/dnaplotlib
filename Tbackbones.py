import math
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches


import matplotlib as mpl
from matplotlib import cm
from colorspacious import cspace_converter
from collections import OrderedDict



cmaps = OrderedDict()



__author__  = 'Alberto Scarampi <as2945@cam.ac.uk>, Howe Lab, Cambridge'
__license__ = 'MIT'
__version__ = '1.0'


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

            
GreenCmap = cm.get_cmap('Reds', 4)
Cmap = cm.get_cmap('Pastel2', 4)



def sbol_recombinase1 (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	""" SBOL recombinase site renderer - forward direction
	"""
	# Default parameters
	color = (0,0,0)
	color2 = (0,0,0)
	start_pad = 0.0
	end_pad = 0.0
	x_extent = 6.0
	y_extent = 6.0
	linestyle = '-'
	# Update default parameters if provided
	if opts != None:
		if 'start_pad' in list(opts.keys()):
			start_pad = opts['start_pad']
		if 'end_pad' in list(opts.keys()):
			end_pad = opts['end_pad']
		if 'x_extent' in list(opts.keys()):
			x_extent = opts['x_extent']
		if 'y_extent' in list(opts.keys()):
			y_extent = opts['y_extent']
		if 'linestyle' in list(opts.keys()):
			linestyle = opts['linestyle']
		if 'linewidth' in list(opts.keys()):
			linewidth = opts['linewidth']
		if 'scale' in list(opts.keys()):
			scale = opts['scale']
		if 'color' in list(opts.keys()):
			color = opts['color']
		if 'color2' in list(opts.keys()):
			color2 = opts['color2']
	# Check direction add start padding
	final_end = end
	final_start = prev_end
	y_lower = -1 * y_extent/2
	y_upper = y_extent/2
	if start > end:
		start = prev_end+end_pad+x_extent+linewidth
		end = prev_end+end_pad
		final_end = start+start_pad
		color = color2
	else:
		start = prev_end+start_pad+linewidth
		end = start+x_extent
		final_end = end+end_pad
	# Draw the site
	p1 = Polygon([(start, y_lower), 
		          (start, y_upper),
		          (end,0)],
		          edgecolor=(0,0,0), facecolor=color, linewidth=linewidth, zorder=11, 
		          path_effects=[Stroke(joinstyle="miter")])		
	ax.add_patch(p1)
	# Add a label if needed
	if opts != None and 'label' in list(opts.keys()):
		if final_start > final_end:
			write_label(ax, opts['label'], final_end+((final_start-final_end)/2.0), opts=opts)
		else:
			write_label(ax, opts['label'], final_start+((final_end-final_start)/2.0), opts=opts)
	# Return the final start and end positions to the DNA renderer
	if final_start > final_end:
		return prev_end, final_start
	else:
		return prev_end, final_end

# Function to calculate darker colour
def dark (col, fac=2.0):
	return (col[0]/fac, col[1]/fac, col[2]/fac)


# Global line width
lw = 1.0


# Create the DNAplotlib renderer
dr = dpl.DNARenderer()

# Use default renderers and append our custom ones for recombinases
reg_renderers = dr.std_reg_renderers()
part_renderers = dr.SBOL_part_renderers()
part_renderers['RecombinaseSite'] = sbol_recombinase1


# Create the construct programmably to plot
sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':1}}


Pspec = {'type':'Promoter', 'name':'Pspec', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.25), 'label':'Pspec', 'label_y_offset':-8}}
Psac = {'type':'Promoter', 'name':'Psac', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.5), 'label':'Psac', 'label_y_offset':-8}}
Plac = {'type':'Promoter', 'name':'Plac', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.25), 'label':'Plac', 'label_y_offset':-8}}

RBS = {'type':'RBS', 'name':'rbs', 'fwd':True, 'opts':{'color':col_map['gray']}}
T = {'type':'Terminator', 'name':'term', 'fwd':True}

lacZ = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':Cmap(0.25), 'label':'lacZ', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
sacB = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':Cmap(0.5), 'label':'sacB', 'label_x_offset':-0.3, 'label_y_offset':-0.5, 'label_style':'italic'}}
specR= {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':Cmap(0.75), 'label':'specR', 'label_x_offset':-0.3, 'label_y_offset':-0.5, 'label_style':'italic'}}
scar= {'type':'Scar', 'name':'scarry', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'EL3', 'label_y_offset':-8}}


ACGG= {'type':'5StickyRestrictionSite', 'name':'3SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'ACGG', 'label_x_offset':-1, 'label_y_offset':-8}}
GGGA= {'type':'3StickyRestrictionSite', 'name':'5SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'GGGA', 'label_x_offset':1, 'label_y_offset':6}}

TTAC= {'type':'5StickyRestrictionSite', 'name':'5SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'TTAC', 'label_x_offset':-1, 'label_y_offset':6}}
CCCT= {'type':'3StickyRestrictionSite', 'name':'3SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'CCCT', 'label_x_offset':1, 'label_y_offset':-8}}


lac = [ACGG, Plac, RBS, lacZ, T, GGGA]
pEL3 = [TTAC, scar, CCCT]
sac = [Psac, RBS, sacB, T]
spec = [Pspec, RBS, specR, T]


# Create the figure
fig = plt.figure(figsize=(15,4))
gs = gridspec.GridSpec(2, 2)
ax_dna1 = plt.subplot(gs[0])
ax_dna2 = plt.subplot(gs[1])
ax_dna3 = plt.subplot(gs[2])
ax_dna4 = plt.subplot(gs[3])


# Redender the DNA to axis
start, end = dr.renderDNA(ax_dna1, lac, part_renderers)
ax_dna1.set_xlim([start, end])
ax_dna1.set_ylim([-18,20])
ax_dna1.set_aspect('equal')
ax_dna1.set_xticks([])
ax_dna1.set_yticks([])
ax_dna1.axis('off')
start, end = dr.renderDNA(ax_dna2, pEL3, part_renderers)
ax_dna2.set_xlim([start, end])
ax_dna2.set_ylim([-18,20])
ax_dna2.set_aspect('equal')
ax_dna2.set_xticks([])
ax_dna2.set_yticks([])
ax_dna2.axis('off')
start, end = dr.renderDNA(ax_dna3, sac, part_renderers)
ax_dna3.set_xlim([start, end])
ax_dna3.set_ylim([-18,20])
ax_dna3.set_aspect('equal')
ax_dna3.set_xticks([])
ax_dna3.set_yticks([])
ax_dna3.axis('off')
start, end = dr.renderDNA(ax_dna4, spec, part_renderers)
ax_dna4.set_xlim([start, end])
ax_dna4.set_ylim([-18,20])
ax_dna4.set_aspect('equal')
ax_dna4.set_xticks([])
ax_dna4.set_yticks([])
ax_dna4.axis('off')



# Update subplot spacing
plt.subplots_adjust(hspace=0.01, left=0.04, right=0.9, top=0.6, bottom=0.01)



# Save the figure
fig.savefig('Tback.png', transparent=True,  dpi = 600)

# Clear the plotting cache
plt.close('all')