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

            
GreenCmap = cm.get_cmap('Greens', 5)
prqRmap = cm.get_cmap('winter', 2)
prqAmap = cm.get_cmap('summer', 2)


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

def sbol_recombinase2 (ax, type, num, start, end, prev_end, scale, linewidth, opts):
	""" SBOL recombinase site renderer - reverse direction
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
		else:
			if 'color' in list(opts.keys()):
				r2 = float(color[0]) / 2
				g2 = float(color[1]) / 2
				b2 = float(color[2]) / 2
				color2 = (r2,g2,b2)
	# Check direction add start padding
	final_end = end
	final_start = prev_end
	y_lower = -1 * y_extent/2
	y_upper = y_extent/2
	if start > end:
		start = prev_end+end_pad+x_extent+linewidth
		end = prev_end+end_pad
		final_end = start+start_pad
		temp = color
		color = color2
		color2 = temp
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
	midpoint = (end + start) / 2
	hypotenuse = math.sqrt( (y_extent/2)**2 + (x_extent)**2 )
	hypotenuse2 = hypotenuse / 2
	cosineA = (y_extent/2) / hypotenuse
	f = hypotenuse2 * cosineA
	p2 = Polygon([(midpoint, -1*f), 
		          (midpoint, f),
		          (end,0)],
		          edgecolor=(0,0,0), facecolor=color2, linewidth=linewidth, zorder=12, 
		          path_effects=[Stroke(joinstyle="miter")]) 	
	ax.add_patch(p1)
	ax.add_patch(p2)	
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
part_renderers['RecombinaseSite2'] = sbol_recombinase2


# Create the construct programmably to plot
sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':1}}
hrtR = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':5, 'label':'prqR HRT', 'label_y_offset':-5}}
hrtRA = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':5, 'label':'prqRA HRT', 'label_y_offset':-5}}
hrtA = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':5, 'label':'prqA HRT', 'label_y_offset':-5}}

J105 = {'type':'Promoter', 'name':'J23105', 'fwd':True,'opts':{'linewidth':lw, 'color':GreenCmap(0.2), 'label':'J23105', 'label_y_offset':-5}}
J114 = {'type':'Promoter', 'name':'J23114', 'fwd':True,'opts':{'linewidth':lw, 'color':GreenCmap(0.4), 'label':'J23114', 'label_y_offset':-8}}
J117 = {'type':'Promoter', 'name':'J12117', 'fwd':True,'opts':{'linewidth':lw, 'color':GreenCmap(0.6), 'label':'J23117', 'label_y_offset':-8}}
J119 = {'type':'Promoter', 'name':'J23119', 'fwd':True,'opts':{'linewidth':lw, 'color':GreenCmap(0.8), 'label':'J23119', 'label_y_offset':-8}}
Plac = {'type':'Promoter', 'name':'Plac', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['green'], 'label':'Plac', 'label_y_offset':-8}}

rbs = {'type':'RBS', 'name':'rbs', 'fwd':True, 'opts':{'color':col_map['gray']}}
TluxI = {'type':'Terminator', 'name':'term', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'TluxI', 'label_y_offset':-8}}
TcpcG1 = {'type':'Terminator', 'name':'term', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'TcpcG1', 'label_y_offset':-8}}


cas12a = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['purple'], 'label':'cas12a', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
gRNA_R = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['red'], 'label':'prqR gRNA', 'label_x_offset':-0.3, 'label_y_offset':-0.5, 'label_style':'italic'}}
gRNA_A = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['orange'], 'label':'prqA gRNA', 'label_x_offset':-0.3, 'label_y_offset':-0.5, 'label_style':'italic'}}
prqR = {'type':'CDS', 'name':'prqR', 'fwd':True, 'opts':{'color':col_map['red'], 'label':'prqR', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
prqA = {'type':'CDS', 'name':'prqA', 'fwd':True, 'opts':{'color':col_map['orange'], 'label':'prqA', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}

recRf = {'type':'RecombinaseSite',  'name':'a1',  'fwd':True,   'opts':{'color':prqRmap(0), 'color2':prqRmap(0), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}
recRr = {'type':'RecombinaseSite',  'name':'a2',  'fwd':False,  'opts':{'color':prqRmap(1), 'color2':prqRmap(1), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}
recAf = {'type':'RecombinaseSite',  'name':'a3',  'fwd':True,   'opts':{'color':prqAmap(0), 'color2':prqAmap(0), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}
recAr = {'type':'RecombinaseSite',  'name':'a4',  'fwd':False,  'opts':{'color':prqAmap(1), 'color2':prqAmap(1), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}
recRAf = {'type':'RecombinaseSite',  'name':'a1',  'fwd':True,   'opts':{'color':prqRmap(0), 'color2':prqRmap(0), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}
recRAr = {'type':'RecombinaseSite',  'name':'a2',  'fwd':False,  'opts':{'color':prqAmap(1), 'color2':prqAmap(1), 'x_extent':16, 'y_extent':12, 'start_pad':3, 'end_pad':3}}

TGCC= {'type':'5StickyRestrictionSite', 'name':'5SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'TGCC', 'label_x_offset':-1, 'label_y_offset':6}}
CGTT= {'type':'3StickyRestrictionSite', 'name':'3SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'CGTT', 'label_x_offset':1, 'label_y_offset':-8}}

GCAA= {'type':'5StickyRestrictionSite', 'name':'5SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'GCAA', 'label_x_offset':-1, 'label_y_offset':6}}
TGAT= {'type':'3StickyRestrictionSite', 'name':'3SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'TGAT', 'label_x_offset':1, 'label_y_offset':-8}}

ACTA= {'type':'5StickyRestrictionSite', 'name':'5SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'ACTA', 'label_x_offset':-1, 'label_y_offset':6}}
AATG= {'type':'3StickyRestrictionSite', 'name':'3SRS', 'fwd':True, 'opts':{'color':col_map['black'], 'label':'AATG', 'label_x_offset':1, 'label_y_offset':-8}}


pJML1017 = [TGCC,  J105, rbs, cas12a, TcpcG1, CGTT]
pJML1016 = [TGCC,  J114, rbs, cas12a, TcpcG1, CGTT]
pJML1015 = [TGCC,  J117, rbs, cas12a, TcpcG1,  CGTT]
pAS104 = [GCAA,  J119, gRNA_R, TluxI,  TGAT]
pAS109 = [GCAA,  J119, gRNA_A, TluxI,  TGAT]
pAS101 = [ACTA, recRf, hrtR, recRr, AATG]
pAS102 = [ACTA,  recRAf, hrtRA, recRAr,  AATG]
pAS106 = [ACTA,  recAf, hrtA, recAr,  AATG]




# Create the figure
fig = plt.figure(figsize=(12,4))
gs = gridspec.GridSpec(2, 4)
ax_dna1 = plt.subplot(gs[0])
ax_dna2 = plt.subplot(gs[1])
ax_dna3 = plt.subplot(gs[2])
ax_dna4 = plt.subplot(gs[3])
ax_dna5 = plt.subplot(gs[4])
ax_dna6 = plt.subplot(gs[5])
ax_dna7 = plt.subplot(gs[6])
ax_dna8 = plt.subplot(gs[7])

# Redender the DNA to axis
start, end = dr.renderDNA(ax_dna1, pJML1017, part_renderers)
ax_dna1.set_xlim([start, end])
ax_dna1.set_ylim([-18,20])
ax_dna1.set_aspect('equal')
ax_dna1.set_xticks([])
ax_dna1.set_yticks([])
ax_dna1.axis('off')
start, end = dr.renderDNA(ax_dna2, pJML1016, part_renderers)
ax_dna2.set_xlim([start, end])
ax_dna2.set_ylim([-18,20])
ax_dna2.set_aspect('equal')
ax_dna2.set_xticks([])
ax_dna2.set_yticks([])
ax_dna2.axis('off')
start, end = dr.renderDNA(ax_dna3, pJML1015, part_renderers)
ax_dna3.set_xlim([start, end])
ax_dna3.set_ylim([-18,20])
ax_dna3.set_aspect('equal')
ax_dna3.set_xticks([])
ax_dna3.set_yticks([])
ax_dna3.axis('off')
start, end = dr.renderDNA(ax_dna4, pAS104, part_renderers)
ax_dna4.set_xlim([start, end])
ax_dna4.set_ylim([-18,20])
ax_dna4.set_aspect('equal')
ax_dna4.set_xticks([])
ax_dna4.set_yticks([])
ax_dna4.axis('off')
start, end = dr.renderDNA(ax_dna5, pAS109, part_renderers)
ax_dna5.set_xlim([start, end])
ax_dna5.set_ylim([-18,20])
ax_dna5.set_aspect('equal')
ax_dna5.set_xticks([])
ax_dna5.set_yticks([])
ax_dna5.axis('off')
start, end = dr.renderDNA(ax_dna6, pAS101, part_renderers)
ax_dna6.set_xlim([start, end])
ax_dna6.set_ylim([-18,20])
ax_dna6.set_aspect('equal')
ax_dna6.set_xticks([])
ax_dna6.set_yticks([])
ax_dna6.axis('off')
start, end = dr.renderDNA(ax_dna7, pAS102, part_renderers)
ax_dna7.set_xlim([start, end])
ax_dna7.set_ylim([-18,20])
ax_dna7.set_aspect('equal')
ax_dna7.set_xticks([])
ax_dna7.set_yticks([])
ax_dna7.axis('off')
start, end = dr.renderDNA(ax_dna8, pAS106, part_renderers)
ax_dna8.set_xlim([start, end])
ax_dna8.set_ylim([-18,20])
ax_dna8.set_aspect('equal')
ax_dna8.set_xticks([])
ax_dna8.set_yticks([])
ax_dna8.axis('off')



# Update subplot spacing
plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)



# Save the figure
fig.savefig('crisprplasmids.png', dpi = 600)

# Clear the plotting cache
plt.close('all')