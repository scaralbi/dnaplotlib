import math
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Polygon, Ellipse, Wedge, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patheffects import Stroke
import matplotlib.patches as patches

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


# Create the construct programmably to plot
sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':1}}
int = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':5}}
J19 = {'type':'Promoter', 'name':'prom1', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'J23119', 'label_y_offset':-5}}
J105 = {'type':'Promoter', 'name':'prom2', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['gray'], 'label':'J23105', 'label_y_offset':-8}}
term = {'type':'Terminator', 'name':'term', 'fwd':True}
tag = {'type':'UNS', 'name':'tag', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['orange'], 'label':'c-myc', 'label_y_offset':-8}}
rbs = {'type':'RBS', 'name':'rbs', 'fwd':True, 'opts':{'color':col_map['gray']}}
NOX = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['green'], 'label':'nox', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
NOXmyc = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['orange'], 'label':'nox+c-myc', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
Pran = {'type':'Promoter', 'name':'pran', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'Prandom', 'label_y_offset':-8}}


n1 = [sp, J19, rbs, NOX, term, sp]
n2 = [sp, J105, rbs, NOX, term, sp]
n3 = [sp, J19, rbs, NOX, tag, term, sp]
n4 = [sp, J105, rbs, NOX, tag, term, sp]

#arc1 = {'type':'Repression', 'from_part':prqR, 'to_part':ope, 'opts':{'color':col_map['black'],'linewidth':lw}}

#reg1 = [arc1]

# Create the figure
fig = plt.figure(figsize=(4,1))
gs = gridspec.GridSpec(2, 2)
ax_dna1 = plt.subplot(gs[0])
ax_dna2 = plt.subplot(gs[1])
ax_dna3 = plt.subplot(gs[2])
ax_dna4 = plt.subplot(gs[3])

# Redender the DNA to axis
start, end = dr.renderDNA(ax_dna1, n1, part_renderers)
ax_dna1.set_xlim([start, end])
ax_dna1.set_ylim([-18,20])
ax_dna1.set_aspect('equal')
ax_dna1.set_xticks([])
ax_dna1.set_yticks([])
ax_dna1.axis('off')
start, end = dr.renderDNA(ax_dna2, n3, part_renderers)
ax_dna2.set_xlim([start, end])
ax_dna2.set_ylim([-18,20])
ax_dna2.set_aspect('equal')
ax_dna2.set_xticks([])
ax_dna2.set_yticks([])
ax_dna2.axis('off')
start, end = dr.renderDNA(ax_dna3, n2, part_renderers)
ax_dna3.set_xlim([start, end])
ax_dna3.set_ylim([-18,20])
ax_dna3.set_aspect('equal')
ax_dna3.set_xticks([])
ax_dna3.set_yticks([])
ax_dna3.axis('off')
start, end = dr.renderDNA(ax_dna4, n4, part_renderers)
ax_dna4.set_xlim([start, end])
ax_dna4.set_ylim([-18,20])
ax_dna4.set_aspect('equal')
ax_dna4.set_xticks([])
ax_dna4.set_yticks([])
ax_dna4.axis('off')



# Update subplot spacing
plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)

# Save the figure
fig.savefig('nox1.pdf', transparent=True)
fig.savefig('nox1.png', dpi=600)

# Clear the plotting cache
plt.close('all')