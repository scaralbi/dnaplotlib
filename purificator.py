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
J100 = {'type':'Promoter', 'name':'prom1', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'J23119', 'label_y_offset':-5}}
Pr = {'type':'Promoter', 'name':'prom2', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'PprqR', 'label_y_offset':-8}}
termP = {'type':'Terminator', 'name':'term', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'TclpPX', 'label_y_offset':-8}}
termL = {'type':'Terminator', 'name':'term', 'fwd':True, 'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'Tlac', 'label_y_offset':-8}}

rbs = {'type':'RBS', 'name':'rbs', 'fwd':True}

prqR = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['red'], 'label':'prqR-6xHis', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
Plac = {'type':'Promoter', 'name':'Plac', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.25), 'label':'Plac', 'label_y_offset':-8}}
lacZ = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':Cmap(0.25), 'label':'lacZ', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}


Neg = [sp, Plac, rbs, lacZ, termL, sp]
Pos = [sp, J100, rbs, prqR,termP, sp]



# Create the figure
fig = plt.figure(figsize=(4,1))
gs = gridspec.GridSpec(1, 2)
ax_dna1 = plt.subplot(gs[0])
ax_dna2 = plt.subplot(gs[1])


# Redender the DNA to axis
start, end = dr.renderDNA(ax_dna1, Neg, part_renderers)
ax_dna1.set_xlim([start, end])
ax_dna1.set_ylim([-18,20])
ax_dna1.set_aspect('equal')
ax_dna1.set_xticks([])
ax_dna1.set_yticks([])
ax_dna1.axis('off')
start, end = dr.renderDNA(ax_dna2, Pos, part_renderers)
ax_dna2.set_xlim([start, end])
ax_dna2.set_ylim([-18,20])
ax_dna2.set_aspect('equal')
ax_dna2.set_xticks([])
ax_dna2.set_yticks([])
ax_dna2.axis('off')



# Update subplot spacing
plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)

# Save the figure
fig.savefig('hispur.png', transparent=True, dpi=600)

# Clear the plotting cache
plt.close('all')