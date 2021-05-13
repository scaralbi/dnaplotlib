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
sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':3}}
int = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':5}}
J19 = {'type':'Promoter', 'name':'prom1', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.5), 'label':'J23119', 'label_y_offset':-5}}
Pr = {'type':'Promoter', 'name':'prom2', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.5), 'label':'PprqR', 'label_y_offset':-8}}
term = {'type':'Terminator', 'name':'term', 'fwd':True}
ope = {'type':'Operator', 'name':'term', 'fwd':True}
ope1 = {'type':'Operator', 'name':'term', 'fwd':True}
rbs = {'type':'RBS', 'name':'rbs', 'fwd':True, 'opts':{'color':col_map['gray']}}
eYFP = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['yellow'], 'label':'eYFP', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
prqR = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':col_map['red'], 'label':'prqR', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
Pran = {'type':'Promoter', 'name':'pran', 'fwd':True,'opts':{'linewidth':lw, 'color':col_map['black'], 'label':'Prandom', 'label_y_offset':-8}}
lacZ = {'type':'CDS', 'name':'cds', 'fwd':True, 'opts':{'color':Cmap(0.25), 'label':'lacZ', 'label_x_offset':-2, 'label_y_offset':-0.5, 'label_style':'italic'}}
Plac = {'type':'Promoter', 'name':'Plac', 'fwd':True,'opts':{'linewidth':lw, 'color':Cmap(0.25), 'label':'Plac', 'label_y_offset':-8}}


pos = [sp, ope, Plac, ope, rbs, prqR, term, sp, J19, rbs, eYFP, term, sp]
neg = [sp, ope, Plac, ope, rbs, prqR, term, sp, Pran, rbs, eYFP, term]
test1 = [sp, ope, Plac, ope, rbs, prqR, term, sp, ope1, Pr, rbs, eYFP, term]
test2 = [sp, ope, Plac, ope, rbs, prqR, term, sp, ope1, Pr, ope, rbs, eYFP, term]


arc1 = {'type':'Repression', 'from_part':prqR, 'to_part':ope, 'opts':{'color':col_map['black'],'linewidth':lw}}
arc2 = {'type':'Repression', 'from_part':prqR, 'to_part':ope1, 'opts':{'color':col_map['black'],'linewidth':lw}}

reg1 = [arc2]

# Create the figure
fig = plt.figure(figsize=(8,2))
gs = gridspec.GridSpec(4, 2)
ax_dna1 = plt.subplot(gs[0])
ax_dna2 = plt.subplot(gs[1])
ax_dna3 = plt.subplot(gs[2])
ax_dna4 = plt.subplot(gs[3])
ax_dna5 = plt.subplot(gs[4])
ax_dna6 = plt.subplot(gs[5])
ax_dna7 = plt.subplot(gs[6])
ax_dna8 = plt.subplot(gs[7])

# Redender the DNA to axis
start, end = dr.renderDNA(ax_dna1, neg, part_renderers)
ax_dna1.set_xlim([start, end])
ax_dna1.set_ylim([-18,20])
ax_dna1.set_aspect('equal')
ax_dna1.set_xticks([])
ax_dna1.set_yticks([])
ax_dna1.axis('off')
start, end = dr.renderDNA(ax_dna2, neg, part_renderers)
ax_dna2.set_xlim([start, end])
ax_dna2.set_ylim([-18,20])
ax_dna2.set_aspect('equal')
ax_dna2.set_xticks([])
ax_dna2.set_yticks([])
ax_dna2.axis('off')
start, end = dr.renderDNA(ax_dna3, pos, part_renderers)
ax_dna3.set_xlim([start, end])
ax_dna3.set_ylim([-18,20])
ax_dna3.set_aspect('equal')
ax_dna3.set_xticks([])
ax_dna3.set_yticks([])
ax_dna3.axis('off')
start, end = dr.renderDNA(ax_dna4, pos, part_renderers)
ax_dna4.set_xlim([start, end])
ax_dna4.set_ylim([-18,20])
ax_dna4.set_aspect('equal')
ax_dna4.set_xticks([])
ax_dna4.set_yticks([])
ax_dna4.axis('off')
start, end = dr.renderDNA(ax_dna5, test1, part_renderers)
ax_dna5.set_xlim([start, end])
ax_dna5.set_ylim([-18,20])
ax_dna5.set_aspect('equal')
ax_dna5.set_xticks([])
ax_dna5.set_yticks([])
ax_dna5.axis('off')
start, end = dr.renderDNA(ax_dna6, test1, part_renderers, regs=reg1, reg_renderers=reg_renderers)
ax_dna6.set_xlim([start, end])
ax_dna6.set_ylim([-18,20])
ax_dna6.set_aspect('equal')
ax_dna6.set_xticks([])
ax_dna6.set_yticks([])
ax_dna6.axis('off')
start, end = dr.renderDNA(ax_dna7, test2, part_renderers)
ax_dna7.set_xlim([start, end])
ax_dna7.set_ylim([-18,20])
ax_dna7.set_aspect('equal')
ax_dna7.set_xticks([])
ax_dna7.set_yticks([])
ax_dna7.axis('off')
start, end = dr.renderDNA(ax_dna8, test2, part_renderers, regs=reg1, reg_renderers=reg_renderers)
ax_dna8.set_xlim([start, end])
ax_dna8.set_ylim([-18,20])
ax_dna8.set_aspect('equal')
ax_dna8.set_xticks([])
ax_dna8.set_yticks([])
ax_dna8.axis('off')

# Update subplot spacing
plt.subplots_adjust(hspace=0.01, left=0.05, right=0.95, top=0.92, bottom=0.01)

# Save the figure
fig.savefig('operators.png', transparent=True,  dpi=600)

# Clear the plotting cache
plt.close('all')
