# Visualization (Fig. 3F)

import numpy as np
import hicstraw
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt 
import seaborn as sns

hic = hicstraw.HiCFile("adult_cd8.hic")

matrix_object_chr4 = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 25000)
numpy_matrix_chr4 = matrix_object_chr4.getRecordsAsMatrix(4000000, 12000000, 4000000, 12000000)
REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
def plot_hic_map(dense_matrix, maxcolor):
    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
    #plt.show()
    plt.savefig("hicmap_plots/res25k_chr4_4mb12mb.pdf", bbox_inches='tight')
plot_hic_map(numpy_matrix_chr4, 30)
