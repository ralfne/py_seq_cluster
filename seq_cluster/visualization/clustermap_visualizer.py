import ntpath
import os
import pandas as pd, seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt
import numpy as np
from logger.Logger import StdOutLogger
from seq_cluster.clustering.pairwise_distance_matrix import PairwiseDistanceMatrix


class ClustermapVisualizer(object):
    _FIGURE_FILENAME_SVG = 'figure.svg'
    _FIGURE_FILENAME_PNG = 'figure.png'
    _KEYS_FILENAME = 'keys.npy'
    _LINKAGE_FILENAME = 'linkage.npy'
    _LINKAGE_AUX_FILENAME = 'linkage_aux.npy'

    def __init__(self, distance_dataframe, linkage_method, distance_dataframe_aux=None,
                 linkage_method_aux=None, scaling=None,
                 fig_size_x=20, fig_size_y=20,fig_size_font=10, logger=StdOutLogger()):
        self._logger = logger
        self._linkage_method = linkage_method
        self._distance_dataframe = distance_dataframe
        self._distance_dataframe_aux = distance_dataframe_aux
        self._linkage_method_aux = linkage_method_aux
        self._scaling = scaling
        self._fig_size_x = fig_size_x
        self._fig_size_y = fig_size_y
        self._fig_size_font = fig_size_font

    def run(self, title=None, filename=None):
        df = self._distance_dataframe
        squareform = sp.distance.squareform(self._distance_dataframe)
        linkage = hc.linkage(squareform, method=self._linkage_method)
        self._possibly_save_keys(filename)
        self._possibly_save_linkage(linkage, filename, use_aux=False)
        self._logger.log('Keys and linkage files saved...', includeTimestamp=True, onlyIfVerbose=False)
        linkage_aux = linkage
        if self._distance_dataframe_aux is not None:
            squareform_aux = sp.distance.squareform(self._distance_dataframe_aux)
            linkage_aux = hc.linkage(squareform_aux, method=self._linkage_method_aux)
            self._possibly_save_linkage(linkage_aux, filename, use_aux=True)
            self._logger.log('Aux linkage file saved...', includeTimestamp=True, onlyIfVerbose=False)
            df = PairwiseDistanceMatrix.combine(self._distance_dataframe, self._distance_dataframe_aux)
        cm = sns.clustermap(df, row_linkage=linkage, col_linkage=linkage_aux,
                            figsize=(self._fig_size_x, self._fig_size_y))

        cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xmajorticklabels(), fontsize=self._fig_size_font)

        custom_title = title
        if custom_title is None:
            if filename is None:
                custom_title = ''
            else:
                custom_title = self._path_leaf(filename)
        self._format_plot(cm, custom_title)
        self._scale(cm.fig)
        self._possibly_save_plot(cm, filename)

    def _format_plot(self, cm, title):
        cm.fig.suptitle(title)
        plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    def _possibly_save_keys(self, filename):
        if filename is not None:
            l = len(self._distance_dataframe)
            array = np.arange(l)
            i = 0
            for key, row in self._distance_dataframe.iterrows():
                value = None
                parts = key.split('_')
                if len(parts)==1:
                    value = key
                else:
                    for p in parts:
                        if unicode(p).isnumeric():
                            value=p
                            break
                    if value is None:
                        value = -1
                array[i] = value
                i += 1
            fn = os.path.join(filename, ClustermapVisualizer._KEYS_FILENAME)
            np.save(fn, array)

    def _possibly_save_linkage(self, linkage, filename, use_aux=False):
        if filename is not None:
            if use_aux:
                fn = os.path.join(filename, ClustermapVisualizer._LINKAGE_AUX_FILENAME)
            else:
                fn = os.path.join(filename, ClustermapVisualizer._LINKAGE_FILENAME)
            np.save(fn, linkage, allow_pickle=False)

    def _possibly_save_plot(self, cm, filename):
        if filename is not None:
            fn_png = os.path.join(filename, ClustermapVisualizer._FIGURE_FILENAME_PNG)
            fn_svg = os.path.join(filename, ClustermapVisualizer._FIGURE_FILENAME_SVG)
            fig = cm.fig
            fig.savefig(fn_png)
            fig.savefig(fn_svg)

    def _scale(self, fig):
        if self._scaling is not None:
            scalting_alt = 1 - self._scaling
            fig.subplots_adjust(left=self._scaling, bottom=self._scaling, right=scalting_alt, top=scalting_alt)

    def _path_leaf(self, path):
        head, tail = ntpath.split(path)
        out = tail or ntpath.basename(head)
        filename, file_extension = os.path.splitext(out)
        return filename
