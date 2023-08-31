# IMPORT STATEMENTS

import math
import numpy as np
from scipy.stats import fisher_exact, ttest_ind
import pandas as pd
import seaborn as sns

import os
import sys
sys.path.append("/Users/rohan/public_html/Hegemon")

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import StepMiner as smn
import HegemonUtil as hu
import bone

acolor = ["#00CC00", "#D8A03D","#EC008C",
          'cyan', "#B741DC", "#808285",
          'blue', 'black', 'green', 'red',
          'orange', 'brown', 'pink', 'purple']

def getDataset(dbid, cf="/Users/rohan/public_html/Hegemon/explore.conf"):
    db = hu.Database(cf)
    h = hu.Hegemon(db.getDataset(dbid))
    h.init()
    h.initPlatform()
    h.initSurv()
    return(h)

cf = "/Users/rohan/public_html/Hegemon/explore.conf"
    
def getJeffreySamples(tn=1):
    h = getDataset("MACV52", cf)
    atype = h.getSurvName("c Title")
    ahash = {}
    if tn == 1:
        atypes = ['Eosinophils', 'B cells', 'Basophils', 'T cells', 'mast cells', 'Dendritic cells', 'Macrophages', \
                  'NK cells', 'Neutrophils', 'Th1 cells', 'Th2 cells']
        for title in atype[2:]:
            for celltype in atypes:
                if celltype in title:
                    ahash[title] = celltype
                if celltype == "Dendritic cells":
                    if "endritic cells" in title:
                        ahash[title] = "Dendritic cells"
    else:
        atypes = ['Active Eos.', 'Inactive Eos.', 'B cells', 'Basophils', 'T cells', 'mast cells', 'Dendritic cells', 'Macrophages', \
                  'NK cells', 'Neutrophils', 'Th1 cells', 'Th2 cells']
        for title in atype[2:]:
            for celltype in atypes:
                if celltype in title:
                    ahash[title] = celltype
                if celltype == "Dendritic cells":
                    if "endritic cells" in title:
                        ahash[title] = celltype
                if celltype == 'Active Eos.':
                    if 'Eosinophils_PMA' in title:
                        ahash[title] = celltype
                if celltype == 'Inactive Eos.':
                    if 'Eosinophils_control' in title:
                        ahash[title] = celltype
    return h, atype, atypes, ahash

#h = getDataset("G19", cf)

def getNorvershternSamples(tn=1):
    h = getDataset("G19", cf)
    atype = h.getSurvName("c cell type")
    ahash = {}
    if (tn==1):
        atypes =  ['Hematopoietic stem cell', 'Granulocyte','NK','Memory','T-cell', 'B-cell', 'Erythroid',\
                   'Neutrophil','Eosinophil','Basophil','Megakaryocyte','Metamyelocyte','Dendritic Cell',\
                   'myeloid progenitor','Granulocyte/monocyte progenitor','Monocyte']
        for title in atype[2:]:
            for celltype in atypes:
                if ("Megakaryocyt" in title) and not ("erythroid" in title):
                    ahash[title] = "Megakaryocyte"
                if celltype in title:
                    ahash[title] = celltype
    else:
        atypes =  ['HSC','MP','GMP','Granulocyte','Metamyelocyte','Neutrophil','Eosinophil','Basophil','NK',\
                   'Memory','T-cell','Monocyte','B-cell','Erythroid','Megakaryocyte','Dendritic Cell']
        for title in atype[2:]:
            for celltype in atypes:
                if ("Hematopoietic" in title):
                    ahash[title] = "HSC"
                if ("myeloid progenitor" in title):
                    ahash[title] = "MP"
                if ("Megakaryocyt" in title) and not ("erythroid" in title):
                    ahash[title] = "Megakaryocyte"
                if ("Granulocyte/monocyte" in title):
                    ahash[title] = "GMP"
                if ('Colony Forming Unit-Granulocyte' in title):
                    ahash[title] = "Granulocyte"
                if ('Granulocyte (Neutrophilic Metamyelocyte)' in title):
                    ahash[title] = "Metamyelocyte"
                if ('Granulocyte (Neutrophil)' in title):
                    ahash[title] = "Neutrophil"
                if celltype in title:
                    ahash[title] = celltype        
    return h, atype, atypes, ahash 

#h = getDataset("NEU6", cf)

def getAllantazSamples(tn=1):
    h = getDataset("NEU6", cf)
    atype = h.getSurvName("c cell type (ch1)")
    ahash = {}
    if tn == 1:
        atypes = ['CD19+ B cells', 'CD14+ monocytes', 'CD4+ T cells', 'CD8+ T cells',\
                  'Eosinophils', 'NK cells', 'Neutrophils']    
    else:
        atypes = ['B cells', 'monocytes', 'T cells', 'Eosinophils', 'NK cells', 'Neutrophils']
    for title in atype[2:]:
        for celltype in atypes:
            if celltype in title:
                ahash[title] = celltype
    return h, atype, atypes, ahash

def getMonacoSamples(tn=1):
    h = getDataset("NEU12", cf)
    atype = h.getSurvName("c cell type (ch1)")
    ahash = {}
    if tn == 0:
        atypes = ['Naive CD8 T cells', 'Central memory CD8 T cell', 'Effector memory CD8 T cells', \
                  'Terminal effector CD8 T cells', 'MAIT cells', 'Vd2 gd T cells', 'Non-Vd2 gd T cells', \
                  'Follicular helper T cells', 'T regulatory cells', 'Th1 cells', 'Th1/Th17 cells', 'Th17 cells', 'Th2 cells', \
                  'Naive CD4 T cells', 'Progenitor cells', 'Naive B cells', 'Non-switched memory B cells', 'Exhausted B cells', \
                  'Switched memory B cells', 'Plasmablasts', 'Classical monocytes', 'Intermediate monocytes', \
                  'Non classical monocytes', 'Natural killer cells', 'Plasmacytoid dendritic cells', 'Myeloid dendritic cells', \
                  'Low-density neutrophils', 'Low-density basophils', 'Terminal effector CD4 T cells', 'PBMCs']
        for title in atype[2:]:
            for celltype in atypes:
                if celltype in title:
                    ahash[title] = celltype    
    else:
        atypes = ['T cells', 'Progenitor', 'B cells', 'Plasmablasts', 'Monocytes', "NK cells", "Dendritic cells","Neutrophils",\
                 "Basophils", "PBMCs"]
        for title in atype[2:]:
            for celltype in atypes:
                if ("Th" in title) or ("T" in title and "cell" in title):
                    ahash[title] = 'T cells'
                if "Natural killer" in title:
                    ahash[title] = 'NK cells'
                if celltype.upper() in title.upper():
                    ahash[title] = celltype 
                if celltype in title:
                    ahash[title] = celltype 
    return h, atype, atypes, ahash

def getData(geneName, samples, tn=1):
    h, atype, atypes, ahash = samples
    geneData = h.getExprData(geneName)
    data = []
    for celltype in atypes:
        cellrange = [i for i in h.aRange() if ahash[atype[i]] == celltype]
        cellexpr = [float(geneData[i]) for i in cellrange]
        data.append(cellexpr)
    return data

def getDataDf(geneName, samples, tn=1):
    h, atype, atypes, ahash = samples
    df = pd.DataFrame(getData(geneName, samples, tn=1))
    df.index = atypes
    return df.T

def ttest(gene, samples, g1, g2=None, tn=1):
    data = getDataDf(gene, samples, tn=1)       
    expr1 = [i for i in data[g1] if i != math.nan]
    if g2 is None:
        expr2 = []
        othercols = list(data.columns)
        othercols.remove(g1)
        for col in othercols:
            vals = [i for i in data[col] if i != math.nan]
            expr2 += vals
    else:
        expr2 = [i for i in data[g2] if i != math.nan]
    res = ttest_ind(expr1, expr2, equal_var=False, nan_policy='omit')
    return res.statistic, res.pvalue
    
def ttestGenes(geneList, samples, g1, g2=None, tn=1, alternative="two-sided", alpha=0.05, printres=False):
    significant = {}
    insignificant = {}
    for i in geneList:
        try:
            stat, pval = ttest(i, samples, g1, g2, tn=1)
            if alternative=="two-sided":
                if pval < alpha:
                    significant[i] = stat, pval
                else:
                    insignificant[i] = stat, pval
            elif alternative=="less":
                if (pval/2 < alpha) and (stat < 0):
                    significant[i] = stat, pval/2
                else:
                    insignificant[i] = stat, pval/2                
            elif alternative=="greater":
                if (pval/2 < alpha) and (stat > 0):
                    significant[i] = stat, pval/2
                else:
                    insignificant[i] = stat, pval/2 
            else:
                print("Alternative must be two-sided, less or greater")            
        except Exception as e:
            #print(e)
            print("Error or could not find %s" %i)
    if printres:
        print(g1+" vs. "+str(g2), len(significant)/len(geneList))
        print(" ".join(list(significant.keys())))
    return len(significant)/len(geneList), significant, insignificant

def saveViolin(geneName, acolor, samples, save=False, tn=1, w=1.4, h=10, highlight=None):
    atypes = samples[2]
    try:
        geneData = getData(geneName, samples, tn)
        if highlight is not None:
            acolor = ["#39cdfa" if i != highlight else "#2a07f0" for i in atypes]
        params = {"dpi":150, "w":len(atypes)*w, "h":h, "acolor":acolor, "vert":1}
        if save:
            plot = bone.plotViolin2(geneData, atypes, params, geneName)
        else:
            plot = bone.plotViolin(geneData, atypes, params, geneName)
        return(plot)
    except Exception as e:
        #print(e)
        print("Error or could not find %s" % geneName)
        
def getGenes(resfile):
    res_table = pd.DataFrame(pd.read_csv(resfile, sep = '\t', header=None))
    genes = []
    for i in res_table[3]:
        gene = i.split(":")[0] 
        if " /// " in gene:
            gene = gene.split(" /// ")[0]
        genes.append(gene)
    return genes
    
def testGenes(geneList, samples, tn=1, w=1.4, h=10, acolor=acolor, highlight=None):
    for gene in hu.uniq(geneList):
        try:
            print(gene)
            saveViolin(gene, acolor, samples, False, tn, w, h, highlight)
        except:
            try:
                saveViolin(gene.upper(), acolor, samples, False, tn, w, h, highlight)
            except Exception as e:
                #print(e)
                print("Could not find", gene)
                
def testGenesFromFile(resfile, samples, tn=1, w=1.4, h=10, acolor=acolor, highlight=None):
    testGenes(getGenes(resfile), samples, tn, w, h, acolor, highlight)

def TtestGenesFromFile(resfile, samples, g1, g2=None, tn=1, alternative="two-sided", alpha=0.05, printres=False):
    return ttestGenes(getGenes(resfile), samples, g1, g2, tn, alternative, alpha, printres)


# Saving plots to pdf for figures

def getPDF(cfile):
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(cfile)
    return pdf

def closePDF(pdf):
    import datetime
    d = pdf.infodict()
    d['Title'] = 'Plots'
    d['Author'] = 'Rohan Subramanian'
    d['Subject'] = "BECC"
    d['Keywords'] = 'Neutrophils'
    d['CreationDate'] = datetime.datetime(2023, 7, 27)
    d['ModDate'] = datetime.datetime.today()
    pdf.close()

def testGenesToPDF(geneList, cfile, samples, tn=1, w=1.4, h=10, highlight=None):
    plotPDF = getPDF(cfile)
    for gene in hu.uniq(geneList):
        try:
            print(gene)
            plotPDF.savefig(saveViolin(gene, acolor, samples, True, tn, w, h, highlight))
        except:
            plotPDF.savefig(saveViolin(gene.upper(), acolor, samples, True, tn, w, h, highlight))
    closePDF(plotPDF)

def testFileGenesToPDF(resfile, cfile, samples, tn=1, w=1.4, h=10, acolor=acolor, highlight=None):
    testGenes(getGenes(resfile), cfile, samples, tn, w, h, acolor, highlight)
    
print("ran ViolinPlot.py")