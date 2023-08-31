# Class to get data and generate plots for analyzing neutrophils

import sys
sys.path.append("/Users/rohan/public_html/Hegemon")
import StepMiner as smn
import HegemonUtil as hu
import bone
print("ran NeuAnalysis")
acolor = ["#00CC00", "#D8A03D","#EC008C",
          'cyan', "#B741DC", "#808285",
          'blue', 'black', 'green', 'red',
          'orange', 'brown', 'pink', 'purple']

pathdir = "/Users/rohan/public_html/Hegemon/Data2/Data/Analysis/Neutrophil"
cf = "/Users/rohan/public_html/Hegemon/explore.conf"

class NeuAnalysis(bone.IBDAnalysis):
    
    def __init__(self):
        bone.IBDAnalysis.__init__(self)
    
    def getGlobal(self, tn=1):
        self.prepareData("GL1")
        df = bone.pd.read_csv(pathdir+"/Data/annotation.csv", index_col=0)
        atype = ["", ""] + [str(df['tissue_state'][k]) if k in df['tissue_state']
                            else "nan" for k in self.h.headers[2:]]
        atypes = ["nan", 'solid', 'mix', 'liquid', 'cellline', 'unknown', 'lymphoma',
           'lymph', 'cellline_p3', 'bonemarrow', 'saliva']
        ahash = {}
        if tn == 2:
            atypes = ['solid', 'liquid', 'cellline']
        self.initData(atype, atypes, ahash)
        return

    def getGlobalNeu(self, tn=1):
        self.prepareData("GL1")
        cfile = "/Data/human-gpl570-tissue.txt"
        df = bone.pd.read_csv(cfile, sep="\t", index_col=0)
        atype = ["", ""] + [str(df['c Tissue'][k]) if k in df['c Tissue']
                            else "" for k in self.h.headers[2:]]
        atypes = ["Other", 'Neu']
        ahash = {}
        for title in atype[2:]:
            if "neutrophil" in title.lower():
                ahash[title] = 1
            else:
                ahash[title] = 0
        self.initData(atype, atypes, ahash)
        return

    def getNeuAtypes(self, atype, tn):
        ahash = {}
        if tn == 1:
            atypes = ['Other', 'Eosinophil']
            for title in atype[2:]:
                if 'Eosinophil' in title:
                    ahash[title] = 1
                else:
                    ahash[title] = 0
        elif tn == 2:
            atypes = ['Other', 'Neutrophil']
            for title in atype[2:]:
                if 'Neutrophil' in title:
                    ahash[title] = 1
                else:
                    ahash[title] = 0
        elif tn == 3:
            atypes = ['Other', 'Eosinophil', 'Neutrophil']
            for title in atype[2:]:
                if 'Neutrophil' in title:
                    ahash[title] = 2
                elif 'Eosinophil' in title:
                    ahash[title] = 1    
                else:
                    ahash[title] = 0
        elif tn == 4:
            atypes = ['Other', 'Neutrophil', 'Eosinophil']
            for title in atype[2:]:
                if 'Neutrophil' in title:
                    ahash[title] = 1
                elif 'Eosinophil' in title:
                    ahash[title] = 2    
                else:
                    ahash[title] = 0
        elif tn == 5:
            atypes = ['Other', 'Neutrophil', 'Basophil']
            for title in atype[2:]:
                if 'neutrophil' in title:
                    ahash[title] = 1
                elif 'basophil' in title:
                    ahash[title] = 2    
                else:
                    ahash[title] = 0            
        else:
            atypes = ['Other', 'Neutrophil']
            for title in atype[2:]:
                if 'neutrophil' in title:
                    ahash[title] = 1
                else:
                    ahash[title] = 0    
        return ahash, atypes
    
    def getJeffrey2006(self, tn=1):
        self.prepareData("MACV52")
        atype = self.h.getSurvName("c Title")
        ahash, atypes = self.getNeuAtypes(atype, tn)
        self.initData(atype, atypes, ahash)
        
    def getAllantaz2011(self, tn=1):
        self.prepareData("NEU6", cfile = "/Users/rohan/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName("c cell type (ch1)")
        ahash, atypes = self.getNeuAtypes(atype, tn)
        self.initData(atype, atypes, ahash) 
    
    def getNorvershtern2011(self):
        self.prepareData("G19")
        atype = self.h.getSurvName('c cell type')
        atypes = atype
        ahash = {atype[i]:i for i in range(len(atype))}
        self.initData(atype, atypes, ahash)         
        
    def getNorvershtern2011(self, tn=1):
        self.prepareData("G19")
        atype = self.h.getSurvName('c cell type')
        ahash, atypes = self.getNeuAtypes(atype, tn)
        self.initData(atype, atypes, ahash)  
    
    def getMonaco2017(self, tn=1):
        self.prepareData("NEU12", cfile = "/Users/rohan/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName("c cell type (ch1)")
        ahash, atypes = self.getNeuAtypes(atype, tn)
        self.initData(atype, atypes, ahash)
        
    def getSippel2018(self, tn=1):
        self.prepareData("MACV53")
        atype = self.h.getSurvName('c treatment')
        atypes = ['IL5', 'Vehicle', 'dkPGD2']
        ahash = {'IL5':0, 'Vehicle':1, 'dkPGD2':2}
        if (tn == 2):
            atypes = ['Normal', 'Type 2 Active']
            ahash = {'Vehicle':0, 'IL5':1, 'dkPGD2':1}
        self.initData(atype, atypes, ahash)  
        
    def getNelson2019(self, tn=1):
        self.prepareData("MACV63")
        atype = self.h.getSurvName('c treatment')
        atypes = ['GMCSF', 'IL3', 'IL5', 'NegCtrl']
        ahash = {'GMCSF':0, 'IL3':1, 'IL5':2, 'NegCtrl':3}
        if (tn == 2):
            atypes = ['Normal', 'Active']
            ahash = {'GMCSF':1, 'IL3':1, 'IL5':1, 'NegCtrl':0}
        else:
            atypes = ['Normal', 'Type 2 Active', 'Other Active']
            ahash = {'GMCSF':1, 'IL3':1, 'IL5':1, 'NegCtrl':0}           
        self.initData(atype, atypes, ahash)

    def plotViolinBar(ana, desc=None):
        fig = plt.figure(figsize=(4,4), dpi=100)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        ax1 = plt.subplot2grid((4, 1), (0, 0))
        ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
                  'genes': [], 'ax': ax1, 'acolor': acolor}
        ax = ana.printTitleBar(params)
        res = ana.getROCAUC()
        ax.text(len(ana.cval[0]), 4, res)
        if desc is not None:
            ax.text(-1, 2, desc, horizontalalignment='right',
                        verticalalignment='center')
        params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
                'genes': [], 'ax': ax2, 'acolor': acolor, 'vert': 0}
        ax = ana.printViolin(None, params)
        return fig

    def plotViolin2(data, atypes, params,geneName):
        fig = plt.figure(figsize=(4,4), dpi=100)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        ax1 = plt.subplot2grid((4, 1), (0, 0))
        ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        df = pd.DataFrame()
        df[geneName] = [k for i in range(len(data)) for k in data[i]]
        df['category'] = [atypes[i] for i in range(len(data)) for k in data[i]]
        m1 = []
        pvals = []
        for i in range(1, len(data)):
            if len(data[i]) <= 0:
                m1 += [0]
                pvals += [""]
                continue
            m1 += [max(data[i]) + (max(data[i]) - min(data[i])) * 0.1]
            t, p = ttest_ind(data[0],data[i], equal_var=False)
            if (p < 0.05):
                pvals += ["p=%.3g" % p]
            else:
                pvals += [""]
        dpi = 100
        if 'dpi' in params:
            dpi = params['dpi']
        w,h = (1.5 * len(atypes), 4)
        if 'w' in params:
            w = params['w']
        if 'h' in params:
            h = params['h']
        color_sch1 = acolor
        if 'acolor' in params:
            color_sch1 = params['acolor']
        sns.set()
        sns.set_style("white")
        sns.set_style({'text.color': '.5', 
            'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
        sns.set_context("notebook")
        sns.set_palette([adj_light(c, 1.5, 1) for c in color_sch1])
        ax = None
        if 'ax' in params:
            ax = params['ax']
        if ax is None:
            fig,ax = plt.subplots(figsize=(w,h), dpi=dpi)
        width = 1
        height = 1
        if 'width' in params:
            width = params['width']
        if 'vert' in params and params['vert'] == 1:
            ax = sns.violinplot(x="category", y=geneName, inner='quartile',
                    linewidth=0.5, scale="width", width=0.3, ax = ax, data=df,
                    order = atypes)
            ax = sns.swarmplot(x="category", y=geneName, color = 'blue', alpha=0.2,
                    ax=ax, data=df, order = atypes)
            ax.set_xlabel("")
            pos = range(len(atypes))
            for tick,label in zip(pos,ax.get_xticklabels()[1:]):
                #print(tick, label)
                ax.text(pos[tick], m1[tick - 1], pvals[tick - 1],
                        horizontalalignment='center', size=12,
                        color='0.3')
            ax.yaxis.grid(True, clip_on=False)
        else:
            ax = sns.violinplot(x=geneName, y="category", inner='quartile',
                    linewidth=0.5, width=width, ax = ax, data=df,
                    order = atypes)
            ax = sns.swarmplot(x=geneName, y="category", color = 'blue', alpha=0.2,
                    ax=ax, data=df, order = atypes)
            ax.set_ylabel("")
            pos = range(len(atypes))
            for tick,label in zip(pos[1:],ax.get_yticklabels()[1:]):
                ax.text(m1[tick - 1], pos[tick]-0.5, pvals[tick - 1],
                        horizontalalignment='center', size=12,
                        color='0.3')
            ax.xaxis.grid(True, clip_on=False)
        return fig

    def plotDensityBar(ana, desc=None):
        fig = plt.figure(figsize=(4,4), dpi=100)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        ax1 = plt.subplot2grid((4, 1), (0, 0))
        ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
                  'genes': [], 'ax': ax1, 'acolor': acolor}
        ax = ana.printTitleBar(params)
        res = ana.getMetrics(ana.cval[0])
        ax.text(len(ana.cval[0]), 4, ",".join(res))
        if desc is not None:
            ax.text(-1, 2, desc, horizontalalignment='right',
                        verticalalignment='center')
        ax = ana.densityPlot(ax2, acolor)
        return fig

    def processData(ana, l1, wt1, desc=None, violin=1):
        ana.orderData(l1, wt1)
        if (violin == 1):
            return plotViolinBar(ana, desc)
        return plotDensityBar(ana, desc)