# import ROOT in batch mode
import sys,getopt

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
import math
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv


def makeEfficiencies(inputfiles,num,den,title,legends,outname):
    ROOT.gStyle.SetOptStat(0)
    
    ROOT.gROOT.Reset()
    #make canvas to save plots to
    c = ROOT.TCanvas('c')

    _files=[]
    for i in inputfiles:
        _files.append(ROOT.TFile(i))
    
    _effs=[]
    _num = ROOT.TH1F()
    _den = ROOT.TH1F()
    _num.SetDirectory(0)
    _den.SetDirectory(0)
    ROOT.TH1.AddDirectory(ROOT.kFALSE)   
    for f in _files:
        numname = "muonDebugger/%s" %(num)
        denname = "muonDebugger/%s" %(den)
        print "Trying to get %s and %s from %s" %(numname,denname, f.GetName())
        _num = f.Get(numname).Clone()
        _den = f.Get(denname).Clone()
        
        eff = ROOT.TEfficiency(_num,_den)
        _effs.append(eff)
        
    leg = ROOT.TLegend(0.30,0.15,0.70,0.3);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0)
    
    k=0
    for eff in _effs:
        if k==0 :
            eff.SetLineWidth(2)
            eff.SetLineColor(k+1)
            eff.SetTitle(title)
            eff.Draw()
            
            ## SET CORRECT AXIS VALUES
            ROOT.gPad.Update()
            graph = eff.GetPaintedGraph()
            graph.SetMinimum(0.0)
            graph.SetMaximum(1.02)
            ROOT.gPad.Update()
            
        else:
            eff.SetLineWidth(2)
            eff.SetLineColor(k+1)
            eff.SetTitle(title)
            eff.Draw("SAME")
        
        leg.AddEntry(eff,legends[k],"l")
        k+=1
    
    leg.Draw("SAME")
    
    c.SaveAs(outname+".C");
    c.SaveAs(outname+".pdf");
    c.SaveAs(outname+".png");
  
    
#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeEffs.py [options]",description="Extract Efficiencies from file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--inputfiles", type=str, help='The list of input files, comma separated if more than one file',required=True,nargs=1)
    parser.add_argument("-num",default="numerator",help='the numerator to be used',required=True)
    parser.add_argument("-den",default="denominator",help='the numerator to be used',required=True)
    parser.add_argument("-title",dest="title", type=str, help='title of the histo;xaxis;yaxis')    
    parser.add_argument("-leg", type=str,help='the list of legends to be used, comma separated',required=True,nargs=1)
    parser.add_argument("-outname",dest="outname", default="efficiency", help='name of the outputfile')
    
#    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')



    args = parser.parse_args()
    files = args.inputfiles[0].split(",")
    legends = args.leg[0].split(",")
    
    makeEfficiencies(files,args.num,args.den,args.title,legends,args.outname)
    print "DONE"


