# Import a few matplotlib things 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
# Import the lhapdf library 
import lhapdf

# We first choose the PDF set that we are going to use for the proton:
pset = lhapdf.getPDFSet("CT18Anlo") 
# Define a function to extract the proton pdf
#   - xval corresponds to Björken-x
#   - Q2val to the virtuality squared
#   - index is the number of the PDF set. The maximum is reported in the corresponding .info file
def proton_pdf(pdgid, xval, Q2val,index):
    proton = pset.mkPDF(index)
    return proton.xfxQ2(pdgid, xval, Q2val)

# We are going to plot just two Q^2 values and the gluon. Easy to generalize for more options
nuclei = ['Pb208']
max_index_p = pset.size
Q2values = [3,10]
pdgid = 21 # we plot just the gluon

# In principle we could have an array of nuclei
for nucleus in nuclei:

    Aset = lhapdf.getPDFSet("EPPS21nlo_CT18Anlo_"+nucleus)
    max_index_A = Aset.size
    # Define a function for the Pb nucleus
    def nuclear_pdf(pdgid, xval, Q2val,index):
        nucleus = Aset.mkPDF(index)
        return nucleus.xfxQ2(pdgid, xval, Q2val)
    
    # Now we do the plot for a given nucleus
    with PdfPages('R_'+nucleus+'.pdf') as pdf:
        fig = plt.figure(figsize=(5,3.8))
        gs = gridspec.GridSpec(1,1)
        ax = plt.gca()
        ax = plt.subplot(gs[0])
        ax.grid(True, lw=0.5, ls=':', zorder=0)
        ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
        # Labels and title
        ax.set_ylim(0,2)
        ax.set_xscale("log")
        ax.set_title('Nucleus='+nucleus+' , flavour='+str(pdgid))
        ax.set_xlabel(r'$x$',fontsize=16)
        ax.set_ylabel(rf'$R(x,Q^2)$',fontsize=16)

        # We define the x-values where we evaluate the PDFs
        xs = [x for x in np.logspace(-4, 0, 50)]
        # Loop over Q^2 values
        for Q2val,c in zip(Q2values,["crimson","navy"]):
            # To compute the uncertainty band we need to find the envelope accounting
            # for all PDF members
            # Initialize arrays
            xgp = np.zeros(pset.size)
            xgA = np.zeros(Aset.size)
            R_max_val = np.zeros(np.size(xs))
            R_min_val = np.zeros(np.size(xs))
            R_val = np.zeros(np.size(xs))
            # Loop over x values
            for x,ix in zip(xs,range(np.size(xs))):
                # all proton PDFs
                for imem in range(pset.size):
                    xgp[imem] = proton_pdf(pdgid,x,Q2val,imem)
                # all nucleus PDFs
                for imem in range(Aset.size):
                    xgA[imem] = nuclear_pdf(pdgid,x,Q2val,imem)
                
                # Find the uncertainties
                unc_p = pset.uncertainty(xgp)
                unc_A = Aset.uncertainty(xgA)

                # Upper and lower value of the (n)PDFs
                upper_value_p = unc_p.central + unc_p.errplus
                lower_value_p = unc_p.central - unc_p.errminus
                upper_value_A = unc_A.central + unc_A.errplus
                lower_value_A = unc_A.central - unc_A.errminus

                # We consider 4 combinations to define the envelope
                R_max_val[ix] = max(upper_value_A/upper_value_p, upper_value_A/lower_value_p, lower_value_A/upper_value_p, lower_value_A/lower_value_p)
                R_min_val[ix] = min(upper_value_A/upper_value_p, upper_value_A/lower_value_p, lower_value_A/upper_value_p, lower_value_A/lower_value_p)
                R_val[ix] =  unc_A.central/unc_p.central
            # Draw the central value     
            plt.plot(xs,R_val,color=c,label=rf'$Q^2={Q2val}$ GeV$^2$')
            # Draw the error band
            plt.fill_between(xs,R_min_val,R_max_val,color=c,alpha=0.2)
        ax.legend(loc='best')
        # Line for reference
        ax.hlines(1,0,1,color='gray',ls='--')
        pdf.savefig(bbox_inches='tight')
    plt.close()

