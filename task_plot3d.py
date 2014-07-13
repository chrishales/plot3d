from taskinit import *
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.cm as cm
import numpy as np
from pylab import ion,ioff

# Copyright (c) 2014, Christopher A. Hales
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * The NAMES OF ITS CONTRIBUTORS may not be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL C. A. HALES BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# HISTORY:
#   1.0  12Jul2014  Initial version.
#

def plot3d(vis,fid,datacolumn,corr,plotall,spw,timecomp,chancomp,clipamp,outpng):

    #
    # Task plot3d
    #
    #    Quickly inspect data for RFI by plotting time vs frequency vs amplitude
    #    Christopher A. Hales
    #
    #    Version 1.0 (tested with CASA Version 4.2.1)
    #    12 July 2014
    
    casalog.origin('plot3d')
    
    # get number of baselines
    tb.open(vis+'/ANTENNA')
    atble=tb.getcol('NAME')
    tb.close
    nant=atble.shape[0]
    nbaselines=nant*(nant-1)/2
    
    # channel to frequency conversion
    tb.open(vis+'/SPECTRAL_WINDOW')
    vtble=tb.getcol('CHAN_FREQ')
    tb.close
    nspw=vtble.shape[1]
    
    # calculate number of effective channels per spw
    # I assume that the end channels of each spw have been flagged.
    # Force individual channels to remain at either end of spw,
    # in order to ensure amplitudes are zero in between
    # non-contiguous spw's. This will also ensure that it is
    # easier to see spw boundaries in between contiguous spw's.
    nchan = int(np.floor((vtble.shape[0]-2)/float(chancomp)))+2
    # guard against the user inputting infinite chancomp
    if nchan == 2:
        nchan = 3
    
    if plotall:
        # later we will gather scan information by looking at
	# a single spw and assuming it represents all spw's
	spw = 0
	
	# I don't make any effort to set the amplitude to
        # zero in the gaps between spw's (ie if spw's are not
        # contiguous) because I will assume that flagging of
        # spw edge channels has already taken place. Thus
        # there is no need to pad spw's with extra channels
        # if they happen to sit next to a gap in frequency
        # coverage. For a more general code this would not
        # be appropriate.
        N=np.zeros(nchan*nspw)
        t=0
        for i in range(nspw):
            # the following copied from single-spw "else" part below
            k=0
            
            # 1st channel in spw
            N[t] = vtble[k,i]/1e6
            t += 1
            k += 1
            
            # middle channels
            # check if we are in the last block
            while k+2*chancomp-1 <= vtble.shape[0]-2:
                for h in range(chancomp):
                    N[t] = N[t] + vtble[k+h,i]
                N[t] = N[t]/1e6/chancomp
                t += 1
                k += chancomp
            # for the last block, just combine everything remaining    
            for h in range(k,vtble.shape[0]-1):
                N[t] = N[t] + vtble[h,i]
            N[t] = N[t]/1e6/len(range(k,vtble.shape[0]-1))
            t += 1
            
            # last channel in spw
            N[t] = vtble[vtble.shape[0]-1,i]/1e6
            t += 1
        
        ## TESTING: get regular channel data to compare
        #Q=np.zeros([vtble.shape[0]*nspw])
        #t=0
        #for i in range(nspw):
        #    for k in range(vtble.shape[0]):
        #        Q[t] = vtble[k,i]/1e6
        #        t += 1
    else:
        N=np.zeros(nchan) 
        t=0
        k=0
        
        # 1st channel in spw
        N[t] = vtble[k,spw]/1e6
        t += 1
        k += 1
        
        # middle channels
        # check if we are in the last block
        while k+2*chancomp-1 <= vtble.shape[0]-2:
            for h in range(chancomp):
                N[t] = N[t] + vtble[k+h,spw]
            N[t] = N[t]/1e6/chancomp
            t += 1
            k += chancomp
        
        # for the last block, just combine everything remaining    
        for h in range(k,vtble.shape[0]-1):
            N[t] = N[t] + vtble[h,spw]
        N[t] = N[t]/1e6/len(range(k,vtble.shape[0]-1))
        t += 1
        
        # last channel in spw
        N[t] = vtble[vtble.shape[0]-1,spw]/1e6
        
        ## TESTING: get regular channel data to compare
        #Q=np.zeros(vtble.shape[0])
        #t=0
        #for k in range(vtble.shape[0]):
        #    Q[t] = vtble[k,spw]/1e6
        #    t += 1
    
    ms.open(vis)
    
    # assume time is same for each spw
    # this is not the most efficient place in the code for this bit, meh
    ms.reset()
    ms.msselect({'field':str(fid),'spw':str(spw)})
    
    # get the raw timestamps
    Z=ms.getdata('time')['time']
    
    # Get scan summary in prep for calculating time steps.
    # Note that CASA currently reports all spw's in the
    # scan summary, rather than the 1 selected above. meh
    scan_summary = ms.getscansummary()
    scan_list = []
    for scan in scan_summary:
        if scan_summary[scan]['0']['FieldId'] == fid:
            scan_list.append(int(scan))
    
    scan_list.sort()
    
    # get integration time; assume it doesn't change in any
    # way throughout the observation, ie between spw's, etc
    inttime=scan_summary[str(scan_list[0])]['0']['IntegrationTime']
    
    # Calculate number of true time steps per scan.
    # In the code below, a dummy timestep will be added at each
    # end of each scan to ensure amplitudes are zero in between
    # non-contiguous scans. This will also ensure that it is
    # easier to see scan boundaries in between contiguous
    # scans. The 1st and last timestamp do not contribute to
    # the time compression stuff.
    # Also calculate effective time steps per scan, so that
    # I can call the variable effntime...!
    scan_ntime    = []
    scan_effntime = []
    for scan in scan_list:
        scan_ntime.append(scan_summary[str(scan)]['0']['nRow']/nbaselines/nspw)
        tempvar=int(np.floor(scan_summary[str(scan)]['0']['nRow']/nbaselines/nspw/float(timecomp)))+2
        # guard against the user inputting infinite timecomp
        if tempvar == 2:
            scan_effntime.append(tempvar+1)
        else:
            scan_effntime.append(tempvar)
    
    ntime = sum(scan_effntime)
    
    # go through each scan and add a dummy timestep before
    # and after each one, with time difference equal to
    # one ten thousandth of the integration time (make this
    # small so that a slope doesn't show up in the plot)
    intdividefactor=10000.0
    M=np.zeros(ntime)
    t=0
    for d in range(len(scan_list)):
        checkfirst=True
        k=0
        while k+2*timecomp-1 <= scan_ntime[d]-1:
            for h in range(timecomp):
                if checkfirst:
                    t+=1
                
                M[t] += Z[(sum(scan_ntime[:d])+k+h)*nbaselines]
                if checkfirst:
                    M[t-1] = M[t]-inttime/intdividefactor
                    checkfirst=False
            
            M[t] = M[t]/timecomp
            t += 1
            k += timecomp
        
        for h in range(scan_ntime[d]-k):
            if checkfirst:
                t+=1
            
            M[t] += Z[(sum(scan_ntime[:d])+k+h)*nbaselines]
            if checkfirst:
                M[t-1] = M[t]-inttime/intdividefactor
                checkfirst=False
        
        M[t] = M[t]/len(range(scan_ntime[d]-k))
        t+=1
        M[t] = M[t-1]+inttime/intdividefactor
        t+=1
        
        # TESTING: get regular time data to compare
        # will need to manually increment by nbaselines. meh
        # Z[nbaselines*1]
        # compare before subtracting to make comparison easier
    
    # time is in seconds from zero modified Julian date...not very aesthetic
    # subtract off the starting time and convert to minutes
    M=(M-M.min())/60
    
    # go through each spectral window and extract amplitude data
    if plotall:
        for i in range(nspw):
            ms.reset()
	    ms.msselect({'field':str(fid),'spw':str(i)})
	    
            # visibility data (X,Y,Z) where
            #   X=4 (RR,RL,LR,LL)
            #   Y=number of channels
            #   Z=number of rows (visibilities/4)
            tempdata=ms.getdata(datacolumn)
            tempflag=ms.getdata('flag')
            # true flag means I should flag it, so switch to ensure good points remain
            tempflag=np.invert(tempflag['flag'][corr])
            
            # select amplitude data associated with requested correlation
            # and ensure any existing flagged points are set to zero
            P1 = np.multiply(abs(tempdata[datacolumn][corr]),tempflag)
            
            # time + baseline compression
            P2=np.zeros([P1.shape[0],ntime])
            # loop over channels
            # yes, this is inefficient, but hopefully easier to understand
            for s in range(P1.shape[0]):
                t=0
                for d in range(len(scan_list)):
                    checkfirst=True
                    k=0
                    while k+2*timecomp-1 <= scan_ntime[d]-1:
                        if checkfirst:
                            t+=1
                        
                        P2[s,t] = max(P1[s,(sum(scan_ntime[:d])+k)*nbaselines:(sum(scan_ntime[:d])+k+timecomp)*nbaselines])
                        if clipamp>=0:
                            P2[s,t] = min(clipamp,P2[s,t])
                        
                        if checkfirst:
                            P2[s,t-1] = 0.0
                            checkfirst=False
                        
                        t += 1
                        k += timecomp
                    
                    if checkfirst:
                        t+=1
                    
                    tempvar=len(range(scan_ntime[d]-k))
                    P2[s,t] = max(P1[s,(sum(scan_ntime[:d])+k)*nbaselines:(sum(scan_ntime[:d])+k+tempvar)*nbaselines])
                    if clipamp>=0:
                        P2[s,t] = min(clipamp,P2[s,t])
                    
                    if checkfirst:
                        P2[s,t-1] = 0.0
                        checkfirst=False
                    
                    t+=1
                    P2[s,t] = 0.0
                    t+=1
            
            # channel compression
            # for clarity, don't combine this step with the
            # time+baseline compression above
            P3=np.zeros([nchan,ntime])
            
            # 1st channel in spw
            t=0
            k=0
            P3[t] = P2[t]
            t += 1
            k += 1
            
            # middle channels
            while k+2*chancomp-1 <= P2.shape[0]-2:
                for h in range(chancomp):
                    P3[t] = np.maximum(P3[t],P2[k+h])
                t += 1
                k += chancomp
            for h in range(k,P2.shape[0]-1):
                P3[t] = np.maximum(P3[t],P2[h])
            t += 1
            
            # last channel in spw
            P3[t] = P2[P2.shape[0]-1]
            
            if i == 0:
                P=P3
            else:
                P=np.concatenate((P,P3),axis=0)
            
            # not needed because of selection above
            # spectral window, with same number of rows as Z above
            #sdata=ms.getdata('data_desc_id')
            
            # not needed because of selection above
            # field ID
            #fdata=ms.getdata('field_id')
    else:
        # just copy the important steps from above
        ms.reset()
	ms.msselect({'field':str(fid),'spw':str(spw)})
	tempdata=ms.getdata(datacolumn)
        tempflag=ms.getdata('flag')
        tempflag=np.invert(tempflag['flag'][corr])
        P1=np.multiply(abs(tempdata[datacolumn][corr]),tempflag)
        
        # time + baseline compression
        P2=np.zeros([P1.shape[0],ntime])
        # loop over channels
        # yes, this is inefficient, but hopefully easier to understand
        for s in range(P1.shape[0]):
            t=0
            for d in range(len(scan_list)):
                checkfirst=True
                k=0
                while k+2*timecomp-1 <= scan_ntime[d]-1:
                    if checkfirst:
                        t+=1
                    
                    P2[s,t] = max(P1[s,(sum(scan_ntime[:d])+k)*nbaselines:(sum(scan_ntime[:d])+k+timecomp)*nbaselines])
                    if clipamp>=0:
                        P2[s,t] = min(clipamp,P2[s,t])
                    
                    if checkfirst:
                        P2[s,t-1] = 0.0
                        checkfirst=False
                    
                    t += 1
                    k += timecomp
                
                if checkfirst:
                    t+=1
                
                tempvar=len(range(scan_ntime[d]-k))
                P2[s,t] = max(P1[s,(sum(scan_ntime[:d])+k)*nbaselines:(sum(scan_ntime[:d])+k+tempvar)*nbaselines])
                if clipamp>=0:
                    P2[s,t] = min(clipamp,P2[s,t])
                
                if checkfirst:
                    P2[s,t-1] = 0.0
                    checkfirst=False
                
                t+=1
                P2[s,t] = 0.0
                t+=1
        
        # channel compression
        # for clarity, don't combine this step with the
        # time+baseline compression above
        P=np.zeros([nchan,ntime])
        
        # 1st channel in spw
        t=0
        k=0
        P[t] = P2[t]
        t += 1
        k += 1
        
        # middle channels
        while k+2*chancomp-1 <= P2.shape[0]-2:
            for h in range(chancomp):
                P[t] = np.maximum(P[t],P2[k+h])
            t += 1
            k += chancomp
        for h in range(k,P2.shape[0]-1):
            P[t] = np.maximum(P[t],P2[h])
        t += 1
        
        # last channel in spw
        P[t] = P2[P2.shape[0]-1]
    
    ms.close()
    
    # clear memory, not needed any more
    vtble=[]
    Z=[]
    P1=[]
    P2=[]
    P3=[]
    tempdata=[]
    tempflag=[]
    
    # M=time, N=frequency , P=amplitude
    M2D,N2D=np.meshgrid(M,N)
    ion()
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('time (mins)')
    ax.set_ylabel('frequency (MHz)')
    ax.set_zlabel('amplitude')
    if plotall:
        plot_title='column='+str(datacolumn)+' fid='+str(fid)+' corr='+str(corr)
        figname=vis.strip('.ms')+'_plot3d_fid'+str(fid)+'_corr'+str(corr)+'_t'+str(timecomp)+'_c'+str(chancomp)
    else:
        plot_title='column='+str(datacolumn)+' fid='+str(fid)+' corr='+str(corr)+' spw='+str(spw)
        figname=vis.strip('.ms')+'_plot3d_fid'+str(fid)+'_corr'+str(corr)+'_spw'+str(spw)+'_t'+str(timecomp)+'_c'+str(chancomp)
    
    ax.set_title(plot_title)
    #ax.set_zscale('log')
    ax.plot_surface(M2D, N2D, P, rstride=1, cstride=1, cmap=cm.jet)
    #if isinstance(plotfig,str):
    #    figname=plotfig
    #    plotfig=1
    
    if outpng:
        fig.savefig(figname)
    
    ioff()
