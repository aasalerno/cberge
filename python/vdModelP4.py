# % ### vdModelP4.py ###     03.15.16  C. Bergevin
'''
# # Edited version of C. Bergevin's code written to run in Python 2.7.12 

# % -----------------------------------------------------------------------

# % **********************
# % [original A. Salerno; re-hacked by C. Bergevin]
# % [NOTE: also requires modelP#.m and visualizeVDP#.m]
# % This code is built to solve the a modified version of Vilfan and Duke (BJ 2008)
# % using the global "panp.pilla" coupling of Bergevin & Shera (JASA 2010). In short:
# % - array of coupled oscillators w/ some form of user-specified CF distribution (P.CFtype)
# % - oscillators can be active (P['e']_base>0) and/or nonlinear (P.B_base)
# % - coupling can have two forms: nearest-neighbor (P.dI_base and P.dR_base) and/or "global" via the
# % panp.pilla (P.alpha and P.beta)
# % - solver (P.solver) uses either a straight RK4 ingegrator (via ode4.m) w/ fixed step-size
# % or Matlab's ode45 (adaptive step-size; interpolates final time series to
# % get the spacing you want)

# % **********************
# % v4 UPDATES
# % - (03.15.16) options to "chunkify" the tonotonp.pic map so to make more
# % consistent w/ Anolis morphology (assuming that each of the 4-5 hair cells at a given
# % longitudinal location have roughly the same CF)
# % - (03.15.16) cleaned up the organization a bit
# % - (04.01.16) option to "freeze" bundle roughness now implemented


# % **********************
# % TO DO (coming off 2016 ARO)
# % ** coding-wise **
# % - b. add option for tone burst stimuli so to compare to TB-SOAE data
# % - c. set up a SFOAE-like suppression framework, where the response is the difference
# % between either a high and (scaled-up) low level (same freq.) or 2TS
# % - d. in the analysis code, setup timecourse plot of all oscillators (as
# % Anthony did for ASA 2015 slide 16 --> extract some meaningful
# % quantification re "waves"?
# % - e. inclusion of the middle ear (see B&S10 for how things were handled
# % there)

# % ** analysis-wise **
# % - 1. check the "coherence" calculation and verify that the 'tightness' of
# % a cluster can be qualified/quantified via such
# % - 2. examine timecourse of oscillators not in coherent (i.e., decoherent?
# % weakly coherent?) cluster --> multi-freq. at any inst. of 'hopnp.ping' in
# % time?
# % - 3. w/ frozen roughness, check how cluster coherence is affected for a
# % high-ish level tone at several diff. freqs. --> is a freq.-dependent
# % breakup apparent?
# % - 5. can responses to high-level tones be quantified via a cohrence measure
# %   (--> presumably a 'decohered' cluster/peak would see it's coherence
# %   value lowered)
# % - 6. check height of coherent vs decoherent clusters to see what the
# % reduction is --> better yet, develop a quantitative measure of energy (re
# % the peak/clusters)
# % - 7. What if the underlying map is not smoothly tonotonp.pic? THat is, not
# % just some rough exponential.  What variations could we add to the
# % tonotonp.pic map so to make the output here more SOAE-like (e.g., CFs are
# % clustered already?)

# % **********************
# % Some (roughly) answered questions
# % - Q. examine response to very high-level tones --> says something re
# % "dynamic range" of model?
# % ANS: P.A=50 is too high (i.e., everything moves only at tone) whereas
# % clusters far away still apparent for P.A=10; so in this sense, the
# % dynamic range is limited....
# % - Q. systematically explore effects of nearest-neighbor (NN) vs panp.pilla
# % (Global) coupling to compare/contrast --> for example, if NN-coupling is
# % off, can Global coupling allow for clustering?
# % ANS: It seems like I can get some form of clustering with predominantly
# % panp.pilla coupling. Start w/ ARO base params ({#}) and change the following:
# % P.kap_base=0.2, P.dRP= -1.0, P.dIP= 1.5
# % ANS2: If P.kap_base= 0.2, no clustering if P.alpha= P.beta=0. But if they
# % are made non-zero, clustering. So Global coupling can be sufficient to
# % cause clustering in cases of weak NN coupling. Interestingly, seems like this may need
# % to be resistive in nature (rather than purely reactive)
# % ANS3: If the oscillator CFs are "chunkified", clustering can still happen at higher
# % freqs. w/ panp.pilla coupling only (i.e., NN set to zero). However having
# % both forms of coupling (P.kap_base= 1.0) results in too few clusters
# % (higher freqs. are weakly entrained to lowest freqs.)
# % - Q. With bundle roughness "frozen", are resulting spectra identical
# % across runs?
# % ANS: Almost, but not quite. Some minor differences are apparent, presumably
# % related to the randomized initial conditions


# % **********************
# % ===========
# % [as background, the "unmodified" V&D08 eqns. of motion are described as...]
# % The general eqn. of motion (see eqn.11 of V&D 2008; eqn.1 in Wit & van Dijk 2012) is:
# % z'(j) = (i*w(j) + e(j))*z(j) + kap(j)(dR(j) + i*dI(j))*(z(j+1)+z(j-1)-2*z(j)) + L - B*abs(z(j))^2*z(j)
# % - z is the complex state variable [Re= position (x), Im= (scaled) velocity (v)]
# % - L is the non-autonomous external driving term
# % [Note: args[' listed below are static, but vary as a function of j]
# % - j is the oscillator # (i is sqrt(-1))
# % - w is the natural/characteristic frequency (of the j'th oscillator) (=CF)
# % - e is the damnp.ping coefficient (of the j'th oscillator)
# % - dR and dI are the resistive and reactive coupling terms respectively (of the j'th oscillator)
# % - B is the coefficient for the nonlinear (cubic) term

# % ===========
# % NOTES:
# % o Output of this code can then be visualized with visualizeVDP#.m
# % o This code uses the subfunction ./modelP1.m or ./modelP1mod.m (where the eqns. of motion are specified)
# % o User specifies solver: either ode4 or ode45 (see note above)
# % o On 2015 Powerbook w/ default V&D08 params, ode4 takes 53 s to solve and ode45 takes 30 s
# % o P.CFtype (i.e., linear vs log map) has a big effect...
# % o sans global coupling and w/ "default" args[' ([] values in comments), to 1st order,
# %   simulations appear consistent with Wit et al. JASA 2012
# % o oscillator #1 is the panp.pilla, with the "bundles" being oscillators 2-P['N']
# % o Since first oscillator is panp.pilla and that such is not accounted
# % for when defining the active oscillator parameters (e.g., P.w), there
# % will be a slight offset in the initial values (i.e., the first active
# % oscillator will have a CF>P.frange(1))
# % o Initial conditions (ICs) are either uniform [e.g., all oscillators at
# % x=1 (z=0 is the equlibrium) with zero velocity, or noisily distributed]
# % o P['e']FORM allows for several options to handle the active term (assumes P['e']_base>0)
# % - 0: uniform distribution of active oscillators (though possibly noisy if P['e']_noiseF~=0)
# % - 1: flips sign (i.e., makes passive), and makes a smaller subset (user-specified below) active
# % - 2: Same as Fig.2B from Wit et al. 2012 (i.e., tamped down to zero at both ends)
# % o Tried replacing structure P with better substructure hierarchy (e.g.,
# % P.bundlestuff --> P.b.bundlestuff) to improve "bookeenp.ping", but such
# % greatly slowed down ode45 for some reason...

# % ===========
# % To do (re 12/12/14)
# % --> are there 'waves' propagating along the 'chain'?
# % - test w/ basic 'reality check' conditions
# % - optimal value for P.stepSize?? (1/128 too small?)

# % -----------------------------------------------------------------------
'''
# Here we're going to define some funcitons that we're going to need
def quantFREQ(f,SR,N):
    df = SR/N
    return np.ceil(f/df)*df

def VDModel(t,z,args):
    dR= args['dR'];
    dI= args['dI'];
    e= args['e']
    B= args['B'];
    N= args['N'];
    w= args['w'];
    stepSize= args['stepSize'];
    L0= args['L'];
    kap= args['kap'];
    w_p = args['w_p']; #% Papillar frequency
    I_p = args['I_p']; #% Moment of inertia
    R = args['R']; #% Stiffness
    A_p = args['A_p'];
    h = args['h'];
    r_p = args['r_p'];
    pv = args['pv'];
    pt = args['pt'];
    K_p = args['K_p'];
    kval = np.ones(N-1)*K_p/N;
    dz1 = np.zeros(N); #% a column vecor - for the solution
    ti = np.floor(t/stepSize);
    if ti >= L.size:
        ti = L.size-1;

    dz1 = np.zeros(N)
    for j in range(N):
        L1 = L0(ti)

        if j==0:
            dz1[j] = (1j*w[j] + e[j])*z[j] + kap[j]*(dR[j] + 1j*dI[j])*(z[j+1]-z[j]) - B[j]*abs(z[j])**2*z[j] + L1 
        elif j == N-2:
            dz1[j] = (1j*w[j] + e[j])*z[j] + kap[j]*(dR[j] + 1j*dI[j])*(z[j-1]-z[j]) - B[j]*abs(z[j])**2*z[j] + L1
        elif j == N-1:
            dz1[j] = (-w_p/I_p)*(1+1j*R/w_p)*z[j].imag + 1j*K_p/(I_p*w_p)*z[j].real \
                   + 1j*h**2/(I_p*w_p)*np.sum(kval[:N]*z[N].real+np.sum(z[:N].real))\
                   - 1j*h*r_p*A_p/(I_p*w_p)*(pv-pt);
        else:
            dz1[j] = (1j*w[j] + e[j])*z[j] + kap[j]*(dR[j] + 1j*dI[j])*(z[j+1]+z[j-1]-2*z[j]) - B[j]*abs(z[j])**2*z[j] + L1

    return dz1

from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import ode
import matplotlib.pyplot as plt
plt.ion() # More "MATLAB" style plots (interactive)
# Also done to work with WSL (Windows Subsystem for Linux)

# for reproducability, include set seed
np.random.seed(123)

# For now, we won't worry about the counter...
# global nR;   % dummy variable for ODE solver update counter
# % ---------------------------------------
# % User Parameters (P is a structure containing these values)
# % NOTE: {##} indicates default 2016 ARO value
# % +++

# To make the shift as simple as possible, we will use dictionaries as opposed to structures...
P = {} # Initialize the dictionary
P['N'] = 100; # number of active oscillators (really P['N']-1) [80] {100}
P['CFtype']= 1;  # tonotonp.pic map type: 0-linear, 1-log (determines w) [0] {1}
P['frange']= np.array([1, 4.5]);    # min/max freqs. for tonotonp.pic map (determines w) [1 5] {[1 4.5]}
P['CFnoise']= 1;   # add noise to CF distribution? (determines w) [0] {1}
P['CFnoiseFact']= 2*10**(-2);    # factor for noisiness of CF distrib. [10^(-2)] {2*10^(-2)}
P['chunkify']= 0; # boolean re making the tonotonp.pic map "chunky" or not 0=no, 1=yes [0]
P['rows']= 29; # if P.chunkify=1, approx. how many groups? (since it can be noisy)
P['rowNoise']= 0.2; # adds a bit of variability as to the integer # of oscillators per CF
# % +++

# % Nearest-neighbor coupling properties (passive)
P['dR_base']= 0.15; # dR (dissipative coupling) base value {positive = dissipative} [0.15] {0.15}
P['dR_noiseF']= 0.0;  # dR irregularity factor (% re base) [0; though try small values like 0.05] {0}
P['dI_base']= -1.0; # d (reactive coupling) base value  [-1] {-1}
P['dI_noiseF']= 0.05;  # d irregularity factor (% re base) [0;...] {0.05}
P['kap_base']= 1;       # kap, coupling term [1] (a la Wit & van Dijk 2012) {1}
# % +++ 

# % "Active" properties for bundles
P['eFORM']= 0;      # form of 'active' term (see comments above) [0] {0}
P['e_base']= 1; #  e (individual oscillator damnp.ping) base value [1] (positive=active,negative=passive) {1}
P['e_noiseF']= 0.05;  # e irregularity factor (% re base) [0;...] {0.05}
P['B_base']= 1;    # B (strength of {cubic} nonlinearity) base value [1] {1}
P['B_noiseF']= 0.05;  # B irregularity factor (% re base) [0;...] {0.05}
# % +++

# % Panp.pilla parameters (P.alpha=P.beta=0 means no panp.pilla coupling)
P['w_p']= 2;      # panp.pilla CF [kHz] [2] {2}
P['eP']= -1;    # panp.pilla damnp.ping {positive=active,negative=passive} [-1] {-1}
# % +++

# % Panp.pilla coupling properties (passive)
P['dRP']= -0.15;     # (individual) viscous coupling factor between panp.pilla and bundles {negative val= positive damnp.ping??} [?] {-0.15}
P['dIP']= 1;     # (individual) elastic coupling factor between panp.pilla and bundles [?] {1}
P['alpha']= 1;   # (total) viscous coupling factor between panp.pilla and bundles [0] {1}
P['beta']= 1;    # (total) elastic coupling factor between panp.pilla and bundles [0] {1}
# % +++

# % Initial conditions
P['ICnoise']= 1;  # add noise to initial condition (IC) distribution? [1] {1}
P['ICbase']= 1;  # amplitude factor for IC distrib. [1] {1}
P['ICnoiseFact']= 1;  # amplitude factor for IC noisiness (reqs. P.ICnoise=1) [1] {1}
# % +++

# % External drive (i.e., non-autonomous term)
P['extT']= 0;   # 0-none, 1-single tone, 2-noise (uniform to all bundles) [0] {0,1}
rampDown = 0;
P['extD']= 1;   # 0-drives panp.pilla only, 1-drives panp.pilla and all bundles [0] {1}
P['A'] = 5.2; # amplitude  (Note: 'dynamic range' of model is relatively small [try 0.1-5] {5.2,...}
P['fe'] = 3.4; # tone freq. (requires extT=1) {3.4,...}
# % ----

# % Integration properties
P['lengthT']= 330;    # length ('seconds') to run integration [660] {a bit kludgy; stems from Wit et al. 2012?} {330}
P['stepSize'] = 1/128; # integration step-size (1/'seconds') [1/128??] {1/128}
# % ----

# % solver type: 0 - fixed-step RK4 (ode4), 1 - Matlab's ode45 (adaptive step-size)
P['solver']= 1; # {1}
# % +++

# % "Freeze" the roughness?? (only bundles have roughness, currently via: P.w, P['e'], P.dR, P.dI, P.B, P.kap)
# % Note: USE CAUTION (args[' other than those noted above are used as specified here, so ensure self-consist.)
# % NOTE: ICs are not frozen, nor is any dynamic noise (i.e., Brownian drive term)
P['freeze']= 0;
P['freezeFL']= 'yy.npy';  # name of parameter file to load (reqs. P.freeze=1)
P['freezeFS']= 'yy.npy';  # filename to save params too (reqs. P.freeze=1)
# % ----
# % Save data to file?
P['saveit'] = 1; # save data structure to file?? {1}
P['saveName'] = 'test';  # filename to save to (req. P.saveit=1)
# % ---------------------------------------

# % -------------------------------------------------------------------------
# % minor bookeenp.ping

# % +++
# % ensure the number of oscillators is even (for simplicity in further calculations) (why??)
if P['N']%2: P['N'] += 1

# % -------------------------------------------------------------------------
# % deal w/ various parameters

# % =====
# % deal w/ potentially "freezing" the roughness
# % NOTE: relevant variables - P.w, P['e'], P.dR, P.dI, P.B, P.kap

# % ^^^^^
# % no freeze (i.e., create new set of bundles props. w/ args[' specified above)
if P['freeze']==0:
    print('Generating fresh roughness pattern for bundle properties \n');
    # ++++ build parameter arrays, allowing for (normally-distributed) irregularity if specified
    P['dR']= P['dR_base']*np.ones([P['N']])*(1+P['dR_noiseF']*np.random.np.random.randn(P['N'])); # dissipative coupling
    P['dI']= P['dI_base']*np.ones([P['N']])*(1+P['dI_noiseF']*np.random.np.random.randn(P['N'])); # reactive coupling
    P['B']= P['B_base']*np.ones([P['N']])*(1+P['B_noiseF']*np.random.np.random.randn(P['N']));    # strength of (cubic) nonlinearity
    P['kap']= P['kap_base']*np.ones(P['N']);  # coupling term (a la Wit & van Dijk 2012)

    P['e'] = np.zeros([P['N']])
    # % ++++ deal w/ 'active' term (i.e., individual oscillator damnp.ping) and various options
    if P['eFORM']==0:  # uniform (possibly noisy)
        P['e']= P['e_base']*np.ones([P['N']])*(1+P['e_noiseF']*np.random.randn(P['N']));
    elif P['eFORM']==1:   # all are passive, save for a central subset which are active
        P['e_subset']= np.array([29, 35]);   # ** specify subset of oscillators to 'flip'
        P['e']= -P['e_base']*np.ones([P['N']])*(1+P['e_noiseF']*np.random.randn(P['N'])); 
        P['e'][P['e_subset'][0]:P['e_subset'][-1]]= -P['e'][P['e_subset'][0]:P['e_subset'][-1]]; # flip sign
    elif P['eFORM']==2:   #% Wit et al. 2012 Fig.2B (i.e., tamped at ends)
        if (P['N']!=80): 
            raise NameError('Wit et al. assumed N=80');
        P['e']= P['e_base']*np.ones([P['N']])*(1+P['e_noiseF']*np.random.randn(P['N']));
        P['e'][:10]= np.linspace(0,P['e_base'],10);  
        P['e'][-10:]= np.linspace(P['e_base'],0,10);

    # % ++++ Determine angular frequency (CF) distribution
    if P['CFtype']==0:
        P['w']= np.linspace(P['frange'][0],P['frange'][1],P['N'])*2*np.pi
    elif (P['CFtype']==1):
        P['w']= np.logspace(np.log10(P['frange'][0]),np.log10(P['frange'][1]),P['N'])*2*np.pi


    # % ++++ chunkify the tonotonp.pic map
    if (P['chunkify']==1):
        P['wc'] = P['w'].copy() # python requires prealloaction
        nn=0; mm= 0;    #indexers
        while nn < P['N']:
            chunks[mm]= np.floor(P['N']/P['rows']+ P['rowNoise']*np.random.randn(1));  # # of oscillators for this CF
            if (nn + 1 + chunks[mm])>=P['N']: # +1 because of 0-based counting in python
                break     #% if near the end, break out of while loop
            else:
                P['wc'][nn:nn+chunks[mm]+1]= np.mean(P['w'][nn:nn+chunks[mm]+1]);  # avg. the CF for that "chunk"
            nn += chunks(mm); 
            mm += 1;


        # take care of any final stragglers at the end
        if (P['wc'].size<=P['N']):
            P['wc'][P['wc'].size:P['N']]= np.mean(P['w'][P['wc'].size:P['N']])

        P['w']= P['wc'];  # rename (kludge?)

    # % ++++ Add noise to oscillator CFs?
    if P['CFnoise']==1:
        dw = (8*np.random.randn(P['N']) - 4)*P['CFnoiseFact']
        P['w'] += dw

elif (P['freeze']==1):
    print(['Frozen roughness: Loading in previously-saved bundle properties \n'])
    Fp= np.load(P['freezeFL']); # if you want to load in the .mat file, use scipy.loadmat -- however, .npy files will load faster
    P['w'] = Fp['P']['w']; 
    P['e'] = Fp['P']['e'];
    P['dR'] = Fp['P']['dR']; 
    P['dI'] = Fp['P']['dI'];
    P['B'] = Fp['P']['B']; 
    P['kap'] = Fp['P']['kap'];


# % -------------------------------------------------------------------------
# % ++++ IC (initial conditions)
if (P['ICnoise']==0):
    IC= P['ICbase']*np.ones([P['N']]);   # % uniform ICs
else:
    IC = P['ICnoiseFact']* np.sin(2*np.pi*np.random.randn(P['N'])); #% noisy IC (semi-bounded)

# % +++ Array of time values
P['T']= np.arange(0,P['lengthT'],P['stepSize'])

# % +++ External drive term
P['L'] = np.zeros(P['T'].size)
if P['extT']==0 : # % none
    None # Just as a holder, we use the zeros case as a method to preallocate
elif P['extT']==1:  #% tone
    P['fe'] = quantFREQ(P['fe'],P['stepSize']**(-1),P['T'].size);
    tRamp = 2; #%Ramnp.ping time in seconds
    nRamp = tRamp/P['stepSize'];
    P['L'][:nRamp] = P['A']*0.5*(1-np.cos(2*np.pi/(2*tRamp)*P['T'][:nRamp]))*np.sin(2*np.pi*P['fe']*P['T'][:nRamp]);
    if rampDown == 1:
        P['L'][nRamp:-nRamp] = P['A']*np.sin(2*np.pi*P['fe']*P['T'][nRamp:-nRamp]);
        P['L'][-nRamp:]= P['A']*0.5*(1-np.cos(2*np.pi/(2*tRamp)*P['T'][:nRamp]+np.pi))*np.sin(2*np.pi*P['fe']*P['T'][-nRamp:]);
    else:
        P['L'][nRamp:] = P['A']*np.sin(2*np.pi*P['fe']*P['T'][nRamp:]);
elif P['extT']==2:  #% noise
    P['L'] = P['A']*np.random.randn(P['T'].size);


# % -------------------------------------------------------------------------
nR= 1;


# % *** run the integration routine ***
# if P['solver']==0:

#     Z=ode4(@modelP3,P['T'],IC,P);   % explicit 4th order RK
# elseif P.solver==1
#     %[t Z] =ode45('modelP1mod',P['T'],IC,[],P);
#     %[T,Z]=ode45(@(t,z) modelP1v45(t,z,P),P['T'],IC); % note that ode45 solves, then interpolated back to time-spacing re P['T']
#     [T,Z]=ode45(@(t,z) modelP3(t,z,P),P['T'],IC); % note that ode45 solves, then interpolates back to time-spacing re P['T']
# end
# % -------------------------------------------------------------------------
# % package up various quantities as outputs
# info.params= P;
# data.info= info;
# data.IC= IC;
# data.Z = Z; % main output

# % ---
# % save data structure
# if (P.saveit==1), fprintf(['Saving waveform data and overwriting (',P.saveName,'.mat) \n']); save(P.saveName,'data'); end
 
# % ---
# % save model args[' to file (this should be turned "ON" by default in case
# % there is a decent spectra
# if 1==1
#     input(['Saving args[' and overwriting ',P.freezeFS,' (hit Cntl+C to cancel) \n']);
#     save(P.freezeFS,'P');
# end


# % ---
# clearvars -except data



# % **********************
# % v2 UPDATES
# % - (02.05.16) now can use ode45 as a solver option; updated overview/comments
# % - (02.10.16) now uses modelP3.m which allows user to specify what drive
# % term *directly* affects (via P['e']xtD) and also handles both ode4.m and
# % ode45.m solver routines

