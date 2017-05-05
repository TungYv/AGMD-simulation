
This notebook is used for AGMD programming. Attention! For membrane's
porosity needs to change both in Global and in Class.

For Alsaadi's 2013, with a spacer --> Change: 1.func.CalculateH 2.spacer
properties.

.. code:: python

    import numpy as np
    from matplotlib import pyplot
    % matplotlib inline

.. code:: python

    # Universal constant
    Mv = 18 #(kg/kmol) water 18kg/kmole
    R = 8.314*1000 #(J/(kmol K)) kmol, J = Pa * m3
    G = 9.80 # (m/s2) gravity
    Ptotal = 101325 # (Pa) Total pressure in the model, p = 101325Pa.
    ## For spacer
    SPACER_df = 0.4e-3 # (m) filament diameter
    SPACER_hs = 0.8e-3 #(m) thickness of spacer
    SPACER_epsilon = .5 # 77.85% spacer porosity from [28]. unsure...
    
    ## Each of the following variables is considered as Constant with slight deflection from its present value.
    #air gap
    Kab = 2.90e-2*3.6 #kJ/(m K hr) = kW/(m K)/3600 thermal conductivity of air mixture in air gap. From "heat transer"--dry air's property.
    
    # For hot side
    rho_hotwater = 980. #(kg/m3)
    Cp_hot = 4.185 #(kJ/kg/C)
    
    # For falling liquid film
    rho_l = 993.5 #(kg/m3) density of falling liquid water. From 'heat transfer'
    rho_av = 1.05 # (kg/m3) density of gas mixture in air gap. From both 'heat transfer and thermal dynamics'
    
    # For cold side
    rho_coldwater = 993. #(kg/m3), From "heat transer"
    Cp_cold = 4.179 # (kJ/kg/C)
    
    ## for the membrane
    Kair = 2.93e-2*3.6 #kJ/(m K hr) = kW/(m K)/3600 From 'Heat Transfer' dry air property
    Km = 0.29*3.6 #kJ/(m K hr) = kW/(m K)/3600 From A.G Fane.2003, at 75C.
    epsilon = 0.8 # porosity of membrane
    K = epsilon*Kair + (1-epsilon)*Km

.. code:: python

    # Compute C
    def CalculateC(d,delta_m,tau,T, P, Pv):
        '''Params
        -----------
        Some are global variables
        d: pore size of membrane(m)
        delta_m: thickness of membrane(m)
        tau: torousity of membrane
        T: Average temperature in membrane.(C)
        P: Total pressure in membrane.(Pa)
        Pv: Vapour pressure in membrane.(Pa)
        
        Return
        --------
        C: membrane mass transfer coefficient.(kg/m2/hr/Pa)
        '''
        #T = T + 273.15 # C to K
        return 7.81*d*epsilon/delta_m/tau*Mv*T**1.072/(5.685*(2*np.pi*R*Mv*T)**.5+4*d*R*(P-Pv))

.. code:: python

    # Compute Thm
    def VerifyThm(Q, Hh, Thb):
        '''
        Using Q(eq.37) compute Thm, then compare it with guess Thm. 
        Params
        --------
        Q: (kJ/(m2 hr)) 
        Hh: (kJ/(m2 K hr)) heat transfer coefficent in hotside channel.
        Returns
        ------
        Thm (verified) C 
        '''
        return Thb - Q/Hh 

.. code:: python

    # Dab
    def CalculateDab(T, P):
        '''
        Dab is diffusivity between air and water vapour.
        Equation of Dab both used in membrane & air-gap mass transfer. 
        Dab is estimated by T&p, but the dimension is m2/s(or said m2/hr).
        ## Attention! Dab influences Jv profoundly!
        Params
        -------
        T: average temperature in the membrane.(C)
        P: total pressure in the membrane.(Pa)
        
        Returns
        ---------
        Dab: (m2/hr)
        
        '''
        T = T+273.15 # C to K
        #return 1.895e-5*T**2.072/P*3600 # according to Alsaadi
        return 0.302e-2/25*T**1.75/P*3600 # according to Stephan
        #return 3e-5*3600

.. code:: python

    # Compute A
    def CalculateA(Tavg, P,air_gap_width):
        '''
        Params
        ------
        Tavg: (C) Average temperature inside the membrane.
        P: (Pa) Total pressure inside the membrane.
        Returns
        -------
        A : (kg/hr/m2/Pa) Jv/A = Pma-Pf 
        '''
        Tavg = Tavg + 273.15 # C to K
        Dab = CalculateDab(Tavg,P) #(m2/hr)
        return Dab*Mv/air_gap_width/R/Tavg

.. code:: python

    # Water vapour partial pressure
    def VaporPressure(T):
        '''
        Params
        ------
        T: (K) Temperature of water vapour(or mixture of water vapour and air).
        Attention! T in K!
        
        Returns
        -------
        P: (Pa) water vapour Partial pressure.
        '''
        b = [-7.8889166,2.5514255,-6.716169,33.239495,-105.38479,174.35319,-148.39348,48.631602]
        Tcrit = 647.25 # K 
        Pcrit = 22.093e6 # Pa
        # sum_b = np.zeros(len(T))
        sum_b = 0.
        for i in range(8):
            sum_b += b[i]*(1-T/Tcrit)**(.5*(i+1)+.5)
        Pd = Pcrit*np.exp(Tcrit/T * sum_b)
        return Pd 

.. code:: python

    # miu
    def Calculatemiu(T, S):
        '''
        Compute miu, dynamic viscosity of water.
        Parameters
        ---------
        T: (C) water temperature (20<T<180)
        S: (g/kg) salinity (0<S<130)
        Returns
        --------
        miu: (kg/m/s)
        
        '''
        uw = np.exp(-3.7942+604.129/(139.18+T))
        A = 1.474e-3 + 1.5e-6*T - 3.927e-8*T**2
        B = 1.0734e-5 - 8.5e-8*T + 2.23e-10*T**2
        uR = 1 + A*S + B*S**2
        return uw*uR*0.001 

.. code:: python

    # Compute specific heat capacity
    ## not used yet.
    def CalculateCp(T,S):
        '''
        Compute heat capacity for water.
        Parameters
        ----------
        T: (C) water temperature (10,180)
        S: (g/kg) Salinity (20,160)
        Returns 
        ----------
        Cp: (kJ/kg/C) at constant pressure.
        
        '''
        A = 4206.8 - 6.6197*S + 1.2288e-2*S**2
        B = -1.1262 + 5.4178e-2*S - 2.2719e-4*S**2 
        C = 1.2026e-2 - 5.3566e-4*S + 1.8906e-6*S**2 
        D = 6.8777e-7 + 1.517e-6*S - 4.4268e-9*S**2 
        return (A+B*T+C*T**2+D*T**3)*0.001 

.. code:: python

    # enthalpy for saturated water vapour
    def CalculateEnthalpy(T):
        '''
        Params
        ------
        T: (C) range of (0.01-200 C) Temperature of water vapour.
        
        Returns
        ----------
        hg: (kJ/kg) Enthalpy in kJ/kg.
        '''
        return 2501.689845 + 1.80692*T + 5.08772e-4*T**2

.. code:: python

    # Compute K of saline water.
    def CalculateK_water(T, S):
        '''
        This func is for liquid salt water.
        Params
        -------
        T: (C) 20-180 C Temperature. (20-80)C is acceptable.
        S: (g/kg) 0-160 Salinity.
        
        Returns
        --------
        K: (kJ/(m C hr)) Thermal conductivity.
        
        ---------------------------------------------------------------------
        For K, thermal conductivity of membrane.
        For Kab, thermal conductivity of gas mixture of air and water vapour.
        For Kf, thermal conductivity of condensed water film.
        For Kw, thermal conductivity of cooling plate. 
        
        Is K varies too small in Temperature range?
        Range of temperature: 20-80 C
        20-80C, Kab varies about 16%.
        ## Premise, Kab = 2.9e-2(50C dry air)(W/(m K)), Kw = 40(char-steel)(W/(m K)).
        # Attention! This func. has problem that results is 1000 times larger. I just make results /1000. 
        # in the range of (20-80)C, the error is acceptable. 
        '''
        A = 2e-4
        B = 3.7e-2 
        C = 3e-2 
        Log10_K = np.log10(240+A*S) + 0.434*(2.3-(343.5+B*S)/(T+273.15))*(1-(T+273.15)/(647.3+C*S))**(1./3)
        Kwater = 10**Log10_K*.001 # (W/(m K))
        return Kwater*3.6

.. code:: python

    # Compute H 
    def CalculateH(width,height,k_liquid,rho,miu,V,Pr):
        '''
        This is used for compute H in the hotside or coldside channel.
        Params
        ----------
        k_liquid: (kJ/(m K hr)) thermal conductivity of hot or cold body.
        rho: (kg/m3) use different rho when compute hot and cold water.
        miu: (kg/m/s)dynamic viscosity.
        V: (m/hr) velocity of water.
        Pr = niu/a = miu*Cp/K # varies a lot!
        
        Returns
        ---------
        H: (kJ/(m2 K hr))=(3600*kW/(m2 K)) heat transfer coefficient.
        '''
        Ks = 1.904*(SPACER_df/SPACER_hs)**(-.039)*SPACER_epsilon**.75*np.sin((np.pi/4)**.086)
        #Ks = 1.
        dh = 4.*SPACER_epsilon*SPACER_df*SPACER_hs/(2*SPACER_df + 4*(1-SPACER_epsilon)*SPACER_hs)
        #dh = 4.*width*height/(width+height)/2
        Re = rho*dh*V/3600./miu # params of the body.
        Nu = 0.029*Ks*Re**.8*Pr**.33
        '''
        print "Ks: %r" %Ks
        print "Nu: %r" %Nu
        print "Re: %r" %Re
        '''
        return Nu*k_liquid/dh

.. code:: python

    # Compute Pr
    def CalculatePr(miu,Cp,lambda_liquid):
        '''
        This func has 2 method to compute Pr, so that Nu can be calculated.
        One is linear function.
        The other is number sections.
        
        Params
        ------ 
        1.
        miu: (kg/m/s)
        cp: (kJ/kg/K)
        lambda_liquid: (kJ/(m K hr)) # attention! kJ/(...hr)
        
        2. T:(C) Temperature of water.
        
        Returns
        ---------
        Pr 
        '''
        
        return miu*Cp/(lambda_liquid/3600)
        # if T>=20 and T<40:
            # Pr = 5.5 
        # elif T>=40 and T<60:
            # Pr = 3.6 
        # elif T>=60 and T<=80:
            # Pr = 2.6
        # else:
            # raise ValueError("Input is out of range!")
        # return Pr 

.. code:: python

    # Checking nominal Numbers
    def CheckNumber(name):
        '''
        This func is for check: vapour pressure,miu,enthalpy,K_water,heat capacity,Re,Pr,diffusivity.
        '''
        # Input values.
        T_input = np.array([20,30,40,50,60,70,80]) # in C
        # Test values.
        # testvalue = np.empty_like(inputvalue,dtype=np.ndarray)
        print "While the T_inputs are:[20,30,40,50,60,70,80]C."
        print "The test values are:"
        if name == 'VapourPressure':
            Pd = VaporPressure(T_input+273.15)
            print "Pd, partial vapour pressure: %r (Pa)" %Pd # Pd, VaporPressure() is validated.
        elif name == 'miu':
            # pure water S = 0
            S_purewater = np.zeros(len(T_input))
            miu_purewater = Calculatemiu(T_input,S_purewater)
            print "miu, dynamic viscosity of pure water: %r (Pa s)" %miu_purewater # miu, Validated.
            S_salinewater = np.ones(len(T_input))*35
            miu_salinewater = Calculatemiu(T_input,S_salinewater)
            print "miu, dynamic viscosity of saline water: %r (Pa s)" %miu_salinewater # miu, Validated.
        elif name == 'enthalpy':
            enthalpy = CalculateEnthalpy(T_input)
            print "Enthalpy of water vapour: %r (kJ/kg)" %enthalpy # Enthalpy validated.
        elif name == 'K_water':
            K_water = CalculateK_water(T_input,0)
            print "K_water for pure water: %r (kJ/(m C hr))" %K_water # Validated!
        elif name == 'Pr':
            # for pure water:
            S_purewater = np.zeros(len(T_input))
            miu = Calculatemiu(T_input,S_purewater)
            #Cp = 4.19 #kJ/kg/C
            #K = .63 #W/m/C
            K_hotwater = CalculateK_water(T_input,S_purewater)
            Pr = CalculatePr(miu,Cp_hot,K_hotwater)
            print "Pr for pure hot water is: %r" %Pr # Pr is Validated.
        elif name == 'diffusivity':
            P = 101325 # Pa
            T = np.array([0.,25.,40.,60.]) # C
            D = CalculateDab(T,P) # m2/hr
            D = D/3600. # m2/s
            print "T at %r C" %T
            print "Diffusivity between water vapour and air: %r*10^(-5)" %(D*10**5) # NOT! but in the range(20-80)is ok.
        elif name == 'H':
            V = 1.5*3600 #(m/hr) average velocity of water
            miu = Calculatemiu(40.,0)
            Pr = 4.102 #CalculatePr(miu_hot,Cp_hot,K_hotwater)
            rho = 991. #(kg/m3)
            k_l = 0.638*3.6 #(W/(m K))
            width = 1.
            height = 1.
            H = CalculateH(width,height,k_l,rho,miu,V,Pr)
            print "H, for Heat Transfer a example: %r (kJ/(m2 K hr))= %r (W/(m2 K))" %(H,H/3.6) # Validated, but whether dh is right?
            print "miu: %r(Pa s)" %miu
        else:
            raise NameError("My Error! No such function. Please input a valid string!")

.. code:: python

    CheckNumber('H')


.. parsed-literal::

    While the T_inputs are:[20,30,40,50,60,70,80]C.
    The test values are:
    H, for Heat Transfer a example: 64561.367537425169 (kJ/(m2 K hr))= 17933.713204840325 (W/(m2 K))
    miu: 0.0006553616086282997(Pa s)
    

.. code:: python

    # Computation for Jv
    def ComputationForJv(d,delta_m,tau,Tave,Thm,Tf,air_gap_width):
        '''
        For neat code in the Loop.
        '''
        Pv = VaporPressure(Tave+273.15) # water vapour pressure inside membrane
        C = CalculateC(d,delta_m,tau,Tave,Ptotal,Pv)
        A = CalculateA(.5*(Thm+Tf),Ptotal,air_gap_width)
        Phm = VaporPressure(Thm+273.15)
        Pf = VaporPressure(Tf+273.15)
        return C,A,Phm,Pf

.. code:: python

    # Computation for Q
    def ComputationForQ(width,height,mh,Thm,Tf,Thb,S_hot,K_hotwater):
        '''
        Forc neat code in the Loop.
        '''
        Vh = mh/rho_hotwater/(width*height)#(m/hr) average velocity of water
        hg = CalculateEnthalpy(.5*(Thm+Tf))
        miu_hot = Calculatemiu(Thb,S_hot)
        Pr_hot = CalculatePr(miu_hot,Cp_hot,K_hotwater)
        Hh = CalculateH(width,height,K_hotwater,rho_hotwater,miu_hot,Vh,Pr_hot)
        return hg,Hh

.. code:: python

    # Computation for Tf
    def ComputationForTf(width,height,mc,Sc,Tcb,Tf):
        '''For neat code.'''
        Vc = mc/rho_coldwater/(width*height)#(m/hr) average velocity of water
        miu_cold = Calculatemiu(Tcb,Sc)
        K_coldwater = CalculateK_water(Tcb,Sc)
        K_fallingwater = CalculateK_water(Tf,0)
        Pr_cold = CalculatePr(miu_cold,Cp_cold,K_coldwater)
        Hc = CalculateH(width,height,K_coldwater,rho_coldwater,miu_cold,Vc,Pr_cold)
        return Hc,K_fallingwater

.. code:: python

    # Checking
    def Checking(Jv,delta_f,Q,Tf,temp):
        '''
        This func is only for checking inter-loop results by print them.
        
        '''
        print 
        print "This is checking function:"
        print 'Jv, water vapor flux: %r(kg/m2/hr)' %Jv
        print "delta_f: %r (m)" %delta_f 
        print "Q, transport energy: %r(kJ/m2/hr)" %Q 
        print "Tf,temperature of falling water flim: %r(C)" %Tf 
        print "The criterion value is: %r " %temp

.. code:: python

    class AGMD:
        '''Using for computing AGMD process. Single input and equal output.'''
        def __init__(self, _L,_W,_H,_delta_a,_delta_c,_Kw,_d,_delta_m,_tau,_mhi,_mci,_thi,_tci,_shi,_sc):
            '''
            Setting up the module features and operation condition.
            Params
            -------
            _L,_W: (m) Length scale of the module
            _H: (m) height of the channel
            _delta_a: (m) distance of air gap
            _delta_c: (m) thickness of membrane
            _Kw: (kJ/m2 K hr) thermal conductivity of cooling plate
            _d: (m) diameter of pore size inside membrane
            _delta_m: (m) thickness of membrane
            _tau: torsion of membrane
            _mhi: (kg/hr) mass of hot water inlet
            _mci: (kg/hr) mass of cooling water inlet(never change)
            _thi: (C) temperature of hot water inlet
            _tci: (C) temperature of cooling water inlet
            _shi: (g/kg) salinity of hot water inlet
            _sc: (g/kg) salinity of cooling water(never change)
            '''
            # Module characters.
            self.length_vertical_effective = _L
            self.width_effective = _W
            self.height_channel = _H
            self.delta_a = _delta_a
            #self.thermal_conductivity_air = _thermal_conductivity_air # It is a Global params.
            self.delta_c = _delta_c
            self.Kw = _Kw
            
            ## For membrane
            self.diameter_pore_size = _d
            self.delta_m = _delta_m
            self.tau_membrane = _tau
            #self.porosity_membrane = _porosity_membrane # Generally it is about 78-80%, assuming a Global parameter.
            
            # Operation conditions.
            self.mass_hot_inlet = _mhi
            self.mass_cold_inlet = _mci
            self.temperature_hot_inlet = _thi
            self.temperature_cold_inlet = _tci
            self.salinity_hot_inlet = _shi
            self.salinity_cold = _sc
            print "Next step: Setup mesh please."
            
        def get_mesh_setup(self, Nx):
            '''
            Mesh setup and initializes Numpy array.
            Params
            ------
            Nx, numerber of grids
            
            Generates
            -----------
            x: 1-D array of Nx floats, stores dimension info.
            Thb,Tcb,Tf,Thm,Tma: 1-D array of Nx floats, stores temperature info.
            mh,S_hot: 1-D array of Nx floats, stores feed water info.
            '''
            
            # Nx # Mesh grid number. value from input
            self.x = np.linspace(0,self.length_vertical_effective,Nx)
            self.dx = self.length_vertical_effective/(Nx-1)
            
            self.Thb = np.ones(Nx)*self.temperature_hot_inlet # the body temperature of hot water inlet.
            self.Tcb = np.ones(Nx)*self.temperature_cold_inlet # the body temperature of cooling water inlet.
            self.Tf = np.ones(Nx) # The T of condensing film interface.
            self.Thm = np.ones(Nx) # The T of the interface of hotfeed and membrane.
            self.Tma = np.ones(Nx) # The T of the surface of membrane facing the air channel.
            #Tfw = np.ones(Nx) # The T of film water in contact with cooling plate.
            #Tcw = np.ones(Nx) # The T of the wall in cooling channel.
            self.mh = np.ones(Nx)*self.mass_hot_inlet # mass flow of hot water at every grid.
            self.S_hot = np.ones(Nx)*self.salinity_hot_inlet # salinity of hot water.
            ## store output values wanted to show.
            self.Jv_flux_condensed_water_along = np.zeros(Nx) # condensed water flux along the membrane
            self.Q_heattransfer_along = np.zeros(Nx) # heat transfer along the membrane
            self.delta_f_condensed_water_along = np.zeros(Nx) #delta_f along the membrane
            print "Mesh setup done."
            
            
        def get_co_current(self, Nx):
            '''To calculate co-current regime.
            Params
            ------
            Nx: grids of the mesh.
            
            Generates
            ---------
            Jv_flux_condensed_water_along: kg/m2/hr. 1-D array, with Nx floats.
            Q_heattransfer_along: kJ/m2/hr. 1-D array, with Nx floats.
            delta_f_condensed_water_along: m. Water film along the cooling plate. 1-D array, with Nx floats.
            '''
            print "This is co-current flow pattern."
            print "Hot water inlet at %d C" %self.temperature_hot_inlet
            print "Cooling water inlet at %d C" %self.temperature_cold_inlet
            print "Hot feed at %.2f LPM" %(self.mass_hot_inlet/60.)
            print "Cooling mass inlet: %.2f LPM"%(self.mass_cold_inlet/60)
            print "Air gap distance: %.2f mm"%(self.delta_a*1000)
            print 
    
            for i in range(Nx):
                #print "Step %r of %r." %(i,Nx)
                error_Tf = 10.
                self.Tf[i] = .5*(self.Thb[i]+self.Tcb[i]) # supposed
                Tf_step = 0
                while error_Tf > 0.1 or Tf_step <= 20:
                    Tf_step += 1
                    error_Thm = 10.
                    self.Thm[i] = .5*(self.Thb[i]+self.Tf[i]) # supposed
                    error_Thm_threshold = 1e-1 # threshold for error(while)
                    Thm_step = 0
                    while error_Thm > error_Thm_threshold or Thm_step <= 30 :
                        Thm_step += 1
                        Tave = .5*(self.Thm[i]+self.Tma[i]) # average temperature inside membrane
                        C,A,Phm,Pf = ComputationForJv(self.diameter_pore_size,self.delta_m,self.tau_membrane,\
                                                      Tave,self.Thm[i],self.Tf[i],self.delta_a)
                        Jv = (1./C + 1./A)**(-1)*(Phm-Pf) #Compute C,A,Pf,Pmh
    
                        K_feedwater = CalculateK_water(self.Thb[i],self.S_hot[i])
                        hg,Hh = ComputationForQ(self.width_effective,self.height_channel,self.mh[i],self.Thm[i],\
                                                self.Tf[i],self.Thb[i],self.S_hot[i],K_feedwater)
                        Q = (self.Thb[i] - self.Tf[i] + Jv*hg*(self.delta_m/K+self.delta_a/Kab))/(1./Hh+self.delta_m/K+self.delta_a/Kab) 
                        #Compute hg,Hh. K,delta_a/m are Const. T depends on grid.
    
                        Thm_temp = VerifyThm(Q,Hh,self.Thb[i])
                        error_Thm = abs(Thm_temp - self.Thm[i])
                        self.Thm[i] = Thm_temp
                        #Checking
                        #Checking(Jv,"no deltaf",Q,Tf[i],Thm[i])
    
                    #print 
                    #print "## 1st Loop ends here. Use kepboard to interrupt."
                    #print 'Difference of Thm between steps is: %5.5f' %error_Thm
                    #print 
                    #raw_input()
                    miu_f = Calculatemiu(self.Tf[i],0)
                    delta_f = (3.*Jv/3600*miu_f/(rho_l*(rho_l - rho_av)*G)*self.x[i])**(1./3) 
                    #rho_av,rho_l,g is Const. Compute miu. x depends on grid.
    
                    Hc,Kf = ComputationForTf(self.width_effective,self.height_channel,\
                                             self.mass_cold_inlet,self.salinity_cold,self.Tcb[i],self.Tf[i])
                    Tf_temp = Q*(delta_f/Kf+self.delta_c/self.Kw+1./Hc) + self.Tcb[i] 
                    # Kw,delta_c/f are Const. Hc,Kf needs computed. T depends on grid.
    
                    error_Tf = abs(Tf_temp - self.Tf[i])
                    self.Tf[i] = Tf_temp
                    self.Tma[i] = (Kab/self.delta_a*self.Tf[i] + K/self.delta_m*self.Thm[i])/(Kab/self.delta_a + K/self.delta_m)
                    # Checking 
                    #Checking(Jv,delta_f,Q,Tf[i],"Tf[i]")
    
                #print 
                #print "## 2nd Loop ends here. Use keyboard to interrupt."
                #print 'Difference of Thm between steps is: %r' %error_Thm
                #print 'Difference of Tf between steps is: %r' %error_Tf
                #print 
                # raw_input()
                # Next state
                if i < (Nx-1):
                    self.mh[i+1] = self.mh[i]-Jv*self.dx*self.width_effective
                    #print "mh[%r]: %r" %(i,mh[i])
                    self.S_hot[i+1] = self.mh[i]*self.S_hot[i]/self.mh[i+1]
                    #print "S_hot[%r]: %r" %(i,S_hot[i])
                    # compute Cph,Cpc separately
                    self.Thb[i+1] = (self.mh[i]*Cp_hot*self.Thb[i] - Q*self.dx*self.width_effective)/Cp_hot/self.mh[i+1]
                    #print "Thb[%r]: %r" %(i,Thb[i])
                    self.Tcb[i+1] = self.Tcb[i] + Q*self.dx*self.width_effective/(Cp_cold*self.mass_cold_inlet)
                    #print "Tcb[%r]: %r" %(i,Tcb[i])
                else:
                    print "The end. i = %d."%i
    
                # Outputs Jv,Q,delta_f arrays
                self.Jv_flux_condensed_water_along[i] = Jv # condensed water flux along the membrane
                self.Q_heattransfer_along[i] = Q # heat transfer along the membrane
                self.delta_f_condensed_water_along[i] = delta_f #delta_f along the membrane
    
            #Checking output
            print "Results are: "
            print "Jv along the length: %r (kg/(m2 hr))" %self.Jv_flux_condensed_water_along
            #print "Q along the length: %r (kJ/m2/hr)" %Q_N
            #print "delta_f along the length: %r (m)" %delta_fN
            print "Average permeate water: %-5.5f kg/hr/m2" %((sum(self.Jv_flux_condensed_water_along)-\
                                                           self.Jv_flux_condensed_water_along[-1])*self.dx/self.length_vertical_effective) 
            # shape(Jv_N)=number(dx)+1
            print "Average of sum(flux)/Nx: %.2f " %(sum(self.Jv_flux_condensed_water_along)/Nx)
            
        def get_counter_current(self,Nx):
            '''
            This func is to calculate counter-current regime.
            Params
            -------
            Nx: grids of the mesh.
            
            Generates
            ---------
            Jv_flux_condensed_water_along: kg/m2/hr. 1-D array, with Nx floats.
            Q_heattransfer_along: kJ/m2/hr. 1-D array, with Nx floats.
            delta_f_condensed_water_along: m. Water film along the cooling plate. 1-D array, with Nx floats.
            '''
            print "This is counter-current flow pattern."
            print "Hot water inlet at %d C" %self.temperature_hot_inlet
            print "Cooling water inlet at %d C" %self.temperature_cold_inlet
            print "Hot feed at %.2f LPM" %(self.mass_hot_inlet/60.)
            print "Cooling water inlet at %d C" %self.temperature_cold_inlet
            print "Hot feed at %.2f LPM" %(self.mass_hot_inlet/60.)
            print "Cooling mass inlet: %.2f LPM"%(self.mass_cold_inlet/60)
            print "Air gap distance: %.2f mm"%(self.delta_a*1000)
            print 
            
            learnRate = 1.
            error_Tcb = 10.
            error_i = 0. # alpha/(error_i+1) to make alpha decrease every step.
            self.Tcb[0] = .5*(self.Thb[0]+self.Tcb[-1]) # supposed
            while error_Tcb > 0.01 or error_i <= 20:
    
                for i in range(Nx):
                    #print "Step %r of %r." %(i,Nx)
                    error_Tf = 10.
                    self.Tf[i] = .5*(self.Thb[i]+self.Tcb[i]) # supposed
                    Tf_step = 0
                    while error_Tf > 0.1 or Tf_step <= 20:
                        Tf_step += 1
                        error_Thm = 10.
                        self.Thm[i] = .5*(self.Thb[i]+self.Tf[i]) # supposed
                        error_Thm_threshold = 1e-1 # threshold for error(while)
                        Thm_step = 0
                        while error_Thm > error_Thm_threshold or Thm_step <= 30:
                            Thm_step += 1
                            Tave = .5*(self.Thm[i]+self.Tma[i]) # average temperature inside membrane
                            C,A,Phm,Pf = ComputationForJv(self.diameter_pore_size,self.delta_m,self.tau_membrane,\
                                                          Tave,self.Thm[i],self.Tf[i],self.delta_a)
                            Jv = (1./C + 1./A)**(-1)*(Phm-Pf) #Compute C,A,Pf,Pmh
    
                            K_feedwater = CalculateK_water(self.Thb[i],self.S_hot[i])
                            hg,Hh = ComputationForQ(self.width_effective,self.height_channel,self.mh[i],self.Thm[i],\
                                                    self.Tf[i],self.Thb[i],self.S_hot[i],K_feedwater)
                            Q = (self.Thb[i] - self.Tf[i] + Jv*hg*(self.delta_m/K+self.delta_a/Kab))/(1./Hh+self.delta_m/K+self.delta_a/Kab) 
                            #Compute hg,Hh. K,delta_a/m are Const. T depends on grid.
    
                            Thm_temp = VerifyThm(Q,Hh,self.Thb[i])
                            error_Thm = abs(Thm_temp - self.Thm[i])
                            self.Thm[i] = Thm_temp
                            #Checking
                            #Checking(Jv,"no deltaf",Q,Tf[i],Thm[i])
    
                        #print 
                        #print "## 1st Loop ends here. Use kepboard to interrupt."
                        #print 'Difference of Thm between steps is: %r' %error_Thm
                        #print 
                        #raw_input()
                        miu_f = Calculatemiu(self.Tf[i],0)
                        delta_f = (3.*Jv/3600*miu_f/(rho_l*(rho_l - rho_av)*G)*self.x[i])**(1./3) 
                        #rho_av,rho_l,g is Const. Compute miu. x depends on grid.
    
                        Hc,Kf = ComputationForTf(self.width_effective,self.height_channel,
                                                 self.mass_cold_inlet,self.salinity_cold,self.Tcb[i],self.Tf[i])
                        Tf_temp = Q*(delta_f/Kf+self.delta_c/self.Kw+1./Hc) + self.Tcb[i] 
                        # Kw,delta_c/f are Const. Hc,Kf needs computed. T depends on grid.
    
                        error_Tf = abs(Tf_temp - self.Tf[i])
                        self.Tf[i] = Tf_temp
                        self.Tma[i] = (Kab/self.delta_a*self.Tf[i] + K/self.delta_m*self.Thm[i])/(Kab/self.delta_a + K/self.delta_m)
                        # Checking 
                        #Checking(Jv,delta_f,Q,Tf[i],"Tf[i]")
    
                    #print 
                    #print "## 2nd Loop ends here. Use keyboard to interrupt."
                    #print 'Difference of Thm between steps is: %r' %error_Thm
                    #print 'Difference of Tf between steps is: %r' %error_Tf
                    #print 
                    # raw_input()
                    # Next state
                    if i < (Nx-1):
                        self.mh[i+1] = self.mh[i]-Jv*self.dx*self.width_effective
                        #print "mh[%r]: %r" %(i,mh[i])
                        self.S_hot[i+1] = self.mh[i]*self.S_hot[i]/self.mh[i+1]
                        #print "S_hot[%r]: %r" %(i,S_hot[i])
                        # compute Cph,Cpc separately
                        self.Thb[i+1] = (self.mh[i]*Cp_hot*self.Thb[i] - Q*self.dx*self.width_effective)/Cp_hot/self.mh[i+1]
                        #print "Thb[%r]: %r" %(i,Thb[i])
                        self.Tcb[i+1] = self.Tcb[i] + Q*self.dx*self.width_effective/(Cp_cold*self.mass_cold_inlet)
                        #print "Tcb[%r]: %r" %(i,Tcb[i])
                    else:
                        print "The end. i = %d."%i
    
                    # Outputs Jv,Q,delta_f arrays
                    self.Jv_flux_condensed_water_along[i] = Jv # condensed water flux along the membrane
                    self.Q_heattransfer_along[i] = Q # heat transfer along the membrane
                    self.delta_f_condensed_water_along[i] = delta_f #delta_f along the membrane
    
                
                temp_move = self.Tcb[-1] - self.temperature_cold_inlet
                print "Tcb[-1]: %.3f C" %self.Tcb[-1]
                error_Tcb = abs(temp_move)
    
                alpha = 1./(error_i+1)*learnRate
                error_i += 1 
                self.Tcb[0] = self.Tcb[0] - temp_move*alpha
            #checking output
            print 'Difference of Tcb between steps is: %r' %error_Tcb
            print 
            print "Results are: "
            print "Jv along the length: %r (kg/(m2 hr))" %self.Jv_flux_condensed_water_along
            #print "Q along the length: %r (kJ/m2/hr)" %Q_N
            #print "delta_f along the length: %r (m)" %delta_fN
            print "Average permeate water: %r kg/hr/m2" %((sum(self.Jv_flux_condensed_water_along)-\
                                                           self.Jv_flux_condensed_water_along[-1])*self.dx/self.length_vertical_effective) 
            # shape(Jv_N)=number(dx)+1
            print "Average of sum(flux)/Nx: %.2f " %(sum(self.Jv_flux_condensed_water_along)/Nx)

.. code:: python

    def plot_along_membrane(plot_title,plot_ylabel,x_values,y_values,y_limits):
        '''
        Pyplot to show T, mh, Jv, delta_f along the membrane.
        Params
        ------
        plot_title,plot_ylabel: strings, to describe Title and Y-Label of the fig plotted.
        x_values,y_values: 1-D array of floats.
        y_limits: 1-D array with 2 elements which defines the lower/upper boundary of the plot y-value
                    i.e. y_limits = [lowerbound, upperbound]
                    if no concerns for boundary, input y_limits = "no boundary" instead.
        '''
        pyplot.figure(figsize=(10,5))
        pyplot.title(plot_title);
        pyplot.xlabel('x(m)');
        pyplot.ylabel(plot_ylabel);
        if type(y_limits) == str:
            print "No specific restrict boundary."
        elif type(y_limits) == list:
            y_max = y_limits[1]
            y_min = y_limits[0]
            pyplot.ylim(y_min,y_max)
        else : 
            raise NameError("My Error. y_limits input is wrong.")
        pyplot.scatter(x_values,y_values,marker='^')
        pyplot.plot(x_values,y_values,color='#003366',ls='-')

.. code:: python

    def plot_temperature(plot_title,plot_ylabel1,plot_ylabel2,plot_ylabel3,plot_ylabel4,plot_ylabel5,\
                         x_values,y_values1,y_values2,y_values3,y_values4,y_values5,y_limits):
        '''
        Pyplot to show T, mh, Jv, delta_f along the membrane.
        Params
        ------
        plot_title,plot_ylabel_N: strings, to describe Title and Y-Label of the fig plotted.
        x_values: x along the membrane. 1-D array of floats.
        y_values_N: (C) Temperature along the membrane. 1-D array of floats.
        y_limits: 1-D array with 2 elements which defines the lower/upper boundary of the plot y-value
                    i.e. y_limits = [lowerbound, upperbound]
                    if no concerns for boundary, input y_limits = "no boundary" instead.
        '''
        pyplot.figure(figsize=(10,5))
        pyplot.title(plot_title);
        pyplot.xlabel('x(m)');
        pyplot.ylabel("Temperature (C)");
        if type(y_limits) == str:
            print "No specific restrict boundary."
        elif type(y_limits) == list:
            y_max = y_limits[1]
            y_min = y_limits[0]
            pyplot.ylim(y_min,y_max)
        else : 
            raise NameError("My Error. y_limits input is wrong.")
        #pyplot.scatter(x_values,y_values,marker='^')
        pyplot.plot(x_values,y_values1,color='#003366',ls='-')
        pyplot.plot(x_values,y_values2,ls='--')
        pyplot.plot(x_values,y_values3)
        pyplot.plot(x_values,y_values4,ls='-.')
        pyplot.plot(x_values,y_values5)
        pyplot.legend((plot_ylabel1,plot_ylabel2,plot_ylabel3,plot_ylabel4,plot_ylabel5))

.. code:: python

    help(AGMD)


.. parsed-literal::

    Help on class AGMD in module __main__:
    
    class AGMD
     |  Using for computing AGMD process. Single input and equal output.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, _L, _W, _H, _delta_a, _delta_c, _Kw, _d, _delta_m, _tau, _mhi, _mci, _thi, _tci, _shi, _sc)
     |      Setting up the module features and operation condition.
     |      Params
     |      -------
     |      _L,_W: (m) Length scale of the module
     |      _H: (m) height of the channel
     |      _delta_a: (m) distance of air gap
     |      _delta_c: (m) thickness of membrane
     |      _Kw: (kJ/m2 K hr) thermal conductivity of cooling plate
     |      _d: (m) diameter of pore size inside membrane
     |      _delta_m: (m) thickness of membrane
     |      _tau: torsion of membrane
     |      _mhi: (kg/hr) mass of hot water inlet
     |      _mci: (kg/hr) mass of cooling water inlet(never change)
     |      _thi: (C) temperature of hot water inlet
     |      _tci: (C) temperature of cooling water inlet
     |      _shi: (g/kg) salinity of hot water inlet
     |      _sc: (g/kg) salinity of cooling water(never change)
     |  
     |  get_co_current(self, Nx)
     |      To calculate co-current regime.
     |      Params
     |      ------
     |      Nx: grids of the mesh.
     |      
     |      Generates
     |      ---------
     |      Jv_flux_condensed_water_along: kg/m2/hr. 1-D array, with Nx floats.
     |      Q_heattransfer_along: kJ/m2/hr. 1-D array, with Nx floats.
     |      delta_f_condensed_water_along: m. Water film along the cooling plate. 1-D array, with Nx floats.
     |  
     |  get_counter_current(self, Nx)
     |      This func is to calculate counter-current regime.
     |      Params
     |      -------
     |      Nx: grids of the mesh.
     |      
     |      Generates
     |      ---------
     |      Jv_flux_condensed_water_along: kg/m2/hr. 1-D array, with Nx floats.
     |      Q_heattransfer_along: kJ/m2/hr. 1-D array, with Nx floats.
     |      delta_f_condensed_water_along: m. Water film along the cooling plate. 1-D array, with Nx floats.
     |  
     |  get_mesh_setup(self, Nx)
     |      Mesh setup and initializes Numpy array.
     |      Params
     |      ------
     |      Nx, numerber of grids
     |      
     |      Generates
     |      -----------
     |      x: 1-D array of Nx floats, stores dimension info.
     |      Thb,Tcb,Tf,Thm,Tma: 1-D array of Nx floats, stores temperature info.
     |      mh,S_hot: 1-D array of Nx floats, stores feed water info.
    
    

.. code:: python

    test_alsaadi = AGMD(_L=.1,_W=.05,_H=2e-3,_delta_a=9e-3,_delta_c=.25e-3,_Kw=40.*3.6,\
                        _d=.2e-6,_delta_m=100e-6,_tau=1.5,_mhi=1.5*60,_mci=1.5*60,_thi=80.,_tci=20,_shi=42.,_sc=0.)
    Nx = 20
    test_alsaadi.get_mesh_setup(Nx)
    test_alsaadi.get_co_current(Nx)


.. parsed-literal::

    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 80 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 6.5655002 ,  6.55948067,  6.55535168,  6.55142388,  6.54759408,
            6.54382493,  6.54009817,  6.5364033 ,  6.53273368,  6.52908482,
            6.52545349,  6.52183731,  6.51823446,  6.51464352,  6.51106334,
            6.50749298,  6.50393168,  6.5003788 ,  6.49683378,  6.49329616]) (kg/(m2 hr))
    Average permeate water: 6.52955 kg/hr/m2
    Average of sum(flux)/Nx: 6.53 
    

.. code:: python

    Jv_exp_fig8 = np.array([0.45,1.04,1.84,3.40,5.93])
    print Jv_exp_fig8[np.argsort(-Jv_exp_fig8)]
    np.argsort(-Jv_exp_fig8)


.. parsed-literal::

    [ 5.93  3.4   1.84  1.04  0.45]
    



.. parsed-literal::

    array([4, 3, 2, 1, 0])



.. code:: python

    # For testing Alsaadi fig.8
    Nx = 20
    Temperature_hot_inlet = np.array([80,70,60,50,40]) # C
    Jv_exp_fig8 = np.array([0.45,1.04,1.84,3.40,5.93])
    Jv_exp_fig8_re = Jv_exp_fig8[np.argsort(-Jv_exp_fig8)]
    Jv_simulate_N = np.zeros(len(Temperature_hot_inlet))
    L1_norm = []
    for T_i,thi_N in enumerate(Temperature_hot_inlet):
        verify_fig8 = AGMD(_L=.1,_W=.05,_H=2e-3,_delta_a=9e-3,_delta_c=.25e-3,_Kw=40.*3.6,\
                        _d=.2e-6,_delta_m=100e-6,_tau=1.5,_mhi=1.5*60,_mci=1.5*60,_thi=thi_N,_tci=20,_shi=42.,_sc=0.)
        verify_fig8.get_mesh_setup(Nx)
        verify_fig8.get_co_current(Nx);
        Jv_simulate_N[T_i] = sum(verify_fig8.Jv_flux_condensed_water_along)/Nx
        L1_norm.append((Jv_simulate_N[T_i] - Jv_exp_fig8_re[T_i])/Jv_exp_fig8_re[T_i])
    print "Single error: L1-norm: %r " %L1_norm
    L2_norm = sum(((Jv_simulate_N - Jv_exp_fig8_re)/Jv_exp_fig8_re)**2)**(.5)
    print "Total L2-norm: %r " %L2_norm


.. parsed-literal::

    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 80 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 5.87354825,  5.86872519,  5.86541922,  5.86227444,  5.85920803,
            5.85619004,  5.85320581,  5.85024692,  5.84730805,  5.84438557,
            5.84147693,  5.83858019,  5.83569391,  5.83281694,  5.82994836,
            5.82708742,  5.82423351,  5.82138612,  5.8185448 ,  5.81570918]) (kg/(m2 hr))
    Average permeate water: 5.84475 kg/hr/m2
    Average of sum(flux)/Nx: 5.84 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 70 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 3.6725779 ,  3.6702479 ,  3.66869893,  3.66723255,  3.66580612,
            3.66440424,  3.66301934,  3.66164711,  3.6602848 ,  3.65893055,
            3.65758304,  3.65624129,  3.65490454,  3.65357221,  3.65224383,
            3.65091901,  3.64959744,  3.64827885,  3.64696302,  3.64564976]) (kg/(m2 hr))
    Average permeate water: 3.65911 kg/hr/m2
    Average of sum(flux)/Nx: 3.66 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 2.16937536,  2.16832477,  2.16764882,  2.16701236,  2.16639493,
            2.16578915,  2.16519139,  2.16459958,  2.16401241,  2.16342898,
            2.16284866,  2.16227098,  2.16169557,  2.16112216,  2.16055053,
            2.15998048,  2.15941187,  2.15884457,  2.15827847,  2.15771348]) (kg/(m2 hr))
    Average permeate water: 2.16351 kg/hr/m2
    Average of sum(flux)/Nx: 2.16 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 50 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 1.18120875,  1.18077875,  1.18051045,  1.18025915,  1.18001604,
            1.17977792,  1.17954322,  1.17931106,  1.17908086,  1.17885225,
            1.17862494,  1.17839874,  1.17817349,  1.17794906,  1.17772537,
            1.17750233,  1.17727987,  1.17705795,  1.17683651,  1.17661552]) (kg/(m2 hr))
    Average permeate water: 1.17889 kg/hr/m2
    Average of sum(flux)/Nx: 1.18 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 40 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 0.56084889,  0.56069924,  0.56060776,  0.5605224 ,  0.56043996,
            0.56035932,  0.5602799 ,  0.56020139,  0.56012358,  0.56004632,
            0.55996954,  0.55989314,  0.55981707,  0.5597413 ,  0.55966578,
            0.55959049,  0.5595154 ,  0.5594405 ,  0.55936576,  0.55929118]) (kg/(m2 hr))
    Average permeate water: 0.56006 kg/hr/m2
    Average of sum(flux)/Nx: 0.56 
    Single error: L1-norm: [-0.014620667154675783, 0.076011800972699023, 0.17566561153555416, 0.1334376078259899, 0.24449098968385319] 
    Total L2-norm: 0.33827711138981004 
    

.. code:: python

    def checkingT(check_instance):
        print "Thb along the membrane %r C" %check_instance.Thb
        print "Thm along the membrane %r C" %check_instance.Thm
        print "Tma along the membrane %r C" %check_instance.Tma
        print "Tf along the membrane %r C" %check_instance.Tf
        print "Tcb along the membrane %r C" %check_instance.Tcb
    checkingT(test_alsaadi)


.. parsed-literal::

    Thb along the membrane array([ 80.        ,  79.98917119,  79.97835207,  79.96753934,
            79.95673265,  79.94593183,  79.93513676,  79.92434739,
            79.91356364,  79.90278549,  79.8920129 ,  79.88124583,
            79.87048426,  79.85972817,  79.84897754,  79.83823235,
            79.82749258,  79.81675823,  79.80602927,  79.79530569]) C
    Thm along the membrane array([ 79.15568776,  79.14557057,  79.13520734,  79.12482342,
            79.11443241,  79.10403919,  79.09364614,  79.0832546 ,
            79.07286541,  79.06247913,  79.05209617,  79.04171682,
            79.0313413 ,  79.02096977,  79.01060238,  79.00023923,
            78.98988041,  78.97952599,  78.96917603,  78.95883058]) C
    Tma along the membrane array([ 78.92828901,  78.91870837,  78.90854756,  78.89833094,
            78.88809036,  78.8778373 ,  78.86757736,  78.85731373,
            78.84704843,  78.83678281,  78.82651786,  78.81625427,
            78.8059926 ,  78.79573326,  78.78547658,  78.77522285,
            78.76497228,  78.75472506,  78.74448135,  78.73424127]) C
    Tf along the membrane array([ 21.45443129,  21.58046061,  21.62146042,  21.65352858,
            21.68130857,  21.70647871,  21.72986087,  21.75192642,
            21.7729738 ,  21.79320542,  21.81276558,  21.83176126,
            21.85027423,  21.86836854,  21.88609547,  21.90349672,
            21.92060678,  21.93745449,  21.9540643 ,  21.97045708]) C
    Tcb along the membrane array([ 20.        ,  20.01238214,  20.02475275,  20.03711559,
            20.04947104,  20.06181931,  20.07416051,  20.08649473,
            20.09882203,  20.11114246,  20.12345606,  20.13576287,
            20.14806291,  20.16035622,  20.17264281,  20.1849227 ,
            20.19719592,  20.20946249,  20.22172241,  20.2339757 ]) C
    

.. code:: python

    # For testing Alsaadi air gap fig9.
    Nx = 20
    delta_a_Nvalues = np.array([0.005,0.009,0.013]) # m
    TEST = [] # List to store class-instance
    Jv_delta_aN = np.zeros(len(delta_a_Nvalues))
    Jv_exp_fig9 = np.array([3.17,1.83,1.53])
    L1_norm = []
    for test_i,airgap_distance in enumerate(delta_a_Nvalues):
        ex_airgap = AGMD(_L=.1,_W=.05,_H=2e-3,_delta_a=airgap_distance,_delta_c=.25e-3,_Kw=40.*3.6,\
                        _d=.2e-6,_delta_m=100e-6,_tau=1.5,_mhi=1.5*60,_mci=1.5*60,_thi=60.,_tci=20.,_shi=42.,_sc=0.)
        ex_airgap.get_mesh_setup(Nx)
        ex_airgap.get_co_current(Nx);
        Jv_delta_aN[test_i] = sum(ex_airgap.Jv_flux_condensed_water_along)/Nx #Jv (average along membrane) for each instance
        #TEST.append(ex_airgap)
        #print test_i,TEST[test_i].Thm
        L1_norm.append((Jv_delta_aN[test_i] - Jv_exp_fig9[test_i])/Jv_exp_fig9[test_i])
    print "Single error: L1-norm: %r " %L1_norm
    #print Jv_varies_delta_a


.. parsed-literal::

    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 5.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 2.84175581,  2.83980004,  2.83856745,  2.83741123,  2.83629189,
            2.83519516,  2.83411404,  2.83304452,  2.83198408,  2.83093098,
            2.829884  ,  2.82884223,  2.82780497,  2.82677169,  2.82574194,
            2.82471537,  2.82369169,  2.82267066,  2.82165205,  2.8206357 ]) (kg/(m2 hr))
    Average permeate water: 2.83110 kg/hr/m2
    Average of sum(flux)/Nx: 2.83 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 2.16937536,  2.16832477,  2.16764882,  2.16701236,  2.16639493,
            2.16578915,  2.16519139,  2.16459958,  2.16401241,  2.16342898,
            2.16284866,  2.16227098,  2.16169557,  2.16112216,  2.16055053,
            2.15998048,  2.15941187,  2.15884457,  2.15827847,  2.15771348]) (kg/(m2 hr))
    Average permeate water: 2.16351 kg/hr/m2
    Average of sum(flux)/Nx: 2.16 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 13.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 1.75327297,  1.75262105,  1.75219433,  1.75179133,  1.75139975,
            1.75101515,  1.75063534,  1.7502591 ,  1.74988561,  1.74951436,
            1.74914495,  1.74877711,  1.74841061,  1.7480453 ,  1.74768103,
            1.74731769,  1.7469552 ,  1.74659348,  1.74623245,  1.74587208]) (kg/(m2 hr))
    Average permeate water: 1.74957 kg/hr/m2
    Average of sum(flux)/Nx: 1.75 
    Single error: L1-norm: [-0.10707404573646438, 0.18209001378438233, 0.14338623822601926] 
    

.. code:: python

    # For plot Jv versus delta_a. See fig10, Banat 1998
    # Reuse plot(not for membrane along)
    plot_along_membrane(plot_title="Exam Jv vs. delta_a",plot_ylabel="J v_ave(kg/m2 hr)",x_values=delta_a_Nvalues,\
                        y_values=Jv_delta_aN,y_limits=[1.,3.])



.. image:: output_30_0.png


.. code:: python

    # For testing Alsaadi fig.10
    Nx = 20
    Temperature_hot_inlet = np.array([80,70,60,50,40]) # C
    Jv_exp_fig10 = np.array([0.55,1.10,2.07,3.17,6.59])
    Jv_exp_fig10_re = Jv_exp_fig8[np.argsort(-Jv_exp_fig8)]
    Jv_simulate_N = np.zeros(len(Temperature_hot_inlet))
    L1_norm = []
    for T_i,thi_N in enumerate(Temperature_hot_inlet):
        verify_fig10 = AGMD(_L=.1,_W=.05,_H=2e-3,_delta_a=9e-3,_delta_c=.25e-3,_Kw=40.*3.6,\
                        _d=.45e-6,_delta_m=100e-6,_tau=1.5,_mhi=1.5*60,_mci=1.5*60,_thi=thi_N,_tci=20,_shi=42.,_sc=0.)
        verify_fig10.get_mesh_setup(Nx)
        verify_fig10.get_co_current(Nx);
        Jv_simulate_N[T_i] = sum(verify_fig10.Jv_flux_condensed_water_along)/Nx
        L1_norm.append((Jv_simulate_N[T_i] - Jv_exp_fig10_re[T_i])/Jv_exp_fig10_re[T_i])
    print "Single error: L1-norm: %r " %L1_norm
    L2_norm = sum(((Jv_simulate_N - Jv_exp_fig10_re)/Jv_exp_fig10_re)**2)**(.5)
    print "Total L2-norm: %r " %L2_norm


.. parsed-literal::

    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 80 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 7.51776013,  7.50969861,  7.50436784,  7.49932758,  7.49442873,
            7.48961736,  7.48486706,  7.48016268,  7.47549461,  7.47085635,
            7.46624325,  7.46165187,  7.45707958,  7.4525243 ,  7.44798439,
            7.4434585 ,  7.43894552,  7.43444451,  7.42995468,  7.42547537]) (kg/(m2 hr))
    Average permeate water: 7.47152 kg/hr/m2
    Average of sum(flux)/Nx: 7.47 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 70 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 4.77119876,  4.76719012,  4.76462205,  4.76220646,  4.75986471,
            4.75756833,  4.75530342,  4.75306199,  4.75083896,  4.74863093,
            4.74643543,  4.74425066,  4.74207524,  4.73990806,  4.73774827,
            4.73559516,  4.73344815,  4.73130673,  4.72917051,  4.72703912]) (kg/(m2 hr))
    Average permeate water: 4.74897 kg/hr/m2
    Average of sum(flux)/Nx: 4.75 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 2.86678731,  2.86491958,  2.86376334,  2.86268215,  2.86163716,
            2.86061433,  2.85960678,  2.85861057,  2.8576232 ,  2.85664298,
            2.8556687 ,  2.85469947,  2.8537346 ,  2.85277357,  2.85181594,
            2.85086136,  2.84990954,  2.84896025,  2.84801328,  2.84706845]) (kg/(m2 hr))
    Average permeate water: 2.85681 kg/hr/m2
    Average of sum(flux)/Nx: 2.86 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 50 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 1.59262557,  1.59183188,  1.59135622,  1.59091403,  1.59048794,
            1.59007168,  1.58966218,  1.58925767,  1.58885704,  1.58845953,
            1.58806461,  1.58767187,  1.58728101,  1.58689179,  1.58650402,
            1.58611754,  1.58573223,  1.58534797,  1.58496468,  1.58458228]) (kg/(m2 hr))
    Average permeate water: 1.58853 kg/hr/m2
    Average of sum(flux)/Nx: 1.59 
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 40 C
    Cooling water inlet at 20 C
    Hot feed at 1.50 LPM
    Cooling mass inlet: 1.50 LPM
    Air gap distance: 9.00 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 0.77501738,  0.77472888,  0.77455974,  0.77440315,  0.77425259,
            0.7741057 ,  0.77396133,  0.77381882,  0.77367775,  0.77353783,
            0.77339886,  0.7732607 ,  0.77312323,  0.77298636,  0.77285001,
            0.77271414,  0.77257868,  0.77244361,  0.77230888,  0.77217447]) (kg/(m2 hr))
    Average permeate water: 0.77356 kg/hr/m2
    Average of sum(flux)/Nx: 0.77 
    Single error: L1-norm: [0.13341686582697454, 0.497751783509988, 0.37986455420786325, 0.44394007764942811, 0.40635473880690781] 
    Total L2-norm: 0.87866998578090516 
    

Following are plotting for Alsaadi test.

.. code:: python

    plot_along_membrane(plot_title="Test_Alsaadi",plot_ylabel="Jv_along(kg/m2 hr)",x_values=test_alsaadi.x,\
                        y_values=test_alsaadi.Jv_flux_condensed_water_along ,y_limits="no")


.. parsed-literal::

    No specific restrict boundary.
    

.. parsed-literal::

    C:\WinPython-32bit-2.7.9.2\python-2.7.9\lib\site-packages\matplotlib\font_manager.py:1282: UserWarning: findfont: Font family [u'WenQuanYi Micro Hei'] not found. Falling back to Bitstream Vera Sans
      (prop.get_family(), self.defaultFamily[fontext]))
    


.. image:: output_33_2.png


.. code:: python

    plot_along_membrane(plot_title="Test1",plot_ylabel="deltaf_along(m)",x_values=test_alsaadi.x,\
                        y_values=test_alsaadi.delta_f_condensed_water_along ,y_limits=[0.,1e-4])



.. image:: output_34_0.png


.. code:: python

    print test_alsaadi.delta_f_condensed_water_along


.. parsed-literal::

    [  0.00000000e+00   1.48771586e-05   1.87327023e-05   2.14326637e-05
       2.35787610e-05   2.53883839e-05   2.69679465e-05   2.83784513e-05
       2.96585091e-05   3.08342897e-05   3.19244471e-05   3.29428019e-05
       3.38999127e-05   3.48040495e-05   3.56618253e-05   3.64786209e-05
       3.72588804e-05   3.80063216e-05   3.87240903e-05   3.94148745e-05]
    

.. code:: python

    plot_along_membrane(plot_title="Test1",plot_ylabel="Q(kJ/m2 hr)",x_values=test_alsaadi.x,\
                        y_values=test_alsaadi.Q_heattransfer_along ,y_limits="no")


.. parsed-literal::

    No specific restrict boundary.
    


.. image:: output_36_1.png


.. code:: python

    help(plot_temperature)


.. parsed-literal::

    Help on function plot_temperature in module __main__:
    
    plot_temperature(plot_title, plot_ylabel1, plot_ylabel2, plot_ylabel3, plot_ylabel4, plot_ylabel5, x_values, y_values1, y_values2, y_values3, y_values4, y_values5, y_limits)
        Pyplot to show T, mh, Jv, delta_f along the membrane.
        Params
        ------
        plot_title,plot_ylabel_N: strings, to describe Title and Y-Label of the fig plotted.
        x_values: x along the membrane. 1-D array of floats.
        y_values_N: (C) Temperature along the membrane. 1-D array of floats.
        y_limits: 1-D array with 2 elements which defines the lower/upper boundary of the plot y-value
                    i.e. y_limits = [lowerbound, upperbound]
                    if no concerns for boundary, input y_limits = "no boundary" instead.
    
    

.. code:: python

    plot_temperature("Test1","Thb","Thm","Tma","Tf","Tcb",test_alsaadi.x,\
                        test_alsaadi.Thb,test_alsaadi.Thm,test_alsaadi.Tma,test_alsaadi.Tf,test_alsaadi.Tcb,[19.5,80.5])



.. image:: output_38_0.png


.. code:: python

    print 

.. code:: python

    a = [1,2]
    def test_func_variables(a):
        if type(a) == str:
            print "a is a string."
        elif type(a)==int:
            print "a + 1 = %r " %(a+1)
        elif type(a) == list:
            print a[0],a[1]
        print "a is : %r." %a

.. code:: python

    test_func_variables([1,2])

.. code:: python

    b =1
    type(b)

.. code:: python

    # pyplot to show T, mh, Jv, delta_f along the membrane
    x = np.linspace(0,10,5)
    y_values = [10,20,5,30,45]
    plot_title = 'delta_f along'
    plot_ylabel = "T hb"
    pyplot.figure(figsize=(10,5))
    pyplot.title(plot_title);
    pyplot.xlabel('x(m)');
    pyplot.ylabel(plot_ylabel);
    #pyplot.ylim(0.0,1.0e-4)
    pyplot.scatter(x,y_values,marker='^')
    pyplot.plot(x,y_values,color='#003366',ls='-')

.. code:: python

    # pyplot to show Jv varies attributed to T(inlet)
