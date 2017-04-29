
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
        T = T+273.15 #K
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
    H, for Heat Transfer a example: 75061.207172190465 (kJ/(m2 K hr))= 20850.335325608463 (W/(m2 K))
    miu: 0.0006553616086282997(Pa s)
    

.. code:: python

    # Computation for Jv
    def ComputationForJv(d,delta_m,tau,Tave,Thm,Tf,air_gap_width):
        '''
        For neat code in the Loop.
        '''
        Pv = VaporPressure(Tave+273.15) # water vapour pressure inside membrane
        C = CalculateC(d,delta_m,tau,Tave,Ptotal,Pv)
        A = CalculateA(Tave,Ptotal,air_gap_width)
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
                print "Tcb[-1]: %.2f C" %self.Tcb[-1]
                error_Tcb = abs(temp_move)
    
                alpha = 1./(error_i+1)*learnRate
                error_i += 1 
                self.Tcb[0] = self.Tcb[0] - temp_move*alpha
            #checking output
            print 'Difference of Tcb between steps is: %r' %error_Tcb
            print 
            print "Results are: "
            print "Jv along the length: %-5.3f (kg/(m2 hr))" %self.Jv_flux_condensed_water_along
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
    Jv along the length: array([ 7.63544621,  7.62735388,  7.62199118,  7.61691882,  7.61198772,
            7.607144  ,  7.60236124,  7.59762431,  7.5929236 ,  7.58825263,
            7.58360673,  7.57898248,  7.57437723,  7.56978893,  7.56521593,
            7.56065688,  7.55611066,  7.55157635,  7.54705316,  7.54254041]) (kg/(m2 hr))
    Average permeate water: 7.58891 kg/hr/m2
    Average of sum(flux)/Nx: 7.59 
    

.. code:: python

    # For testing wrap class-instance into 1 list.
    re_airgap = np.array([1.024,1.232,1.504,1.984,2.863,5.23]) # cm-1
    delta_a_Nvalues = 1./re_airgap/100 # m
    TEST = [] # List to store class-instance
    Jv_varies_delta_a = np.zeros(len(delta_a_Nvalues))
    for test_i,airgap_distance in enumerate(delta_a_Nvalues):
        ex_airgap = AGMD(_L=.215,_W=.165,_H=2e-3,_delta_a=airgap_distance,_delta_c=1.5e-3,_Kw=40.*3.6,\
                        _d=.45e-6,_delta_m=110e-6,_tau=1.33,_mhi=5.5*60,_mci=5.5*60,_thi=60.,_tci=20,_shi=32.,_sc=0.)
        ex_airgap.get_mesh_setup(Nx)
        ex_airgap.get_co_current(Nx);
        Jv_varies_delta_a[test_i] = sum(ex_airgap.Jv_flux_condensed_water_along)/Nx #Jv (average along membrane) for each instance
        TEST.append(ex_airgap)
        print test_i,TEST[test_i].Thm
    print Jv_varies_delta_a


.. parsed-literal::

    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 9.77 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 4.18311817,  4.1767363 ,  4.17259936,  4.16870159,  4.16492045,
            4.16121162,  4.15755339,  4.15393329,  4.15034344,  4.14677847,
            4.14323457,  4.1397089 ,  4.13619929,  4.13270404,  4.12922179,
            4.12575142,  4.12229203,  4.11884283,  4.11540318,  4.11197252]) (kg/(m2 hr))
    Average permeate water: 4.14733 kg/hr/m2
    Average of sum(flux)/Nx: 4.15 
    0 [ 59.06745194  59.05508385  59.04222126  59.02931701  59.01639843
      59.00347529  58.99055232  58.97763223  58.96471668  58.9518068
      58.9389034   58.92600704  58.91311817  58.90023712  58.88736416
      58.8744995   58.86164331  58.84879573  58.83595688  58.82312686]
    [ 4.14556133  0.          0.          0.          0.          0.        ]
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 8.12 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 4.62053003,  4.61251924,  4.60736495,  4.60251553,  4.59781506,
            4.59320712,  4.58866409,  4.5841701 ,  4.5797151 ,  4.57529227,
            4.57089673,  4.56652488,  4.56217395,  4.55784177,  4.55352661,
            4.54922704,  4.5449419 ,  4.54067021,  4.53641113,  4.53216394]) (kg/(m2 hr))
    Average permeate water: 4.57600 kg/hr/m2
    Average of sum(flux)/Nx: 4.57 
    1 [ 58.96693687  58.95343832  58.93930508  58.92511794  58.91091192
      58.89669957  58.88248696  58.86827754  58.85407347  58.83987619
      58.82568672  58.81150581  58.79733401  58.78317175  58.76901939
      58.75487717  58.74074534  58.72662405  58.71251347  58.69841372]
    [ 4.14556133  4.57380858  0.          0.          0.          0.        ]
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 6.65 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 5.09326141,  5.083235  ,  5.07682961,  5.07081146,  5.0649829 ,
            5.0592724 ,  5.05364493,  5.04808035,  5.04256591,  5.03709296,
            5.03165532,  5.0262484 ,  5.02086869,  5.01551343,  5.01018041,
            5.00486784,  4.99957422,  4.9942983 ,  4.98903903,  4.98379548]) (kg/(m2 hr))
    Average permeate water: 5.03800 kg/hr/m2
    Average of sum(flux)/Nx: 5.04 
    2 [ 58.85691669  58.84223551  58.8267409   58.81117661  58.79558754
      58.77998972  58.76439095  58.74879565  58.73320657  58.71762557
      58.70205396  58.68649268  58.67094245  58.65540382  58.63987722
      58.62436301  58.60886145  58.59337279  58.57789721  58.56243486]
    [ 4.14556133  4.57380858  5.0352909   0.          0.          0.        ]
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 5.04 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 5.73288971,  5.71968091,  5.7113122 ,  5.70346263,  5.69586794,
            5.68843253,  5.68110952,  5.67387197,  5.66670284,  5.65959051,
            5.65252672,  5.64550532,  5.6385216 ,  5.63157186,  5.62465313,
            5.61776301,  5.61089948,  5.60406087,  5.59724574,  5.59045287]) (kg/(m2 hr))
    Average permeate water: 5.66082 kg/hr/m2
    Average of sum(flux)/Nx: 5.66 
    3 [ 58.70477941  58.68856172  58.67123613  58.65381465  58.63635835
      58.61888907  58.60141742  58.58394941  58.56648881  58.54903812
      58.53159913  58.51417312  58.49676107  58.47936373  58.46198169
      58.44461541  58.42726529  58.40993161  58.39261465  58.37531462]
    [ 4.14556133  4.57380858  5.0352909   5.65730607  0.          0.        ]
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 3.49 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 6.5123666 ,  6.49448363,  6.48324646,  6.47272484,  6.46255595,
            6.45260843,  6.44281786,  6.43314721,  6.423573  ,  6.41407929,
            6.40465467,  6.39529067,  6.38598078,  6.37671992,  6.36750399,
            6.35832964,  6.3491941 ,  6.34009504,  6.33103047,  6.32199871]) (kg/(m2 hr))
    Average permeate water: 6.41581 kg/hr/m2
    Average of sum(flux)/Nx: 6.41 
    4 [ 58.51086172  58.49286132  58.47328857  58.45357803  58.43381628
      58.41403433  58.3942472   58.37446343  58.35468834  58.3349255
      58.31517743  58.29544596  58.27573248  58.25603803  58.23636346
      58.21670942  58.19707646  58.177465    58.15787541  58.13830798]
    [ 4.14556133  4.57380858  5.0352909   5.65730607  6.41112006  0.        ]
    Next step: Setup mesh please.
    Mesh setup done.
    This is co-current flow pattern.
    Hot water inlet at 60 C
    Cooling water inlet at 20 C
    Hot feed at 5.50 LPM
    Cooling mass inlet: 5.50 LPM
    Air gap distance: 1.91 mm
    
    The end. i = 19.
    Results are: 
    Jv along the length: array([ 7.53266162,  7.50693574,  7.49088617,  7.47588489,  7.46140353,
            7.44725071,  7.43333246,  7.41959477,  7.40600333,  7.39253482,
            7.37917256,  7.36590412,  7.35271997,  7.33961263,  7.32657607,
            7.31360535,  7.30069639,  7.28784573,  7.27505044,  7.26230799]) (kg/(m2 hr))
    Average permeate water: 7.39514 kg/hr/m2
    Average of sum(flux)/Nx: 7.39 
    5 [ 58.22252977  58.2022697   58.17952158  58.15655155  58.13349642
      58.11040521  58.08730165  58.06419924  58.04110641  58.0180288
      57.99497037  57.97193404  57.94892196  57.92593583  57.90297695
      57.88004637  57.85714492  57.8342733   57.81143207  57.78862169]
    [ 4.14556133  4.57380858  5.0352909   5.65730607  6.41112006  7.38849896]
    

.. code:: python

    # For plot Jv versus delta_a. See fig10, Banat 1998
    # Reuse plot(not for membrane along)
    plot_along_membrane(plot_title="Exam Jv vs. delta_a",plot_ylabel="J v_ave(kg/m2 hr)",x_values=re_airgap,\
                        y_values=Jv_varies_delta_a,y_limits="no")


.. parsed-literal::

    No specific restrict boundary.
    


.. image:: output_27_1.png


.. code:: python

    plot_along_membrane(plot_title="Test1",plot_ylabel="Jv_along(kg/m2 hr)",x_values=test1_varify.x,\
                        y_values=test1_varify.Jv_flux_condensed_water_along ,y_limits="no")


.. parsed-literal::

    No specific restrict boundary.
    


.. image:: output_28_1.png


.. code:: python

    plot_along_membrane(plot_title="Test1",plot_ylabel="deltaf_along(m)",x_values=test1_varify.x,\
                        y_values=test1_varify.delta_f_condensed_water_along ,y_limits=[0.,1e-4])



.. image:: output_29_0.png


.. code:: python

    plot_along_membrane(plot_title="Test1",plot_ylabel="Q(kJ/m2 hr)",x_values=test1_varify.x,\
                        y_values=test1_varify.Q_heattransfer_along ,y_limits="no")


.. parsed-literal::

    No specific restrict boundary.
    


.. image:: output_30_1.png


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

    plot_temperature("Test1","Thb","Thm","Tma","Tf","Tcb",test1_varify.x,\
                        test1_varify.Thb,test1_varify.Thm,test1_varify.Tma,test1_varify.Tf,test1_varify.Tcb,"no")


.. parsed-literal::

    No specific restrict boundary.
    


.. image:: output_32_1.png


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

    print delta_fN

.. code:: python

    # pyplot to show Jv varies attributed to T(inlet)
