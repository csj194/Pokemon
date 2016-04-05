''' This is a program for the absorbtion of a 2 component mixture+solvent.
For simplicity we have assumed the solvent to be water, but as a class is created we can
change these values for obtaining solutions to our system.

Assumptions:
1)Density of solvent does not change even after adsorption of components
2)Operating temperature remains constant throughout the column
3)For mass transfer of gas to liquid, resitance only on liquid side, gas side resistance is negligible
4)For mass transfer of gas to liquid, resitance only on gas side, liquid side resistance is negligible 

Parameter Definitions:
Top=Operating temperature(in celsius)
Pop=Operating pressure(Pa)
Ar=area of cross section of Column(m^2)
H=Height of the column(m)
rhomlar=molar density of the solvent(kmol/m^3)
kl1=liquid side mass transfer co-efficient for 1st gas component(m^2/s)
kl2=liquid side mass transfer co-efficient for 2nd gas component(m^2/s)
kg=gas side mass transfer co-efficient for solvent(m^2/s)
a=specific interfacial area of column (m^2/m^3)
hc1=Henry's constant for 1st component(Pa/(kmol/m^3))
hc2=Henry's constant for 2nd component(Pa/(kmol/m^3))
antA,antB,antC=Antoine constants for solvent.
Psat=10^(A-B/T+C);Psat=saturation pressure(bar),T=Top (celsius)
gc10,gc20=gas side inlet molar flow rates of components(kmol/s)
gcl0=gas side inlet solvent(kmol/s)
lc10,lc20=liquid side inlet molar flow rates of components(kmol/s)
lcl0=liquid side inlet solvent(kmol/s)
lc11,lc21=liquid side outlet molar flow rate GUESS VALUES of components(kmol/s)
lcl1=liquid side outlet molar flow rate of solvent GUESS VALUES(kmol/s)


Conclusions:
1) There is very little change in water molar rates at temp less than 100 at 1 atm total pressure
2) On increasing the operating pressure more components get adsorbed.
3) On changing the specific interfacial area the the adsorbed components change
4) Weird Values come above the critical point and at about 350 C showing inaccuracy of Antoine Coefficients
'''

class GasAbsorption:
    R=8314.462 #J/kmol.K        
    Top=25.0#Centigrade
    Pop=101325.0 #Pascal
    Ar=1.0 #m^2 
    H=30.0 #m
    rhomolar=55.55 #kmol/m^3
    kl1=29.39e-7#m^2/s
    kl2= 17.77e-7#m^2/s
    kg= 1.367e-3#m^2/s
    a=330.0#1/m
    hc1=140.0#Pa/(kmol/m^3)
    hc2=130.0#Pa/(kmol/m^3)
    if Top<100:
        antA= 8.07131
        antB= 1730.63
        antC= 233.426
    if Top>100:
        antA= 8.14019
        antB= 1810.94
        antC= 244.485
    
    gc10=5.000000 #kmol/s
    gc20=5.000000 #kmol/s
    gcl0=3.0000000 #kmol/s
    lc10=1.0000000 #kmol/s
    lc20=1.0000000 #kmol/s
    lcl0=10.000000 #kmol/s
    lc11=1.0 #kmol/s
    lc21=1.0 #kmol/s
    lcl1=10.000 #kmol/s
    xinc=0.01 #kmol/s
    def Operation(self):
        import scipy.integrate as sci
        import scipy.optimize as sco
        import numpy as np
        import matplotlib.pyplot as plt
        
        plt.style.use('fivethirtyeight')
        
        
                
#Operating Conditions:--------------------------------------------------------------------------------------------------        
        R=self.R        
        Top=self.Top
        Pop=self.Pop
        Ar=self.Ar
        H=self.H
        rhomolar=self.rhomolar
        kl1=self.kl1
        kl2=self.kl2
        kg=self.kg
        a=self.a
        hc1=self.hc1
        hc2=self.hc2
        antA=self.antA
        antB=self.antB
        antC=self.antC
        gc10=self.gc10
        gc20=self.gc20
        gcl0=self.gcl0
        lc10=self.lc10
        lc20=self.lc20
        lcl0=self.lcl0
        lc11=self.lc11
        lc21=self.lc21
        lcl1=self.lcl1    
        xinc=self.xinc
        
        Gin0= np.array([gc10,gc20,gcl0])
        Lin0= np.array([lc10,lc20,lcl0])
        lc11=1.0
        lc21=1.0
        lcl1=10.0
        Lout= np.array([lc11,lc21,lcl1])
        kl1a=kl1*a
        kl2a=kl2*a
        kga=kg*a
        Psat=((10**(antA-(antB/(Top+antC)))*(10**5)))/760
        
        def set_diff_equat(molar_flow_rates,x,param):
            '''To create an array of differential equations which govern the mass transfer
            in the column'''            
            gc1,gc2,gcl,lc1,lc2,lcl=molar_flow_rates
            Pop,Ar,H,kl1a,kl2a,kga,hc1,hc2,Psat,rhomolar,R=param
            
            y1=(gc1/(gc1+gc2+gcl))
            y2=(gc2/(gc1+gc2+gcl))
            yl=(gcl/(gc1+gc2+gcl))
            
            x1=(lc1/(lc1+lc2+lcl))
            x2=(lc2/(lc1+lc2+lcl))
            
            dgc1dh=-kl1a*Ar*(((Pop*y1)/hc1)-(rhomolar*x1))
            dgc2dh=-kl2a*Ar*(((Pop*y2)/hc2)-(rhomolar*x2))
            dgcldh=(kga*Ar*(Psat-(yl*Pop)))/(R*(Top+273.16))
            
            dlc1dh=-kl1a*Ar*(((Pop*y1)/hc1)-(rhomolar*x1))
            dlc2dh=-kl2a*Ar*(((Pop*y2)/hc2)-(rhomolar*x2))
            dlcldh=(kga*Ar*(Psat-(yl*Pop)))/(R*(Top+273.16))
            
            derivs=[dgc1dh,dgc2dh,dgcldh,dlc1dh,dlc2dh,dlcldh]
            return (derivs)        
        
        
        def opt(Lout):      
            '''To find the initial solution to our initial 
            guess of liquid inlet flowrates and then calculating the error from our
            actual values'''
            xstop=H
            n=xstop/xinc
            x=np.arange(0.,xstop,xinc)
            param=[Pop,Ar,H,kl1a,kl2a,kga,hc1,hc2,Psat,rhomolar,R]
            molar_flow_rates0= np.concatenate((Gin0,Lout),axis=0)
            Soln=sci.odeint(set_diff_equat,molar_flow_rates0,x,args=(param,))
            Lin1=Soln[n-1][3:6]
            err=np.square(np.subtract(Lin1,Lin0))
            sum_err=np.sum(err)
            return (sum_err)
        
        '''Obtaining the exact solution by minimising the function error=0'''        
        Lsoln=sco.fmin(opt,Lout)
                        
        xstop=H
        n=xstop/xinc
        no_iteration=int(n)
        x=np.arange(0.,xstop,xinc)
        param=[Pop,Ar,H,kl1a,kl2a,kga,hc1,hc2,Psat,rhomolar,R]        
        molar_flow_rates_soln= np.concatenate((Gin0,Lsoln),axis=0)
        Soln=sci.odeint(set_diff_equat,molar_flow_rates_soln,x,args=(param,))
                
        for i in range(no_iteration):
            gascomp1=Soln[:,0]
        for i in range(no_iteration):
            gascomp2=Soln[:,1]
        for i in range(no_iteration):
            gascomp3=Soln[:,2]
        for i in range(no_iteration):
            liqcomp1=Soln[:,3]
        for i in range(no_iteration):
            liqcomp2=Soln[:,4]
        for i in range(no_iteration):
            liqcomp3=Soln[:,5]
#Displaying Results--------------------------------------------------------------------        
        print"Gas side inlet conc of component 1(kmol/s):",gc10
        print"Gas side inlet conc of component 2(kmol/s):",gc20
        print"Gas side inlet conc of solvent(kmol/s):",gcl0
        print"Liquid side inlet conc of component 1(kmol/s):",lc10
        print"Liquid side inlet conc of component 2(kmol/s):",lc20
        print"Liquid side inlet conc of solvent(kmol/s):",lcl0
        print"Gas side outlet conc of component 1(kmol/s):",Soln[no_iteration-1,0]
        print"Gas side outlet conc of component 2(kmol/s):",Soln[no_iteration-1,1]
        print"Gas side outlet conc of solvent(kmol/s):",Soln[no_iteration-1,2]
        print"Liquid side outlet conc of component 1(kmol/s):",Lsoln[0]
        print"Liquid side outlet conc of component 2(kmol/s):",Lsoln[1]
        print"Liquid side outlet conc of solvent(kmol/s):",Lsoln[2]
                          
        s11=gc10+lc10
        s12=Soln[no_iteration-1,0]+Lsoln[0]
        p1=(abs(s12-s11)/s11)*100
        print "Error in component 1 Mass Balance(orig)",p1 
        s21=gc20+lc20
        s22=Soln[no_iteration-1,1]+Lsoln[1]
        p2=(abs(s22-s21)/s21)*100
        print "Error in component 2 Mass Balance(orig)",p2
        s31=gcl0+lcl0
        s32=Soln[no_iteration-1,2]+Lsoln[2]
        p3=(abs(s32-s31)/s31)*100
        print "Error in solvent Mass Balance(orig)",p3
#Checking Error in M.B. with step size--------------------------------------------------------------------        
        xinc1=[0.01,0.05,0.1,0.5,1,2,5]
        err_array=[0,0,0,0,0,0,0]
        for i in range(7):
            xstop1=H
            n1=xstop1/xinc1[i]
            no_iteration1=int(n1)
            x1=np.arange(0.,xstop1,xinc1[i])
            param1=[Pop,Ar,H,kl1a,kl2a,kga,hc1,hc2,Psat,rhomolar,R]        
            molar_flow_rates_soln1= np.concatenate((Gin0,Lsoln),axis=0)
            Soln1=sci.odeint(set_diff_equat,molar_flow_rates_soln1,x1,args=(param1,))
            s111=gc10+lc10
            s121=Soln1[no_iteration1-1,0]+Lsoln[0]
            p11=(abs(s121-s111)/s111)*100
            err_array[i]=p11
        
#Plotting Results--------------------------------------------------------------------        
        list_of_legends=['Gc1','Gc2','Gcl','Lc1','Lc2','Lcl']
        fig=plt.figure()
        
        ax=fig.add_subplot(111)
        bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
        t = ax.text(1,11, "Gas Phase Direction", ha="center", va="center", rotation=0,size=6,bbox=bbox_props)
        bb = t.get_bbox_patch()
        bb.set_boxstyle("rarrow", pad=0.6)
        bbox_props1 = dict(boxstyle="larrow,pad=0.3", fc="cyan", ec="b", lw=2)
        t1 = ax.text(4,11, "Liquid Phase Direction", ha="center", va="center", rotation=0,size=6,bbox=bbox_props1)
        bb1 = t1.get_bbox_patch()
        bb1.set_boxstyle("larrow", pad=0.6)
        ax.plot(x,gascomp1)
        ax.plot(x,gascomp2)
        ax.plot(x,gascomp3)
        ax.plot(x,liqcomp1)
        ax.plot(x,liqcomp2)
        ax.plot(x,liqcomp3)
        plt.legend(tuple(list_of_legends))
        fig.suptitle('Moles vs Height', fontsize=20)
        plt.xlabel('Height(m)', fontsize=16)
        plt.ylabel('moles(kmol/s)', fontsize=16)
        
        fig1=plt.figure()
        ay=fig1.add_subplot(111)
        ay.plot(xinc1,err_array)
        fig1.suptitle('%Error vs Step Size', fontsize=20)
        plt.xlabel('Step Size in Height', fontsize=16)
        plt.ylabel('%Error', fontsize=16)
        return fig
        
        
            
                
        
        
        
       

       

                  
            