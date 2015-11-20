from numpy import *
from pylab import *
from datetime import datetime, time, timedelta
import numpy as np
import console_colors as ccl
from scipy.io.netcdf import netcdf_file
from ShiftTimes import *
import os
import matplotlib.patches as patches
import matplotlib.transforms as transforms
import h5py

#from read_NewTable import tshck, tini_icme, tend_icme, tini_mc, tend_mc, n_icmes, MCsig
#from z_expansion_gulisano import z as z_exp

def flags2nan(VAR, FLAG):
        cond            = VAR < FLAG
        VAR             = np.array(VAR)
        VAR[~cond]      = np.nan
        return VAR


def date_to_utc(fecha):
    utc = datetime(1970, 1, 1, 0, 0, 0, 0)
    time = (fecha - utc).total_seconds()
    return time


def dates_from_omni(t):
    time = []
    n = len(t)
    for i in range(n):
        yyyy = t[i][0]
        mm = t[i][1]
        dd = t[i][2]
        HH = t[i][3]
        MM = t[i][4]
        SS = t[i][5]
        uSS = t[i][6]
        time += [datetime(yyyy, mm, dd, HH, MM, SS, uSS)]

    return time


def utc_from_omni(file):
    t = np.array(file.variables['time'].data)
    dates = dates_from_omni(t)
    n = len(dates)
    time = np.zeros(n)
    for i in range(n):
        time[i] = date_to_utc(dates[i])

    return time

def selecc_data(data, tshk):
    time    = data[0]       #[s] utc sec
    rate    = data[1]

    day      = 86400.        # [seg]
    utc      = datetime(1970, 1, 1, 0, 0, 0, 0)
    tshk_utc = (tshk - utc).total_seconds()

    ti      = tshk_utc - 10.*day        # [seg] utc
    tf      = tshk_utc + 30.*day
    cond    = (time > ti) & (time < tf)

    time    = (time[cond] - tshk_utc) / day     # [days] since shock
    rate    = rate[cond]
    return (time, rate)

def selecc_window(data, tini, tend):
    time    = data[0]       #[s] utc sec
    y       = data[1]

    day         = 86400.                # [seg]
    utc         = datetime(1970, 1, 1, 0, 0, 0, 0)
    tini_utc    = (tini - utc).total_seconds()      # [s] utc sec
    tend_utc    = (tend - utc).total_seconds()      # [s] utc sec

    ti      = tini_utc              # [seg] utc
    tf      = tend_utc
    cond    = (time > ti) & (time < tf)

    time    = (time[cond] - tini_utc) / day     # [days] since 'ti'
    y       = y[cond]
    return (time, y)


def enoughdata(var, fgap):
    n       = len(var)
    ngood   = len(find(~isnan(var)))
    fdata   = 1.*ngood/n        # fraccion de data sin gaps
    if fdata>=(1.-fgap):
        return True
    else:
        return False


def averages_and_std(n_icmes, t_shck, ti_icme, dTday, nbin, t_utc, VAR, fgap):
        day = 86400.
        nok=0; nbad=0
        adap = []
        for i in range(n_icmes):
            dT      = (ti_icme[i] - t_shck[i]).total_seconds()/day  # [day]
            if dT>dTday:
                dt      = dT/nbin
                t, var  = selecc_window(
                                [t_utc, VAR],
                                t_shck[i], ti_icme[i]
                                )
                if enoughdata(var, fgap):           # pido q haya mas del 80% NO sean gaps
                    adap    += [adaptar(nbin, dt, t, var)]
                    nok     +=1
                else:
                    continue
            else:
                print " i:%d ---> Este evento es muy chico!, dT/day:%g" % (i, dT)
                nbad +=1

        VAR_adap = zeros(nbin*nok).reshape(nok, nbin)
        for i in range(nok):
                VAR_adap[i,:] = adap[i][1]

        VAR_avrg = zeros(nbin)
        VAR_std = zeros(nbin)
        ndata = zeros(nbin)
        for i in range(nbin):
            cond = ~isnan(VAR_adap.T[i,:])
            ndata[i] = len(find(cond))      # nro de datos != flag
            VAR_avrg[i] = mean(VAR_adap.T[i,cond])  # promedio entre los valores q no tienen flag
            VAR_std[i] = std(VAR_adap.T[i,cond])    # std del mismo conjunto de datos

        tnorm   = adap[0][0]
        return [nok, nbad, tnorm, VAR_avrg, VAR_std, ndata]

def adaptar(n, dt, t, r):
    #n  = int(5./dt)        # nro de puntos en todo el intervalo de ploteo
    tt  = zeros(n)
    rr  = zeros(n)
    for i in range(n):
        tmin    = i*dt
        tmax    = (i+1.)*dt
        cond    = (t>tmin) & (t<tmax)
        tt[i]   = mean(t[cond])
        rr[i]   = mean(r[cond])
    return [tt/(n*dt), rr]


def adaptar(nwndw, dT, n, dt, t, r):
    #n  = int(5./dt)        # nro de puntos en todo el intervalo de ploteo
    tt  = zeros(n)
    rr  = zeros(n)
    _nbin_  = n/(1+nwndw[0]+nwndw[1])   # nro de bins en la sheath
    for i in range(n):
        tmin    = (i-nwndw[0]*_nbin_)*dt
        tmax    = tmin + dt
        cond    = (t>tmin) & (t<tmax)
        tt[i]   = mean(t[cond])#; print "tt:", t[i]; pause(1)
        rr[i]   = mean(r[cond])
    return [tt/dT, rr]          # tiempo normalizado x la duracion de la sheath


def adaptar_ii(nwndw, dT, n, dt, t, r, fgap):
    #n  = int(5./dt)            # nro de puntos en todo el intervalo de ploteo
    tt  = zeros(n)
    rr  = zeros(n)
    _nbin_  = n/(1+nwndw[0]+nwndw[1])   # nro de bins en la sheath/mc
    cc  = (t>0.) & (t<dT)       # intervalo de la sheath/mc
    #print " r[cc]: ", r[cc]
    if len(r[cc])==0:           # no hay data en esta ventana
        rr      = nan*ones(n)
        enough  = False
    else:
        enough  = enoughdata(r[cc], fgap)   # [bool] True si hay mas del 80% de data buena.
    if not(enough): rr  = nan*ones(n)   # si no hay suficiente data, este evento no aporta
    for i in range(n):
        tmin    = (i-nwndw[0]*_nbin_)*dt
        tmax    = tmin + dt
        cond    = (t>tmin) & (t<tmax)
        #tt[i]   = mean(t[cond])#; print "tt:", t[i]; pause(1) # bug
        tt[i]   = tmin + .5*dt     # bug corregido
        if enough:
            cc    = ~isnan(r[cond])     # no olvidemos filtrar los gaps
            rr[i] = mean(r[cond][cc])

    return enough, [tt/dT, rr]          # tiempo normalizado x la duracion de la sheath/mc/etc


def selecc_window_ii(nwndw, data, tini, tend):
    time    = data[0]       #[s] utc sec
    y   = data[1]

    day     = 86400.                # [seg]
    utc     = datetime(1970, 1, 1, 0, 0, 0, 0)
    tini_utc    = (tini - utc).total_seconds()      # [s] utc sec
    tend_utc    = (tend - utc).total_seconds()      # [s] utc sec

    dt  = tend_utc - tini_utc
    ti  = tini_utc - nwndw[0]*dt            # [seg] utc
    tf  = tend_utc + nwndw[1]*dt
    cond    = (time > ti) & (time < tf)

    time    = (time[cond] - tini_utc) / day     # [days] since 'ti'
    y   = y[cond]
    return (time, y)


def averages_and_std_ii(nwndw,
        SELECC, #MCsig, MCwant,
        n_icmes, tini, tend, dTday, nbin, t_utc, VAR):
        day = 86400.
        nok=0; nbad=0
        adap = []
        for i in range(n_icmes):
            dT      = (tend[i] - tini[i]).total_seconds()/day  # [day]
            if ((dT>dTday) & SELECC[i]):# (MCsig[i]>=MCwant)):
                dt      = dT*(1+nwndw[0]+nwndw[1])/nbin
                t, var  = selecc_window_ii(
                            nwndw,              # nro de veces hacia atras y adelante
                            [t_utc, VAR],
                            tini[i], tend[i]
                            )
                adap    += [adaptar(nwndw, dT, nbin, dt, t, var)]       # rebinea usando 'dt' como el ancho de nuevo bineo
                nok     +=1
            else:
                print " i:%d ---> Filtramos este evento!, dT/day:%g" % (i, dT)
                nbad +=1

        VAR_adap = zeros(nbin*nok).reshape(nok, nbin)
        for i in range(nok):
            VAR_adap[i,:] = adap[i][1]

        VAR_avrg = zeros(nbin)
        VAR_medi = zeros(nbin)
        VAR_std  = zeros(nbin)
        ndata    = zeros(nbin)
        for i in range(nbin):
            cond = ~isnan(VAR_adap.T[i,:])
            ndata[i] = len(find(cond))      # nro de datos != flag
            VAR_avrg[i] = mean(VAR_adap.T[i,cond])  # promedio entre los valores q no tienen flag
            VAR_medi[i] = median(VAR_adap.T[i,cond])# mediana entre los valores q no tienen flag
            VAR_std[i] = std(VAR_adap.T[i,cond])    # std del mismo conjunto de datos

        tnorm   = adap[0][0]
        return [nok, nbad, tnorm, VAR_avrg, VAR_medi, VAR_std, ndata]


def mvs_for_each_event(VAR_adap, nbin, nwndw, Enough):
    nok         = size(VAR_adap, axis=0)
    mvs         = np.zeros(nok)            # valores medios por cada evento
    binsPerTimeUnit = nbin/(1+nwndw[0]+nwndw[1])    # nro de bines por u. de tiempo
    start       = nwndw[0]*binsPerTimeUnit  # en este bin empieza la MC
    #print " ----> binsPerTimeUnit: ", binsPerTimeUnit
    #print " ----> nok:  ", nok
    #print " ----> VAR_adap.shape: ", VAR_adap.shape
    #print " ----> VAR_adap: \n", VAR_adap
    #raw_input()

    for i in range(nok):
        aux = VAR_adap[i, start:start+binsPerTimeUnit]  # (*)
        cc  = ~isnan(aux)                   # pick good-data only
        #if len(find(cc))>1:
        if Enough[i]:       # solo imprimo los q tienen *suficiente data*
            print ccl.G
            print "id %d/%d: "%(i+1, nok), aux[cc]
            print ccl.W
            mvs[i] = np.mean(aux[cc])
        else:
            mvs[i] = np.nan
    #(*): esta es la serie temporal (de esta variable) para el evento "i"
    pause(1)
    return mvs


def diff_dates(tend, tini):
    n       = len(tend)
    diffs   = np.nan*np.ones(n)
    for i in range(n):
        ok  = type(tend[i]) == type(tini[i]) == datetime    # ambos deben ser fechas!
        if ok:
            diffs[i] = (tend[i] - tini[i]).total_seconds()
        else:
            diffs[i] = np.nan

    return diffs    #[sec]


def write_variable(fout, varname, dims, var, datatype, comments):
    dummy           = fout.createVariable(varname, datatype, dims)
    dummy[:]        = var
    dummy.units     = comments


def calc_beta(Temp, Pcc, B):
    # Agarramos la definicion de OMNI, de:
    # http://omniweb.gsfc.nasa.gov/ftpbrowser/magnetopause/Reference.html
    # http://pamela.roma2.infn.it/index.php
    # Beta = [(4.16*10**-5 * Tp) + 5.34] * Np/B**2 (B in nT)
    #
    beta = ((4.16*10**-5 * Temp) + 5.34) * Pcc/B**2
    return beta


def thetacond(ThetaThres, ThetaSh):
    if ThetaThres<=0.:
        print ccl.Rn + ' ----> BAD WANG FILTER!!: ThetaThres<=0.'
        print ' ----> Saliendo...' + ccl.Rn
        raise SystemExit
        #return ones(len(ThetaSh), dtype=bool)
    else:
        return (ThetaSh > ThetaThres)


def wangflag(ThetaThres):
    if ThetaThres<0:
        return 'NaN'
    else:
        return str(ThetaThres)


def makefig(medVAR, avrVAR, stdVAR, nVAR, tnorm,
        SUBTITLE, YLIMS, YLAB, fname_fig):
    fig     = figure(1, figsize=(13, 6))
    ax      = fig.add_subplot(111)

    ax.plot(tnorm, avrVAR, 'o-', color='black', markersize=5, label='mean')
    ax.plot(tnorm, medVAR, 'o-', color='red', alpha=.5, markersize=5, markeredgecolor='none', label='median')
    inf     = avrVAR + stdVAR/np.sqrt(nVAR)
    sup     = avrVAR - stdVAR/np.sqrt(nVAR)
    ax.fill_between(tnorm, inf, sup, facecolor='gray', alpha=0.5)
    trans   = transforms.blended_transform_factory(
                ax.transData, ax.transAxes)
    rect1   = patches.Rectangle((0., 0.), width=1.0, height=1,
                transform=trans, color='blue',
                alpha=0.3)
    ax.add_patch(rect1)

    ax.legend(loc='upper right')
    ax.grid()
    ax.set_ylim(YLIMS)
    TITLE = SUBTITLE
    ax.set_title(TITLE)
    ax.set_xlabel('time normalized to MC passage time [1]', fontsize=14)
    ax.set_ylabel(YLAB, fontsize=20)
    savefig(fname_fig, format='png', dpi=180, bbox_inches='tight')
    close()


#--- chekea q el archivo no repita elementos de la 1ra columna
def check_redundancy(fname, name):
    f = open(fname, 'r')
    dummy   = {}
    for line in f:
        ll          = line.split(' ')
        varname     = ll[0]
        dummy[varname] = 0

    dummy_names = dummy.keys()
    dummy_set   = set(dummy_names)
    redundancy  = len(dummy_set)<len(dummy_names)

    overwriting = name in dummy_set

    if redundancy or overwriting:
        return True
    else:
        return False


class general:
    def __init__(self):
        self.name = 'name'


class events_mgr:
    def __init__(self, gral, FILTER, CUTS, bd, nBin, fgap, tb, z_exp):
        #self.fnames     = fnames
        self.data_name  = gral.data_name
        self.FILTER     = FILTER
        self.CUTS       = CUTS
        self.bd         = bd
        self.nBin       = nBin
        self.fgap       = fgap
        self.tb         = tb
        self.z_exp      = z_exp
        self.dir_plots  = gral.dirs['dir_plots']
        self.dir_ascii  = gral.dirs['dir_ascii']
        self.gral       = gral
        self._dirs_     = gral.dirs

        #self.f_sc       = netcdf_file(gral.fnames[gral.data_name], 'r')
        self.f_events   = netcdf_file(gral.fnames['table_richardson'], 'r')
        print " -------> archivos input leidos!"

        self.read_ace   = False # True: si ya lei los archivos input
        self.read_murdo = False # True: si ya lei los archivos input
        self.read_o7o6  = False # True: si ya lei los archivos input
        self.read_auger = False # True: si ya lei los archivos input
        self.data_name_ = str(self.data_name) # nombre de la data input inicial (*1)
        """
        (*1):   si despues cambia 'self.data_name', me voy a dar
                cuenta en la "linea" FLAG_001.
        """


    def run_all(self):
        #----- seleccion de eventos
        self.filter_events()
        print "\n ---> filtrado de eventos (n:%d): OK\n" % (self.n_SELECC)
        #----- load data y los shiftimes "omni"
        self.load_data_and_timeshift()
        #----- rebineo y promedios
        #self.rebine_and_avr()
        self.rebine()
        self.rebine_final()
        #----- hacer ploteos
        self.make_plots()
        #----- archivos "stuff"
        #self.build_params_file()


    def rebine(self):
        """
        rebineo de c/evento
        """
        nvars   = self.nvars #len(VARS)
        n_icmes = self.tb.n_icmes
        bd      = self.bd
        VARS    = self.VARS
        nbin    = self.nBin['total']
        nwndw   = [self.nBin['before'], self.nBin['after']]
        day     = 86400.

        #---- quiero una lista de los eventos-id q van a incluirse en c/promedio :-)
        IDs     = {}
        Enough  = {}
        nEnough = {}
        self.__ADAP__       = ADAP    = []   # conjunto de varios 'adap' (uno x c/variable)
        for varname in VARS.keys():
            IDs[varname]        = []
            Enough[varname]     = []
            nEnough[varname]    = 0     # contador

        # recorremos los eventos:
        nok=0; nbad=0;
        nnn     = 0     # nro de evento q pasan el filtro a-priori
        for i in range(n_icmes):
            #nok=0; nbad=0;
            ok=False
            try: #no todos los elementos de 'tend' son fechas (algunos eventos no tienen fecha definida)
                dT      = (bd.tend[i] - bd.tini[i]).total_seconds()/day  # [day]
                ok = True
            except:
                continue    # saltar al sgte evento 'i'

            ADAP += [ {} ] # agrego un diccionario a la lista
            #np.set_printoptions(4)         # nro de digitos a imprimir cuando use numpy.arrays
            if (ok & self.SELECC[i]):# (MCsig[i]>=MCwant)):  ---FILTRO--- (*1)
                nnn += 1
                print ccl.Gn + " id:%d ---> dT/day:%g" % (i, dT) + ccl.W
                print self.tb.tshck[i]
                nok +=1
                # recorremos las variables:
                for varname in VARS.keys():
                    dt      = dT*(1+nwndw[0]+nwndw[1])/nbin
                    t, var  = selecc_window_ii(
                                nwndw, #rango ploteo
                                [self.t_utc, VARS[varname]['value']],
                                bd.tini[i],
                                bd.tend[i]
                              )

                    if self.data_name=='McMurdo':
                        var = 100.*(var - self.rate_pre[i]) / self.rate_pre[i]

                    elif self.data_name=='Auger':
                        #print " ---> var.size: ", var.size
                        var = 100.*(var - self.rate_pre_Auger[i]) / self.rate_pre_Auger[i]

                    # rebinea usando 'dt' como el ancho de nuevo bineo
                    out       = adaptar_ii(nwndw, dT, nbin, dt, t, var, self.fgap)
                    enough    = out[0]       # True: data con menos de 100*'fgap'% de gap
                    Enough[varname]         += [ enough ]
                    ADAP[nok-1][varname]    = out[1]  # donde: out[1] = [tiempo, variable]

                    if enough:
                        IDs[varname]     += [i]
                        nEnough[varname] += 1

            else:
                print ccl.Rn + " id:%d ---> dT/day:%g" % (i, dT), " -->SELECC: ", self.SELECC[i], ccl.W
                nbad +=1

        print " ----> len.ADAP: %d" % len(ADAP)
        self.__nok__    = nok
        self.__nbad__   = nbad
        self.out = {}
        self.out['nok']     = nok
        self.out['nbad']    = nbad
        self.out['IDs']     = IDs
        self.out['nEnough'] = nEnough
        self.out['Enough']  = Enough


    def rebine_final(self):
        """
        rebineo de c/evento ... PARTE FINAL
        """
        nvars   = self.nvars #len(VARS)
        VARS    = self.VARS
        nbin    = self.nBin['total']
        nwndw   = [self.nBin['before'], self.nBin['after']]
        day     = 86400.
        ## salidas del 'self.rebine()'
        #Enough  = self.__Enough__
        #nEnough = self.__nEnough__
        ADAP    = self.__ADAP__
        Enough  = self.out['Enough']
        nEnough = self.out['nEnough']
        IDs     = self.out['IDs']
        nok     = self.out['nok']
        nbad    = self.out['nbad']

        stuff       = {} #[]
        #nok = len(ADAP)/nvars  # (*)
        # (*) la dim de 'ADAP' es 'nvars' por el nro de eventos q pasaro el filtro en (*1)

        # Hacemos un lugar para la data rebineada (posible uso post-analisis)
        if self.data_name==self.data_name_:
            self.rebined_data = {} # creamos el diccionario UNA sola vez

        for varname in VARS.keys():
            print ccl.On + " -------> procesando: %s" % VARS[varname]['label'] #+ "  (%d/%d)"%(j+1,nvars)
            print " nEnough/nok/(nok+nbad): %d/%d/%d " % (nEnough[varname], nok, nok+nbad) + ccl.W
            VAR_adap = np.zeros((nok, nbin))    # perfiles rebineados (*)
            # (*): uno de estos por variable
            # recorro los 'nok' eventos q pasaron el filtro de arriba:
            for i in range(nok):
                VAR_adap[i,:] = ADAP[i][varname][1] # valores rebineados de la variable "j" para el evento "i"

            self.rebined_data[varname] = VAR_adap

            # valores medios de esta variable para c/evento
            avrVAR_adap = mvs_for_each_event(VAR_adap, nbin, nwndw, Enough[varname])
            #print " ---> (%d/%d) avrVAR_adap[]: \n" % (j+1,nvars), avrVAR_adap
            print " ---> (%s) avrVAR_adap[]: \n" % varname, avrVAR_adap

            VAR_avrg        = np.zeros(nbin)
            VAR_avrgNorm    = np.zeros(nbin)
            VAR_medi        = np.zeros(nbin)
            VAR_std         = np.zeros(nbin)
            ndata           = np.zeros(nbin)

            for i in range(nbin):
                cond = ~np.isnan(VAR_adap.T[i,:])  # filtro eventos q no aportan data en este bin
                ndata[i] = len(find(cond))      # nro de datos != nan
                VAR_avrg[i] = np.mean(VAR_adap.T[i,cond])  # promedio entre los valores q no tienen flag
                VAR_avrgNorm[i] = np.mean(VAR_adap.T[i,cond]/avrVAR_adap[cond])
                VAR_medi[i] = np.median(VAR_adap.T[i,cond])# mediana entre los valores q no tienen flag
                VAR_std[i] = np.std(VAR_adap.T[i,cond])    # std del mismo conjunto de datos
                #--- calculo perfil normalizado por c/variable
                #ii = nwndw[0]*binsPerTimeUnit
                #AvrInWndw = mean(VAR_avrg[ii:ii+binsPerTimeUnit])
            first_varname = ADAP[0].keys()[0]
            tnorm   = ADAP[0][first_varname][0] # tiempo del primer evento (0), usando la 1ra variable
            stuff[varname] = [VAR_avrg, VAR_medi, VAR_std, ndata, avrVAR_adap]
            # NOTA: chekar q 'ADAP[j][varname][0]' sea igual para TODOS los
            #       eventos 'j', y para TODOS los 'varname'.

        self.out['dVARS']    = stuff
        self.out['tnorm']    = tnorm #OUT['dVARS'][first_varname][2] # deberia ser =tnorm


    def load_data_and_timeshift(self):
        if self.data_name=='ACE':
            if not(self.read_ace):
                self.load_data_ACE()
                self.read_ace = True # True: ya lei los archivos input

        elif self.data_name=='Auger':
            if not(self.read_auger):
                self.load_data_Auger()
                self.read_auger = True # True: ya lei los archivos input

        elif self.data_name=='McMurdo':
            if not(self.read_murdo):
                self.load_data_McMurdo()
                self.read_murdo = True # True: ya lei los archivos input

        elif self.data_name=='ACE_o7o6':
            if not(self.read_o7o6):
                self.load_data_o7o6()
                self.read_o7o6 = True # True: ya lei los archivos input

        else:
            print " --------> BAD 'self.data_name'!!!"
            print " ---> self.read_ace: ", self.read_ace
            print " ---> self.read_murdo: ", self.read_murdo
            print " ---> self.read_auger: ", self.read_auger
            print " exiting.... "
            raise SystemExit


    def load_data_o7o6(self):
        tb          = self.tb
        nBin        = self.nBin
        bd          = self.bd
        day         = 86400.
        fname_inp   = self.gral.fnames[self.data_name]
        self.f_sc   = netcdf_file(fname_inp, 'r')
        print " leyendo tiempo..."
        t_utc   = utc_from_omni(self.f_sc)
        print " Ready."

        #++++++++++++++++++++ CORRECCION DE BORDES +++++++++++++++++++++++++++++
        # IMPORTANTE:
        # El shift es necesario, pero ya lo hice en la corrida del
        # primer 'self.data_name'. Si lo hago aqui, estaria hacia
        # doble shift-time.
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        o7o6    = self.f_sc.variables['O7toO6'].data.copy()
        print " -------> variables leidas!"
        #------------------------------------ VARIABLES
        self.t_utc  = t_utc
        self.VARS = VARS = {}
        # variable, nombre archivo, limite vertical, ylabel
        #VARS += [[o7o6, 'o7o6', [0.0, 1.5], 'O7/O6 [1]']]
        VARS['o7o6'] = {
            'value' : o7o6,
            'lims'  : [0.0, 1.5],
            'label' : 'O7/O6 [1]'
        }

        self.nvars  = len(VARS.keys())
        #---------
        self.aux = aux = {}
        aux['SELECC']    = self.SELECC


    def load_data_Auger(self):
        tb          = self.tb
        nBin        = self.nBin
        bd          = self.bd
        day         = 86400.
        fname_inp   = self.gral.fnames[self.data_name]
        #data_murdo  = np.loadtxt(fname_inp)
        f5          = h5py.File(fname_inp, 'r')
        self.t_utc  = t_utc = f5['auger/time_seg_utc'][...].copy() #data_murdo[:,0]
        CRs         = f5['auger/sc_wAoP_wPres'][...].copy() #data_murdo[:,1]
        print " -------> variables leidas!"

        self.VARS   = VARS = {} #[]
        VARS['CRs'] = {
            'value' : CRs,
            'lims'  : [-1.0, 1.0],
            'label' : 'Auger rate [%]'
        }
        self.nvars  = len(VARS.keys())
        #---------

        self.aux = aux = {}
        aux['SELECC']    = self.SELECC


    def load_data_McMurdo(self):
        tb          = self.tb
        nBin        = self.nBin
        bd          = self.bd
        day         = 86400.
        #HOME        = os.environ['HOME']
        #fname_inp_murdo = '%s/actividad_solar/neutron_monitors/mcmurdo/mcmurdo_utc_correg.dat' % HOME
        fname_inp   = self.gral.fnames[self.data_name]
        data_murdo  = np.loadtxt(fname_inp)
        self.t_utc  = t_utc = data_murdo[:,0]
        CRs         = data_murdo[:,1]
        print " -------> variables leidas!"

        self.VARS   = VARS = {} #[]
        #VARS        += [[CRs, 'CRs', [-8.0, 1.0], 'mcmurdo rate [%]']]
        VARS['CRs'] = {
            'value' : CRs,
            'lims'  : [-8.0, 1.0],
            'label' : 'mcmurdo rate [%]'
        }
        self.nvars  = len(VARS.keys())
        #---------

        self.aux = aux = {}
        aux['SELECC']    = self.SELECC


    def load_data_ACE(self):
        tb          = self.tb
        nBin        = self.nBin
        bd          = self.bd
        day         = 86400.

        #----------------------------------------------------------
        self.f_sc       = netcdf_file(self.gral.fnames[self.gral.data_name], 'r')
        print " leyendo tiempo..."
        t_utc   = utc_from_omni(self.f_sc)
        print " Ready."

        #++++++++++++++++++++ CORRECCION DE BORDES +++++++++++++++++++++++++++++
        # IMPORTANTE:
        # Solo valido para los "63 eventos" (MCflag='2', y visibles en ACE)
        # NOTA: dan saltos de shock mas marcados con True.
        if self.FILTER['CorrShift']:
            ShiftCorrection(ShiftDts, tb.tshck)
            ShiftCorrection(ShiftDts, tb.tini_icme)
            ShiftCorrection(ShiftDts, tb.tend_icme)
            ShiftCorrection(ShiftDts, tb.tini_mc)
            ShiftCorrection(ShiftDts, tb.tend_mc)
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        B       = self.f_sc.variables['Bmag'].data.copy()
        Vsw     = self.f_sc.variables['Vp'].data.copy()
        Temp    = self.f_sc.variables['Tp'].data.copy()
        Pcc     = self.f_sc.variables['Np'].data.copy()
        rmsB    = self.f_sc.variables['dBrms'].data.copy()
        alphar  = self.f_sc.variables['Alpha_ratio'].data.copy()
        beta    = calc_beta(Temp, Pcc, B)
        rmsBoB  = rmsB/B
        print " -------> variables leidas!"
        #------------------------------------ VARIABLES
        self.t_utc  = t_utc
        #self.VARS = VARS = []
        self.VARS = VARS = {}
        # variable, nombre archivo, limite vertical, ylabel
        VARS['B'] = {
            'value' : B,
            'lims'  : [5., 18.],
            'label' : 'B [nT]'
        }
        VARS['V'] = {
            'value' : Vsw,
            'lims'  : [380., 650.],
            'label' : 'Vsw [km/s]'
        }
        VARS['rmsBoB'] = {
            'value' : rmsBoB,
            'lims'  : [0.01, 0.2],
            'label' : 'rms($\hat B$/|B|) [1]'
        }
        VARS['rmsB'] = {
            'value' : rmsB,
            'lims'  : [0.05, 2.0],
            'label' : 'rms($\hat B$) [nT]'
        }
        VARS['beta'] = {
            'value' : beta,
            'lims'  : [0.001, 5.],
            'label' : '$\\beta$ [1]'
        }
        VARS['Pcc'] = {
            'value' : Pcc,
            'lims'  : [2, 17.],
            'label' : 'proton density [#/cc]'
        }
        VARS['Temp'] = {
            'value' : Temp,
            'lims'  : [1e4, 4e5],
            'label' : 'Temp [K]'
        }
        VARS['AlphaRatio'] = {
            'value' : alphar,
            'lims'  : [1e-3, 0.1],
            'label' : 'alpha ratio [1]'
        }

        self.nvars = len(VARS.keys())
        #---------
        self.aux = aux = {}
        aux['SELECC']    = self.SELECC
        #---- SALIDA:
        #self.VARS   = VARS
        #self.out    = out
        #self.aux    = aux


    def make_plots(self):
        """
        #---- generar figuras y asciis de los perfiles promedio/mediana
        """
        nBin        = self.nBin
        fgap        = self.fgap
        MCwant      = self.FILTER['MCwant']

        #dt_mc       = self.aux['dt_mc']
        #dt_sh       = self.aux['dt_sh']
        ThetaThres  = self.CUTS['ThetaThres']
        if self.FILTER['vsw_filter']:
            v_lo, v_hi = self.CUTS['v_lo'], self.CUTS['v_hi']
        else:
            v_lo, v_hi = 0.0, 0.0   #estos valores significan q no hay filtro

        if self.FILTER['z_filter_on']:
            z_lo, z_hi = self.CUTS['z_lo'], self.CUTS['z_hi']
        else:
            z_lo, z_hi = 0.0, 0.0

        if self.FILTER['B_filter']:
            B_lo, B_hi = self.CUTS['B_lo'], self.CUTS['B_hi']
        else:
            B_lo, B_hi = 0.0, 0.0   #estos valores significan q no hay filtro

        if self.FILTER['filter_dR.icme']:
            dR_lo, dR_hi = self.CUTS['dR_lo'], self.CUTS['dR_hi']
        else:
            dR_lo, dR_hi = 0.0, 0.0   #estos valores significan q no hay filtro

        nbin        = (1+nBin['before']+nBin['after'])*nBin['bins_per_utime']  # [1] nro de bines q quiero en mi perfil promedio


        #-------------------- prefijos:
        # prefijo para filtro Wang:
        #WangFlag = wangflag(ThetaThres) #'NaN' #wangflag(ThetaThres)
        if self.FILTER['wang']:
            WangFlag = str(ThetaThres)
        else:
            WangFlag = 'NaN'

        # prefijo gral para los nombres de los graficos:
        if self.FILTER['CorrShift']:
            prexShift =  'wShiftCorr'
        else:
            prexShift = 'woShiftCorr'

        #-------------------------------
        # nombres genericos...
        DIR_FIGS    = '%s/MCflag%s/%s' % (self.dir_plots, MCwant['alias'], prexShift)
        DIR_FIGS    += '/' + self._dirs_['suffix']
        DIR_ASCII   = '%s/MCflag%s/%s' % (self.dir_ascii, MCwant['alias'], prexShift)
        DIR_ASCII   += '/' + self._dirs_['suffix']
        os.system('mkdir -p %s' % DIR_FIGS)     # si no existe, lo creamos
        os.system('mkdir -p %s' % DIR_ASCII)    # (bis)
        print ccl.On + " -------> creando: %s" % DIR_FIGS + ccl.W
        print ccl.On + " -------> creando: %s" % DIR_ASCII + ccl.W

        FNAMEs = 'MCflag%s_%dbefore.%dafter_fgap%1.1f' % (MCwant['alias'], nBin['before'], nBin['after'], fgap)
        FNAMEs += '_Wang%s' % (WangFlag)
        if self.FILTER['vsw_filter']:       FNAMEs += '_vlo.%03.1f.vhi.%04.1f' % (v_lo, v_hi)
        if self.FILTER['z_filter_on']:      FNAMEs += '_zlo.%2.2f.zhi.%2.2f' % (z_lo, z_hi)
        if self.FILTER['B_filter']:         FNAMEs += '_Blo.%2.2f.Bhi.%2.2f' % (B_lo, B_hi)
        if self.FILTER['filter_dR.icme']:   FNAMEs += '_dRlo.%2.2f.dRhi.%2.2f' % (dR_lo, dR_hi)

        FNAME_ASCII = '%s/%s' % (DIR_ASCII, FNAMEs)
        FNAME_FIGS  = '%s/%s' % (DIR_FIGS, FNAMEs)

        fname_nro   = DIR_ASCII+'/'+'n.events_'+FNAMEs+'.txt'
        #'w': write mode #'a': append mode
        #---FLAG_001
        if self.data_name==self.data_name_:
            fnro = open(fname_nro, 'w')
        else:
            fnro = open(fname_nro, 'a') # si uso otra data input, voy anotando el nro
                                        # de eventos al final del archivo 'fname_nro'

        #--------------------------------------------------------------------------------
        nvars = len(self.VARS)
        #for i in range(nvars):
        for varname in self.VARS.keys():
            fname_fig = '%s_%s.png' % (FNAME_FIGS, varname) #self.VARS[i][1])
            print ccl.Rn+ " ------> %s" % fname_fig
            #varname = self.VARS[i][1]
            ylims   = self.VARS[varname]['lims'] #self.VARS[i][2]
            ylabel  = self.VARS[varname]['label'] #self.VARS[i][3]
            average = self.out['dVARS'][varname][0]
            mediana = self.out['dVARS'][varname][1] #self.out['dVARS'][i][4]
            std_err = self.out['dVARS'][varname][2]
            nValues = self.out['dVARS'][varname][3] # number of values aporting to each data bin
            #binsPerTimeUnit = nbin  #nbin/(1+nbefore+nafter)
            N_selec = self.out['nok'] #self.out['dVARS'][varname][0]
            N_final = self.out['nEnough'][varname] #nEnough[i]

            SUBTITLE = '# of selected events: %d \n\
                events w/80%% of data: %d \n\
                bins per time unit: %d \n\
                MCflag: %s \n\
                WangFlag: %s' % (N_selec, N_final, nBin['bins_per_utime'], MCwant['alias'], WangFlag)

            makefig(mediana, average, std_err, nValues, self.out['tnorm'], SUBTITLE,
                    ylims, ylabel, fname_fig)

            fdataout = '%s_%s.txt' % (FNAME_ASCII, varname) #self.VARS[i][1])
            dataout = np.array([self.out['tnorm'] , mediana, average, std_err, nValues])
            print " ------> %s\n" % fdataout + ccl.W
            np.savetxt(fdataout, dataout.T, fmt='%12.5f')

            #-------- grabamos nro de eventos selecc para esta variable
            line = '%s %d %d\n' % (varname, N_final, N_selec)
            fnro.write(line)

        print ccl.Rn + " --> nro de eventos seleccionados: " + fname_nro + ccl.W
        fnro.close()

        #--- salidas (a parte de los .png)
        self.DIR_ASCII  = DIR_ASCII
        self.FNAMEs     = FNAMEs


    def build_params_file(self):
        """
        #---- construye archivo q contiene cosas de los eventos seleccionados:
        # - valores medios de los observables (B, Vsw, Temp, beta, etc)
        # - los IDs de los eventos
        # - duracion de los MCs y las sheaths
        """
        DIR_ASCII   = self.DIR_ASCII
        FNAMEs      = self.FNAMEs
        #---------------------------------------------- begin: NC_FILE
        print "\n**************************************** begin: NC_FILE"
        #------- generamos registro de id's de los
        #        eventos q entraron en los promedios.
        #        Nota: un registro por variable.
        fname_out   = DIR_ASCII+'/'+'_stuff_'+FNAMEs+'.nc' #'./test.nc'
        #---FLAG_001
        if self.data_name==self.data_name_:
            fout        = netcdf_file(fname_out, 'w')
        else:
            fout        = netcdf_file(fname_out, 'a')
            # modo 'a': si uso otra data input, voy anotando el nro
            # de eventos al final del archivo 'fname_out'

        print "\n ----> generando: %s\n" % fname_out

        IDs = self.out['IDs']
        #for i in range(len(self.VARS)):
        for varname in self.VARS.keys():
            #varname  = self.VARS[i][1]
            print " ----> " + varname
            n_events = len(IDs[varname])
            dimname  = 'nevents_'+varname
            fout.createDimension(dimname, n_events)
            print " n_events: ", n_events
            prom     = self.out['dVARS'][varname][4]
            cc       = np.isnan(prom)
            print " nprom: ", prom.size
            prom     = prom[~cc]
            print " nprom: ", prom.size
            dims     = (dimname,)
            write_variable(fout, varname, dims,
                    prom, 'd', 'average_values per event')
            #---------- IDs de esta variable
            ids      = map(int, IDs[varname])
            vname    = 'IDs_'+varname
            write_variable(fout, vname, dims, ids, 'i',
                    'event IDs that enter in this parameter average')
            #---------- duracion de la estructura
            dtsh     = np.zeros(len(ids))
            dtmc     = np.zeros(len(ids))
            for i in range(len(ids)):
                id      = ids[i]
                dtsh[i] = self.dt_sh[id]
                dtmc[i] = self.dt_mc[id]

            vname    = 'dt_sheath_'+varname
            write_variable(fout, vname, dims, dtsh, 'd', '[days]')
            vname    = 'dt_mc_'+varname
            write_variable(fout, vname, dims, dtmc, 'd', '[days]')

        fout.close()
        print "**************************************** end: NC_FILE"
        #---------------------------------------------- end: NC_FILE


    def filter_events(self):
        tb              = self.tb
        FILTER          = self.FILTER
        ThetaThres      = self.CUTS['ThetaThres']
        dTday           = self.CUTS['dTday']
        day             = 86400.
        AU_o_km         = 1./(150.0e6)
        sec_o_day       = 86400.
        #------------------------------------ EVENTS's PARAMETERS
        #MCsig  = array(f_events.variables['MC_sig'].data)# 2,1,0: MC, rotation, irregular
        #Vnsh   = array(f_events.variables['wang_Vsh'].data) # veloc normal del shock
        ThetaSh     = np.array(self.f_events.variables['wang_theta_shock'].data) # orientacion de la normal del shock
        #i_V         = np.array(self.f_events.variables['i_V'].data) # velocidad de icme
        i_V         = np.array(self.f_events.variables['mc_V'].data) # velocidad de icme
        #i_B         = np.array(self.f_events.variables['i_B'].data) # B del icme
        i_B         = np.array(self.f_events.variables['mc_B'].data) # B del icme
        #i_dt        = np.array(self.f_events.variables['i_dt'].data) # B del icme
        i_dt        = np.array(self.f_events.variables['mc_dt'].data) # B del icme
        i_dR            = i_dt*(i_V*AU_o_km*sec_o_day)
        self.rate_pre   = self.f_events.variables['rate_pre'].data.copy()
        self.rate_pre_Auger   = self.f_events.variables['rate_pre_Auger'].data.copy()
        self.Afd        = self.f_events.variables['A_FD'].data.copy()
        #------------------------------------

        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #++++++++++++++++++ begin: SELECCION DE EVENTOS ++++++++++++++++++++++
        #------- fechas
        BETW1998_2006   = np.ones(tb.n_icmes, dtype=bool)
        try: # si existe el filtro de fechas, q lo haga
            if FILTER['choose_1998-2006']: 
                for i in range(307, tb.n_icmes)+range(0, 26):
                        BETW1998_2006[i]=False # 'False' para excluir eventos

        except: # sino, aplica el filtro por defecto (backwards compatibility)
            for i in range(307, tb.n_icmes)+range(0, 26):
                    BETW1998_2006[i]=False # 'False' para excluir eventos

        #------- seleccionamos MCs con label-de-catalogo (lepping=2, etc)
        MC_FLAG = np.ones(tb.n_icmes, dtype=bool)
        for i in range(tb.n_icmes):
                MC_FLAG[i]  = tb.MCsig[i] in FILTER['MCwant']['flags']

        #------- excluimos eventos de 2MCs
        EVENTS_with_2MCs= (26, 148, 259, 295)
        MCmultiple      = FILTER['Mcmultiple'] #False #True para incluir eventos multi-MC
        MCmulti         = np.ones(tb.n_icmes, dtype=bool)  # False para eventos multi-MC (SI, escribi bien)
        if(~FILTER['Mcmultiple']):
            for i in EVENTS_with_2MCs:
                MCmulti[i] &= False

        #------- orientacion del shock (catalogo Wang)
        if FILTER['wang']:
            ThetaCond   = thetacond(ThetaThres, ThetaSh)

        #------- duration of sheaths
        self.dt_mc      = diff_dates(tb.tend_mc, tb.tini_mc)/day     # [day]
        self.dt_sh      = diff_dates(tb.tini_mc, tb.tshck)/day     # [day]
        dt              = diff_dates(self.bd.tend, self.bd.tini)/day
        DURATION        = dt > dTday    # sheaths>0

        #------- speed of icmes
        if FILTER['vsw_filter']:
            v_lo        = self.CUTS['v_lo']
            v_hi        = self.CUTS['v_hi']
            SpeedCond   = (i_V>=v_lo) & (i_V<v_hi)

        #------- z expansion (a. gulisano)
        z_exp   = self.z_exp
        if FILTER['z_filter_on']:
            z_lo    = self.CUTS['z_lo']
            z_hi    = self.CUTS['z_hi']
            z_cond  = (z_exp>=z_lo) & (z_exp<z_hi)

        #------- <B> of icmes
        if FILTER['B_filter']:
            B_lo        = self.CUTS['B_lo']
            B_hi        = self.CUTS['B_hi']
            BfieldCond  = (i_B>=B_lo) & (i_B<B_hi)

        #------- size of icmes
        if FILTER['filter_dR.icme']:
            dR_lo        = self.CUTS['dR_lo']
            dR_hi        = self.CUTS['dR_hi']
            """print " ---> i_dR: \n", i_dR
            print " ---> i_dt: \n", i_dt
            raw_input()"""
            dRicmeCond   = (i_dR>=dR_lo) & (i_dR<dR_hi)


        #------- filtro total
        SELECC  = np.ones(tb.n_icmes, dtype=bool)
        SELECC  &= BETW1998_2006    # nos mantenemos en este periodo de anios
        SELECC  &= MCmulti          # nubes multiples
        SELECC  &= MC_FLAG          # catalogo de nubes
        SELECC  &= DURATION         # no queremos sheaths q duran 1hr xq solo aportan ruido
        if FILTER['wang']:           SELECC &= ThetaCond # cerca a 180 es nariz del shock
        if FILTER['vsw_filter']:     SELECC &= SpeedCond
        if FILTER['z_filter_on']:    SELECC &= z_cond
        if FILTER['B_filter']:       SELECC &= BfieldCond
        if FILTER['filter_dR.icme']: SELECC &= dRicmeCond

        self.SELECC     = SELECC
        self.n_SELECC   = len(find(SELECC))
        #+++++++++++++++++ end: SELECCION DE EVENTOS +++++++++++++++++++++++++
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if self.n_SELECC<=0:
            print ccl.Rn + "\n --------> FATAL ERROR!!!: self.n_SELECC=<0"
            print " exiting....... \n" + ccl.W
            raise SystemExit


##
