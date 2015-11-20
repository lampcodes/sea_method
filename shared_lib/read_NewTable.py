from pylab import *
from numpy import *
from datetime import datetime, timedelta
import os

HOME        = os.environ['HOME']
fname_rich	= '%s/ASOC_ICME-FD/icmes_richardson/RichardsonList.csv' % HOME
print " leyendo tabla Richardson: %s" % fname_rich
frich		= open(fname_rich, 'r')
print " archivo leido."

ll 	= []
n 	= 0
for line in frich:
	ll 	+= [line.split(',')]
	n +=1

print " lineas leidas: %d" % n

tshck 		= []
tini_icme	= []
tend_icme	= []
tini_mc		= []
tend_mc		= []
Qicme		= []
MCsig		= []
Dst		= []
for i in range(1,n):
	#------ fecha shock
	tshck += [datetime.strptime(ll[i][1][1:20], "%Y-%m-%d %H:%M:%S")]
	#------ fecha ini icme
	ss	= ll[i][2][1:11].split()	# string de la fecha ini-icme
	HH	= int(ss[1][0:2])
	MM	= int(ss[1][2:4])
	mm	= int(ss[0].split('/')[0])
	dd	= int(ss[0].split('/')[1])
	if mm==tshck[i-1].month:
		yyyy = tshck[i-1].year
	else:
		yyyy = tshck[i-1].year + 1
	tini_icme += [datetime(yyyy, mm, dd, HH, MM)]
	#------ fecha fin icme
	ss      = ll[i][3][1:11].split()
	HH      = int(ss[1][0:2])
	MM      = int(ss[1][2:4])
	mm      = int(ss[0].split('/')[0])
	dd      = int(ss[0].split('/')[1])
	if mm==tshck[i-1].month:
		yyyy = tshck[i-1].year
	else:
		if tshck[i-1].month==12:
			yyyy = tshck[i-1].year + 1

	tend_icme += [datetime(yyyy, mm, dd, HH, MM)]
	#------ fechas MCs
	if ll[i][6]=='':
		tini_mc += [nan]
		tend_mc += [nan]
	else:
		hrs_ini	= int(ll[i][6])			# col6 es inicio del MC
		dummy = ll[i][7].split('(')		# col7 es fin del MC
		ndummy = len(dummy)
		if ndummy==1:
			hrs_end = int(ll[i][7])
		else:
			hrs_end	= int(ll[i][7].split('(')[0][1:])
		tini_mc += [ tini_icme[i-1] + timedelta(hours=hrs_ini) ]
		tend_mc += [ tend_icme[i-1] + timedelta(hours=hrs_end) ]
	# calidad de ICME boundaries
	Qicme 	+= [ ll[i][10] ]		# quality of ICME boundaries
	# flag de MC
	MCsig	+= [ ll[i][15] ]
	#if ll[i][15]=='2H':
	#	MCsig   += [ 2 ]
	#else:
	#	MCsig	+= [ int(ll[i][15]) ]	# MC flag
	#
	Dst	+= [ int(ll[i][16]) ]		# Dst

#--------------------------------------
MCsig 	= array(MCsig)
Dst	= array(Dst)
n_icmes = len(tshck)
#
"""
col0 : id
col1 : disturbance time
col2 : ICME start
col3 : ICME end
col4 : Composition start
col5 : Composition end
col6 : MC start
col7 : MC end
col8 : BDE
col9 : BIF
col10: Quality of ICME boundaries (1=best)
col11: dV --> 'S' indica q incluye shock
col12: V_ICME
col13: V_max
col14: B
col15: MC flag --> '0', '1', '2', '2H': irregular, B-rotation, MC, or MC of "Huttunen etal05" respectively.
col16: Dst
col17: V_transit
col18: LASCO_CME --> time of associated event, generally the CME observed by SOHO/LASCO.
       A veces tiene 'H' por Halo. 
"""
