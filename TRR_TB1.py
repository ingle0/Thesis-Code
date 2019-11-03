#Two-stage Receptor Recruitment (TRR) model for TGF-beta1

#Importing functions for model creation
from __future__ import print_function
from pysb import *
_pysb_doctest_suppress_modelexistswarning = True

#Create model
Model()

Monomer('BG')
Monomer('TBRII')
Monomer('TBRI')
Monomer('TB1_BG')
Monomer('TB1_BG_TBRII')
Monomer('TB1_BG_TBRII_TBRI')
Monomer('TB1_TBRI_TBRI')
Monomer('TB1_TBRII_TBRI')
Monomer('TB1_TBRII_TBRI_TBRI')
Monomer('TB1_TBRII_TBRII_TBRI')
Monomer('TB1_TBRII_TBRII_TBRI_TBRI')
Monomer('TB1_TBRII')
Monomer('TB1_TBRI')
Monomer('TB1_TBRII_TBRII')

#List reactions
#Reaction rates for this system can be found in Supplemental Table 7

#Rxn 13a
Parameter('k13a_on', MF * Z)
Parameter('k13a_off', MR)
#BG + TB1 <-> TB1_BG
Rule('TB1_to_BG', BG() <> TB1_BG(), k13a_on, k13a_off)

#Rxn 14a
Parameter('k14a_on', NF*SEF)
Parameter('k14a_off', NR)
#TB1_BG + TBRII <-> TB1_BG_TBRII
Rule('TB1_BG_to_TB1_BG_TBRII', TB1_BG() + TBRII() <> TB1_BG_TBRII(), k14a_on, k14a_off)

#Rxn 16
Parameter('k16a_on', QF*SEF)
Parameter('k16a_off', QR)
#TB1_BG_TBRII + TBRI <-> TB1_BG_TBRII_TBRI
Rule('TB1_BG_TBRII_to_TB1_BG_TBRII_TBRI', TB1_BG_TBRII() + TBRI() <> TB1_BG_TBRII_TBRI(), k16a_on, k16a_off)

#Rxn 17a # inverted reaction
Parameter('k17a_on', SF*SEF)
Parameter('k17a_off', SR)
#TB1_BG_TBRII_TBRI <-> TB1_TBRII_TBRI + BG
Rule('TB1_BG_TBRII_TBRI_to_TB1_TBRII_TBRI', TB1_TBRII_TBRI() + BG() <> TB1_BG_TBRII_TBRI(), k17a_on, k17a_off)

#Rxn 8a
Parameter('k8a_on', HF*SEF)
Parameter('k8a_off', HR)
#TB1_TBRII_TBRI + TBRII <-> TB1_TBRII_TBRII_TBRI
Rule('TB1_TBRII_TBRI_to_TB1_TBRII_TBRII_TBRI', TB1_TBRII_TBRI() + TBRII() <> TB1_TBRII_TBRII_TBRI(), k8a_on, k8a_off)

#Rxn 9a
Parameter('k9a_on', IF*SEF)
Parameter('k9a_off', IR)
#TB1_TBRII_TBRI + TBRI <-> TB1_TBRII_TBRI_TBRI
Rule('TB1_TBRII_TBRI_to_TB1_TBRII_TBRI_TBRI', TB1_TBRII_TBRI() + TBRI() <> TB1_TBRII_TBRI_TBRI(), k9a_on, k9a_off)

#Rxn 11a
Parameter('k11a_on', KF*SEF)
Parameter('k11a_off', KR)
#TB1_TBRII_TBRII_TBRI + TBRI <-> TB1_TBRII_TBRII_TBRI_TBRI
Rule('TB1_TBRII_TBRII_TBRI_to_TB1_TBRII_TBRII_TBRI_TBRI', TB1_TBRII_TBRII_TBRI() +TBRI() <> TB1_TBRII_TBRII_TBRI_TBRI(), k11a_on, k11a_off)

#Rxn 12a
Parameter('k12a_on', LF*SEF)
Parameter('k12a_off', LR)
#TB1_TBRII_TBRI_TBRI + TBRII <-> TB1_TBRII_TBRII_TBRI_TBRI
Rule('TB1_TBRII_TBRI_TBRI_to_TB_TBRII_TBRII_TBRI_TBRI', TB1_TBRII_TBRI_TBRI() + TBRII() <> TB1_TBRII_TBRII_TBRI_TBRI(), k12a_on, k12a_off)

#Rxn 1
Parameter('k1a_on', AF*Z)
Parameter('k1a_off', AR)
#TB1 + TBRII <-> TB1_TBRII
Rule('TB1_to_TB1_TBRII', TBRII() <> TB1_TBRII(), k1a_on, k1a_off)

#Rxn 4
Parameter('k4a_on', DF*SEF)
Parameter('k4a_off', DR)
#TB1_TBRII + TBRI <-> TB1_TBRII_TBRI
Rule('TB1_TBRII_to_TB1_TBRII_TBRI', TB1_TBRII() + TBRI() <> TB1_TBRII_TBRI(), k4a_on, k4a_off)

#Rxn 3
Parameter('k3a_on', CF*SEF)
Parameter('k3a_off', CR)
#TB1_TBRII + TBRII <-> TB1_TBRII_TBRII
Rule('TB1_TBRII_to_TB1_TBRII_TBRII', TB1_TBRII() +TBRII() <> TB1_TBRII_TBRII(), k3a_on, k3a_off)

#Rxn 7
Parameter('k7a_on', GF*SEF)
Parameter('k7a_off', GR)
#TB1_TBRII_TBRII +TBRI <> TB1_TBRII_TBRII_TBRI
Rule('TB1_TBRII_TBRII_to_TB1_TBRII_TBRII', TB1_TBRII_TBRII() +TBRI() <> TB1_TBRII_TBRII_TBRI(), k7a_on, k7a_off)

#Rxn 2
Parameter('k2a_on', BF*Z)
Parameter('k2a_off', BR)
#TB1 + TBRI <-> TB1_TBRI
Rule('TB1_to_TB1_TBRI', TBRI() <> TB1_TBRI(), k2a_on, k2a_off)

#Rxn 5
Parameter('k5a_on', EF*SEF)
Parameter('k5a_off', ER)
#TB1_TBRI + TBRII <-> TB1_TBRI_TBRI
Rule('TB1_TBRI_to_TB1_TBRII_TBRI', TB1_TBRI() + TBRII() <> TB1_TBRII_TBRI(), k5a_on, k5a_off)

#Rxn 6
Parameter('k6a_on', FF*SEF)
Parameter('k6a_off', FR)
#TB1_TBRI + TBRI <-> TB1_TBRI_TBRI
Rule('TB1_TBRI_to_TB1_TBRI_TBRI', TB1_TBRI() + TBRI() <> TB1_TBRI_TBRI(), k6a_on, k6a_off)

#Rxn 10
Parameter('k10a_on', JF*SEF)
Parameter('k10a_off', JR)
#TB1_TBRI_TBRI + TBRII <-> TB1_TBRII_TBRI_TBRI
Rule('TB1_TBRI_TBRI_to_TB1_TBRII_TBRI_TBRI', TB1_TBRI_TBRI() + TBRII() <> TB1_TBRII_TBRI_TBRI(), k10a_on, k10a_off)

#Rxn 15
Parameter('k15a_on', PF*SEF)
Parameter('k15a_off', PR)
#TB1_TBRI + TBRI <-> TB1_TBRI_TBRI
Rule('TB1_TBRII_to_TB1_BG_TBRII', TB1_TBRII() + BG() <> TB1_BG_TBRII(), k15a_on, k15a_off)

#Receptor Recycling and Endocytosis in one step

#Rxn1001

Parameter('kendo', kendo)
#TB1_BG -> BG + TB1
Rule('TB1_BG_Endo1', TB1_BG() >> BG(), kendo)

#Rxn1002

#TB1_BG_TBRII -> BG +TBRII + TB1
Rule('TB1_BG_TBRII_Endo2', TB1_BG_TBRII() >> BG() + TBRII(), kendo)

#Rxn1003

#TB1_TBRII -> TB1 + TBRII
Rule('TB1_TBRII_Endo3', TB1_TBRII() >> TBRII(), kendo)

#Rxn1004

#TB1_TBRI -> TBRI + TB1
Rule('TB1_TBRI_Endo4', TB1_TBRI() >> TBRI(), kendo)

#Rxn1005
#TB1_BG_TBRII_TBRI ->
Rule('TB1_BG_TBRII_TBRI_Endo5', TB1_BG_TBRII_TBRI() >> BG() + TBRII() + TBRI(), kendo)

#Rxn1006

#TB1_TBRII_TBRI -> TBRII + TBRI + TB1
Rule('TB1_TBRII_TBRI_Endo6', TB1_TBRII_TBRI() >> TBRII() + TBRI(), kendo)

#Rxn1007

#TB1_TBRII_TBRII_TBRI -> TBRII(2) + TBRI + TB1
Rule('TB1_TBRII_TBRII_TBRI_Endo7', TB1_TBRII_TBRII_TBRI() >> TBRII() + TBRII() + TBRI(), kendo)

#Rxn1008

#TB1_TBRII_TBRI_TBRI -> TBRII + TBRI + TBRI + TB1
Rule('TB1_TBRII_TBRI_TBRI_Endo8', TB1_TBRII_TBRI_TBRI() >> TBRII() + TBRI() + TBRI(), kendo)

#Rxn1009

#TB1_TBRII_TBRII_TBRI_TBRI -> TBRII + TBRII + TBRI + TBRI + TB1
Rule('TB1_TBRII_TBRII_TBRI_TBRI_Endo9', TB1_TBRII_TBRII_TBRI_TBRI() >> TBRII() + TBRII() + TBRI() + TBRI(), kendo)

#Rxn1010

#TB1_TBRII_TBRII -> TRBII + TBRII + TB1
Rule('TB1_TBRII_TBRII_Endo10', TB1_TBRII_TBRII() >> TBRII() + TBRII(), kendo)

#Rxn1011
Rule('TB1_TBRI_TBRI_Endo11', TB1_TBRI_TBRI() >> TBRI() + TBRI(), kendo)

#Defining parameters that are not rates
Parameter('BG_0', Z2)
Parameter('TBRII_0', RII)
Parameter('TBRI_0', RI)

#Initial Values
#These values would vary depending on screen performed
Initial(BG(), BG_0)
Initial(TBRII(), TBRII_0)
Initial(TBRI(), TBRI_0)

#Observable Amount of Monomer
Observable('BG', BG());
Observable('TBRII', TBRII());
Observable('TBRI', TBRI());
Observable('TB1_TBRI_TBRI', TB1_TBRI_TBRI());
Observable('TB1_BG', TB1_BG());
Observable('TB1_TBRI', TB1_TBRI());
Observable('TB1_TBRII', TB1_TBRII());
Observable('TB1_TBRII_TBRII', TB1_TBRII_TBRII());
Observable('TB1_BG_TBRII', TB1_BG_TBRII());
Observable('TB1_BG_TBRII_TBRI', TB1_BG_TBRII_TBRI());
Observable('TB1_TBRII_TBRI', TB1_TBRII_TBRI());
Observable('TB1_TBRII_TBRII_TBRI', TB1_TBRII_TBRII_TBRI());
Observable('TB1_TBRII_TBRI_TBRI', TB1_TBRII_TBRI_TBRI());
Observable('TB1_TBRII_TBRII_TBRI_TBRI', TB1_TBRII_TBRII_TBRI_TBRI());

#Importing functions for data analysis
from matplotlib.pyplot import *
from numpy import linspace, array
from pysb.simulator import ScipyOdeSimulator

# We will integrate from t=0 to t=86400
t = linspace(0, 186400, 186401)
y = ScipyOdeSimulator(model, rtol=1e-4, atol=[1e-8, 1e-14, 1e-6]).run(tspan=t).all

# Gather the observables of interest into a matrix
yobs = array([y[obs] for obs in ('BG', 'TBRII', 'TBRI', 'TB1_BG', 'TB1_TBRI', 'TB1_TBRII', 'TB1_TBRII_TBRII', 'TB1_BG_TBRII', 'TB1_BG_TBRII_TBRI', 'TB1_TBRII_TBRI', 'TB1_TBRII_TBRII_TBRI', 'TB1_TBRII_TBRI_TBRI', 'TB1_TBRII_TBRII_TBRI_TBRI', 'TB1_TBRI_TBRI')]).T

#Recording the concentrations of individual species at 186400 seconds
TetramerShot = ((y['TB1_TBRII_TBRII_TBRI_TBRI'][186400]))
Shot_BG = ((y['BG'][186400]))
ShotTBRII = ((y['TBRII'][186400]))
ShotTBRI = ((y['TBRI'][186400]))
ShotTB1_BG = ((y['TB1_BG'][186400]))
ShotTB1_TBRI = ((y['TB1_TBRI'][186400]))
ShotTB1_TBRII = ((y['TB1_TBRII'][186400]))
ShotTB1_TBRI_TBRI = ((y['TB1_TBRI_TBRI'][186400]))
ShotTB1_TBRII_TBRII = ((y['TB1_TBRII_TBRII'][186400]))
ShotTB1_BG_TBRII = ((y['TB1_BG_TBRII'][186400]))
ShotTB1_BG_TBRII_TBRI = ((y['TB1_BG_TBRII_TBRI'][186400]))
ShotTB1_TBRII_TBRI = ((y['TB1_TBRII_TBRI'][186400]))
ShotTB1_TBRII_TBRII_TBRI = ((y['TB1_TBRII_TBRII_TBRI'][186400]))
ShotTB1_TBRII_TBRI_TBRI = ((y['TB1_TBRII_TBRI_TBRI'][186400]))
ShotTB1_TBRII_TBRII_TBRI_TBRI = ((y['TB1_TBRII_TBRII_TBRI_TBRI'][186400]))

#End of TRR model for TGF-beta1
