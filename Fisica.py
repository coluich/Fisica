# v 1.0.1
import math
import time

def helpp():
	print("""
Fg = Forza di gravità
Ka = Coefficiente di attrito
Ke = Coefficiente elastico
dL= Coef delta (differenza di)
V = Velocità
v = volume
D = distanza
d = densità
L = Lunghezza
T = Temperatura
t = tempo
A = Accelerazione
Te = Temperatura di equilibrio
Mo = momento
KDl = coefficente dilatazione
Sup = superficie
fr = frequenza
P = pressione
Pl = pressione liquido
T₀ -> T2
v₀ -> v2
p₀ -> p2
Mmol = Massa
Per altre informazioni contattateci.
""")

Fg=Fper=Fpar=F = E=V=S=D=t=A = Te = Mo=Braccio = Ka=Ke=KDl = dlT=dlV=dlv=dlL=dlE = d=v=v2 = Base=Lato=h=Sup =r=Mmol=N=R=fr=I=W=w=T=T2=M=M2=Cs=N=Cs2=dl=lamda=L=P=Pl=Ohm=Amper=Volt=None
π = 3.1415926535 # pi greco 
g = 9.81 #coefficiente gravitazionale 
G = 6.67*10**-11 #coeff. gravitazionale
C = 300000000   # (m/s) velocità della luce
r = 8,31

list={
"Fg":"N",
"Fper":"N",
"Fpar":"N",
"F":"N",
"E":"J",
"V":"m/s",
"S":"m",
"D":"m",
"t":"s",
"A":"m/s²",
"Te":"°C",
"Mo":"N/m",
"Braccio":"m",
"Ka":"",
"Ke":"N/m",
"dlT":"°C",
"dlV":"m/s",
"dlv":"m³",
"dlL":"m",
"dlE":"J",
"d":"Kg/m³",
"v":"m³",
"v2":"m³",
"Base":"m",
"Lato":"m",
"h":"m",
"Sup":"m²",
"R":"m",
"fr":"hz",
"I":"kg*m²",
"w":"rad/s",
"T":"°C",
"T2":"°C",
"M":"Kg",
"M2":"Kg",
"Cs":"J/Kg°C",
"Cs2":"J/Kg°C",
"dl":"1/°C",
"lamda":"(1/°C) (W/m°C)",
"L":"m",
"P":"Pa",
"Pl":"Pa",
"Ohm":"Ω",
"Volt":"V",
"Amper":"A",
"W":"W",
"N":"mol",
"r":"kPa*L/mol*k",
"Mmol":"g/mol"
}
formulas = [
"E=F*t",
"F=E/t",
"t=E/F",

"E=W*t",
"W=E/t",
"t=E/W",

"V=S/t",
"S=V*t",
"t=S/V",

"A = dlV / t",
"t = sqrt(2*S/A)",
"S = A/2 * t**2",

"v = R**2*π",
"R = sqrt(v/π)",

"d = M/v",
"M = d*v",
"v = M/d",

"Pl = d* g * h",
"d = Pl/ (g * h)",
"g = Pl/ (d * h)",
"h = Pl/ (d * g)",

"P = F / S",
"F = P * S",
"S = F / P",

"E = M * g * h",
"M = E / (g * h)",
"h = E / (M * g)",
"g = E / (M * h)",

"Mo = F * Braccio",
"F = Mo / Braccio",
"Braccio = Mo / F",

"F = Ke * dlL",
"Ke = F / dlL",
"dlL = F / Ke",

"F = Ka * Fper",
"Ka = F / Fper",
"Fper = F / Ka",

"M = F / g",
"F = M * g",
"g = F / M",

"Fg = G * (M * M2)/D**2",
"D = sqrt(F/(G * (M * M2)))",

"F = Ka * Fp",
"Ka = F/Fp",
"Fp = F/Ka",

"Fper = F * (Base/Lato)",
"F = Fper / Base/Lato",
"Base = Fper * Lato / F",
"Lato = Fper / (F * Base)",

"Fpar = F * h/Lato",
"F = Fpar / (h/Lato)",
"h = Fpar * Lato / F",
"Lato = Fpar / (F * h)",

"E = M/2 * **2",
"M = 2*E / v**2",
"V = sqrt(2*E / M)",

"E = Ke/2 * dlL**2",
"dlL = sqrt(2*E / Ke)",
"Ke = 2*E / dlL**2",

"V = 2 * π * R",
"R = V / 2 * π",

"E = M*2* ((π*R*fr) **2)",
"M = 2*E / (π*R*fr) **2",
"R = sqrt(E / (M*2 * π**2 * fr**2))",
"fr = sqrt(E / (M * 2 * π**2 * R**2))",

"E = I/2 * w**2",
"I = 2*E / w**2",
"w = sqrt(2*E / I)",

"I = M * 2 * R",
"R = I / 2*M",
"M = I / 2*R",

"E = π**2 * M * R**2 * fr**2"

"Te = (M*T + M2*T2) / (M + T2)",
"Te = ((Cs*M*T)+(Cs2*M2*T2))/((Cs*M)+(Cs2*M2))",
"Cs2 = ((M+M2)*Cs*(Te-T))/(M2*(T2-Te))",
"Cs = ((M+M2)*Cs2*(Te-T))/(M2*(T2-Te))",

"L = dlL / (lamda * dlT)",
"dlT = dlL / (lamda * L)",
"dlL = L * lamda * dlT",
"lamda = dlL * (dlT  * L)"

"v = dlv / (3 * lamda * dlT)",
"dlv =  v * 3 * lamda * dlT",
"lamda = dlv / (3 * dlT * v)",
"dlT = dlv / (3 * lamda * v)",

"v = dlv / (dl * dlT)",
"dlT = dlv / (dl * v)",
"dl =  dlv / (dlT * v)",
"dlv = dl * dlT * v",

"dlE = (S * lamda * dlT) / D",
"S = (dlE * D) / (lamda * dlT)",
"lamda = (dlE * D) / (S * dlT)",
"dlT = (dlE * D) / (S * lamda)",
"D = dlE / (S * lamda * dlT)",

"T2 = (T1 * v2) / v",
"T1 = (v * T2) / v2",
"v = (T1 * v2) / T2",
"v2 = (v * T2) / T1",

"P = (Sup * lamda * dlT) / D",
"D = (Sup * lamda * dlT) / P",
"lamda = P / (Sup * dlT)",
"Sup = P / (lamda * dlT)",
"dlT = P / (Sup * lamda)",

"v = (T * v2) / T2",
"T = (v * T2) / v2",
"v2 = (v * T2) / T",
"T2 = (T * v2) / v",

"p = (T * p2) / T2",
"T = (p * T2) / p2",
"p2 = (p * T2) / T",
"T2 = (T * p2) / p",

"p = (p2 * v2) / v",
"p2 = (p * v) / v2",
"v = (p2 * v2) / p",
"v2 = (p * v) / p2",


"Amper = Volt / Ohm",
"Volt = Amper /Ohm",
"Ohm = V / Amper",

"W = Volt * Amper",
"Volt = W / Amper",
"Amper = W / Volt",

"W = Volt ** 2 / Ohm",
"Volt = sqrt(W * Ohm)",
"Ohm = Volt ** 2 / W",


"W = Amper ** 2 * Ohm",
"Amper = sqrt(W / Ohm)",
"Ohm = W / Amper ** 2",

"P = r*N*T/v",
"v = r*N*T/P",
"N = P*v/r*T",
"T = P*v/N*r",

"N = M/Mmol",
"M = N*Mmol",
"Mmol = M/N"
]
print ("\n"*100)

while True:
	cont=None
	variabile=None
	while True:
		print("Sviluppato da Siamicas e coluich")
		for i in list.keys():   # stampa list{} con le rispettive chiavi
			print(f"""{i} = {eval(i)} {list.get(i)}""" if eval(i) != None else f"""{i} = """)   
		print("\nPer calcolare digitare q\n")
		print("Per aiuto digitare help\n")
		variabile = input("Cosa vuoi impostare?\n--> ")
		if variabile == "help":
			helpp()
			cont = input("Per continuare digitare un qualsiasi tasto")
			if cont != None:
				print("\n"*100)
				continue
				
		cont=0
		if variabile == "q":    # quit
			print("\n"*2)
			break

		valore = input(f"{variabile} [{list.get(variabile)}] = ")  # raccoglie il valore      
		try:
			exec(f"\ntry:\n    {variabile}=float(valore)\nexcept:\n    pass")  # assegna i valorei raccolti
		except:
			pass
		print("\n"*100)
	
	# variabile  --->   valore  --->   variabile=valore
	
	operation = 0
	print("\u001b[92mOperazioni eseguite:\n")
	for f in formulas:  # elabora le formule
		if '=' in f:		# controlla se il valore da assegnare è diverso da None
			var, formula = f.split('=', 1)
			var = var.strip()
			formula = formula.strip()
			if var not in locals() or locals()[var] is None:
				try:
					exec(f)
					print(f"""\u001b[92m{f}""")
					operation=operation+1
				except:
					pass
	time.sleep(1)
	if operation == 0:
		print("\u001b[31m0")

	print("\n\033[34mDati:\n")
	numero= 0
	for i in list.keys():   # stampa list{} con i rispettivi valori
		if eval(i) is not None:
			print(f"""\033[34m{i} = {eval(i)} {list.get(i)}\033[0m""")
		numero = numero + 1
	if numero == 0:
		print("\n\u001b[31mNessun dato")
	if operation == 0 and numero != 0:
		print("\n\u001b[31mNessun operazione possibile con i dati forniti.")
	Sno = 0
	while True:
		Sno= input("\n\n\033[0mVuoi resettare i dati? [s/N]:")
		if Sno == "S" or Sno == "s":
			Fg=Fper=Fpar=F = E=V=S=D=t=A = Te = Mo=Braccio = Ka=Ke=KDl = dlT=dlV=dlv=dlL=dlE = d=v=v2 = Base=Lato=h=Sup = r=N=Mmol=R=fr=I=W=w=T=T2=M=M2=Cs=Cs2=dl=lamda=L=P=Pl=Ohm=Amper=Volt=None
			break
		elif Sno == "N" or Sno == "n" or Sno == "" or None:
			print("\n"*100)
			break
		else:
			print("Il carattere specificato non rispetta le richieste.")
			
		variabile = 0
