import subprocess
import sys

if __name__ == "__main__":
	#nuclei = ["all"]
	#nuclei = ["Co","Ni","Cu","Zn","Ga","Ge"]
	nuclei = ["Ga"]
	for nucl in nuclei:
		cmd = "anamerger_main -c config/anamerger/config_"+nucl+".yaml -o anamerger_outputs/anamerger_output_"+nucl+"_.root >> logs/anamerger_"+nucl+".log 2>> logs/anamerger_"+nucl+".errlog"
		sp = subprocess.call([cmd], shell=True)
	
