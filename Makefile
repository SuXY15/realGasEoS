

generate:
	python3 src/genAlpha.py

train:
	matlab -nodisplay -nosplash -nodesktop -r "run('src/trainGP.m');exit;"