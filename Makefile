

generate:
	python3 src/genAlpha.py

train:
	matlab -nodisplay -nosplash -nodesktop -r "run('src/trainGP.m');exit;"

test:
	`(source /usr/local/bin/setup_cantera)`
	python3 src/test_PR.py