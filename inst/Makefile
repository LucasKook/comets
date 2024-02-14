dependencies:
	(cd .. && Rscript --vanilla inst/dependencies.R > dependencies.out)

mimic:
	(cd .. && Rscript --vanilla inst/application/mimic.R > mimic.out)

multiomics:
	(cd .. && Rscript --vanilla inst/application/multiomics.R > multiomics.out)

ccle:
	(cd .. && Rscript --vanilla inst/application/ccle.R > ccle.out)

dep:
	bash dispatch.sh run-install.sh 1200

mim:
	bash dispatch.sh run-mimic.sh 21600

mul:
	bash dispatch.sh run-multiomics.sh 3600

sig:
  bash dispatch.sh run-ccle.sh 3600
