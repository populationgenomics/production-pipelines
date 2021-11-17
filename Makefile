default: patch package

.PHONY: patch
patch:
	bump2version patch
	git push

.PHONY: minor
minor:
	bump2version minor
	git push

.PHONY: package
package:
	rm -rf dist/*
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: sleep
sleep:
	sleep 60

.PHONY: seqr
seqr:
	python pipelines/seqr_loader.py \
	-n main \
	--analysis-project seqr \
	--input-project acute-care \
	--input-project perth-neuro \
	--output-project acute-care \
	--ped-file gs://cpg-acute-care-main-upload/cpg_acute_positives_20211003_213917/acute-care-topup-mfranklin.ped \
	--ped-file gs://cpg-acute-care-main-upload/acute-care-sm.ped \
	--ped-file gs://cpg-perth-neuro-main-upload/perth-neuro-sm.ped \
	--last-stage GvcfPedCheck \
	--validate-smdb-analyses \
	--check-smdb-seq-existence \
	-S CPG11783 \
	-S CPG13326 \
	--keep-scratch

.PHONY: pedigree
pedigree:
	python pipelines/pedigree.py \
	-n main \
	--analysis-project seqr \
	--input-project acute-care \
	--input-project perth-neuro \
	--ped-file gs://cpg-acute-care-main-upload/acute-care-sm.ped \
	--ped-file gs://cpg-perth-neuro-main-upload/perth-neuro-sm.ped \
	-S CPG11783 \
	-S CPG13326 \
	--keep-scratch
